function result = minimalReactome(model, varargin)
% Finds the smallest set of reactions that still meets a requirement
%
% USAGE:
%
%    result = minimalReactome(model, varargin)
%
% INPUTS:
%    model:             COBRA model structure to cut down. It is not
%                       modified.
%
% OPTIONAL INPUTS (name/value pairs):
%    media:             Struct array, one entry per medium, with fields:
%
%                         * rxns - the exchange reactions to open
%                         * lb - the lower bound to open each to
%                         * label - a name for the medium (optional)
%
%                       Everything not named is closed, and that matters: a
%                       medium is defined as much by what it withholds as by
%                       what it supplies. Pass several to require that the
%                       one network works in all of them: each becomes its
%                       own condition with its own flux vector, so the
%                       result is the smallest network that grows
%                       aerobically *and* anaerobically rather than the
%                       smaller of the two. (Default: none, meaning leave
%                       the model's own bounds alone)
%    biomassReaction:   The reaction whose flux is growth. (Default: the
%                       model's objective, if there is exactly one)
%    minGrowth:         The growth rate to demand, as an absolute lower
%                       bound on the biomass reaction. Mutually exclusive
%                       with growthFraction.
%    growthFraction:    The growth rate to demand, as a fraction of what the
%                       model achieves in each medium before anything is
%                       removed. Pass 0 to not require growth. (Default: 1,
%                       i.e. lose no growth at all)
%    products:          Struct array with fields id and rate, naming
%                       metabolites or reactions the result must still make
%                       and the flux to demand of each as a lower bound. A
%                       reaction is bounded directly; a metabolite is
%                       resolved to the boundary reaction the model already
%                       has for it. A metabolite the model does not exchange
%                       is refused rather than given one - see the note.
%                       These are required alongside growth, in the same
%                       flux distribution. (Default: none)
%    tasks:             Struct array of metabolic tasks the result must
%                       still be able to perform, as metabolicTask or
%                       parseTaskList returns them. Each becomes its own
%                       condition, with its own flux vector and its own
%                       medium, sharing the one set of binaries.
%                       (Default: none)
%    keepRxns:          Cell array of reactions to retain whatever they do.
%                       They are given no binary, so they cost nothing and
%                       cannot be removed. MinReact's "reactions to be
%                       retained". (Default: none)
%    weights:           What each reaction costs to keep, as a vector in
%                       model order. Non-negative. A cost of 0 still leaves
%                       the reaction removable - use keepRxns to force one
%                       in. (Default: ones)
%    tol:               The flux threshold. (Default: 1e-4)
%    preprocess:        Whether to decide what can be decided before the
%                       MILP starts. Two passes, both exact: reactions that
%                       can carry no flux in *any* condition are deleted
%                       outright, and reactions some condition is infeasible
%                       without are proved essential and given no binary.
%                       (Default: true)
%    nSol:              How many alternative minimal reactomes to return.
%                       Each is forbidden from recurring exactly before the
%                       next solve, so they are distinct but not necessarily
%                       all of them, and only the first is guaranteed
%                       minimal. (Default: 1)
%    solveTime:         Seconds to allow the solve. A solution found before
%                       the limit is used even if minimality was not proven.
%                       (Default: 7200)
%    inclusionCutoff:   The absolute flux at which a reaction is kept
%                       whatever its binary says. (Default: tol*1e-7)
%    seed:              Seed for the consistency check's random
%                       coefficients, making preprocessing reproducible.
%                       (Default: none)
%    check:             Whether to verify the result against every
%                       requirement and report what it misses.
%                       (Default: true)
%
% OUTPUTS:
%    result:            Struct, or struct array when nSol > 1, with fields:
%
%                         * model - the minimal reactome, as a COBRA model
%                         * kept - the reactions it retained
%                         * removed - the reactions dropped, including the
%                           blocked ones preprocessing deleted
%                         * growth - the rate reached in each medium
%                         * requiredGrowth - the rate it had to reach
%                         * production - the flux each required product
%                           carries in each medium
%                         * requiredProduction - the rate each had to carry
%                         * conditions - the labels imposed, in order
%                         * essential - the reactions preprocessing proved
%                           are in every possible answer
%                         * keptButBlocked - per growth condition, the
%                           keepRxns that cannot carry tol of flux in it
%                         * unsatisfied - the conditions the result does not
%                           meet. Empty is the expected outcome
%                         * summary - a printable account of the above
%
% EXAMPLE:
%
%    ::
%
%        model = readCbModel('ecoli_core_model.mat');
%        M9 = {'EX_pi(e)','EX_h2o(e)','EX_h(e)','EX_nh4(e)','EX_o2(e)','EX_co2(e)'};
%        medium = struct('label','glucose', ...
%                        'rxns',{[M9,{'EX_glc(e)'}]}, ...
%                        'lb',{[-1000*ones(1,numel(M9)), -10]});
%        result = minimalReactome(model, 'media', medium, 'growthFraction', 0.5);
%        disp(result.summary)
%
% NOTE:
%    A genome-scale reconstruction is hugely redundant. Most of its
%    reactions are not needed to grow on any one carbon source, and a large
%    subnetwork can be deleted without the model noticing. The minimal
%    reactome is the smallest set that still meets a stated requirement -
%    the original question of Burgard, Vaidyaraman and Maranas (2001), who
%    asked how few reactions Escherichia coli needs to grow and found the
%    answer depends entirely on what is in the medium: 224 reactions on
%    glucose alone, 229 on acetate alone, 122 on a medium engineered to
%    supply as much as possible. A minimal reactome is minimal *for an
%    environment*, never in the abstract.
%
%    That is a minimum-cardinality problem, so it is a MILP. Every reaction
%    gets a binary tied to its flux by lb*y <= v <= ub*y, so that y = 0 pins
%    the flux to zero, and the weighted count is minimised subject to
%    S.v = 0 and whatever the network has to achieve. That is minNetMILP
%    with no core reactions at all: the requirements enter as bounds.
%
%    Growth, production and the medium are one condition, solved as one flux
%    vector. Tasks are not, because a task defines its own medium: the
%    conditions are therefore replicated - one copy of the network each -
%    and tied together by one binary per reaction, so what is minimised is
%    the size of the union while every condition is satisfied separately.
%    The model handed to the solver is k times the size of yours, which is
%    what makes a long task list expensive. Solving them jointly is what
%    makes the answer minimal *across* them; filling them one at a time and
%    taking the union, as gapfillForTasks does, gives a feasible network but
%    not the smallest one.
%
%    **Nothing is added to the network.** Requiring the production of a
%    metabolite the model does not exchange would need a drain for it, and a
%    drain added here would make the result the smallest subnetwork of a
%    model you do not have, feasible only because of the reaction that was
%    added to measure it. Add it to your own model first, where it is
%    visible and counted like any other.
%
%    This is the opposite operation to gapfillForTasks and they use the same
%    vocabulary of tasks. Gap-filling adds reactions until the tasks pass;
%    this removes reactions while they still do.
%
% REFERENCES:
%    Burgard, A. P., Vaidyaraman, S., and Maranas, C. D. (2001). Minimal
%    reaction sets for Escherichia coli metabolism under different growth
%    requirements and uptake environments. Biotechnology Progress, 17(5),
%    791-797.
%
%    Sambamoorthy, G., and Raman, K. (2020). MinReact: a systematic approach
%    for identifying minimal metabolic networks. Bioinformatics, 36(15),
%    4309-4315.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: spectraME, gapfillForTasks, minNet

parser = inputParser();
parser.addParameter('media', [], @(x) isstruct(x) || isempty(x));
parser.addParameter('biomassReaction', '', @ischar);
parser.addParameter('minGrowth', [], @(x) isnumeric(x) && isscalar(x) || isempty(x));
parser.addParameter('growthFraction', [], @(x) isnumeric(x) && isscalar(x) || isempty(x));
parser.addParameter('products', [], @(x) isstruct(x) || isempty(x));
parser.addParameter('tasks', [], @(x) isstruct(x) || isempty(x));
parser.addParameter('keepRxns', {}, @(x) iscell(x) || ischar(x));
parser.addParameter('weights', [], @isnumeric);
parser.addParameter('tol', 1e-4, @isnumeric);
parser.addParameter('preprocess', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('nSol', 1, @isnumeric);
parser.addParameter('solveTime', 7200, @isnumeric);
parser.addParameter('inclusionCutoff', [], @isnumeric);
parser.addParameter('seed', [], @isnumeric);
parser.addParameter('check', true, @(x) islogical(x) || isnumeric(x));
parser.parse(varargin{:});
opt = parser.Results;

if ~isempty(opt.minGrowth) && ~isempty(opt.growthFraction)
    error('minimalReactome:bothGrowthRequirements', ...
        ['Pass minGrowth or growthFraction, not both: one is an absolute ' ...
         'rate and the other a fraction of what the model achieves, so ' ...
         'giving both leaves it ambiguous which applies.']);
end
if opt.nSol < 1
    error('minimalReactome:badSolutionCount', ...
        'nSol must be at least 1, not %d.', opt.nSol);
end

tol = opt.tol;
inclusionCutoff = opt.inclusionCutoff;
if isempty(inclusionCutoff)
    inclusionCutoff = tol * 1e-7;
end
keepRxns = opt.keepRxns;
if ischar(keepRxns)
    keepRxns = {keepRxns};
end
keepRxns = keepRxns(:);

base = model;
originalRxns = base.rxns(:);
nOriginal = numel(originalRxns);

weights = opt.weights;
if isempty(weights)
    weights = ones(nOriginal, 1);
end
if numel(weights) ~= nOriginal
    error('minimalReactome:badWeightCount', ...
        'weights must hold one entry per reaction: got %d for %d reactions.', ...
        numel(weights), nOriginal);
end
if any(weights < 0)
    error('minimalReactome:negativeWeight', ...
        ['A negative weight pays the objective to keep a reaction, which ' ...
         'is not what minimising a count means. Use keepRxns to force one in.']);
end

%% what growth is, and whether any was asked for
growthRequested = (~isempty(opt.minGrowth) && opt.minGrowth > 0) || ...
    (isempty(opt.minGrowth) && (isempty(opt.growthFraction) || opt.growthFraction > 0));
biomass = resolveBiomass(base, opt.biomassReaction, growthRequested);

%% the media, each one condition
media = normaliseMedia(opt.media);
for k = 1:numel(media)
    checkMedium(base, media(k));
end

%% what must be produced, and how fast
products = resolveProducts(base, opt.products);

%% the rate to demand in each medium
requiredGrowth = zeros(numel(media), 1);
for k = 1:numel(media)
    requiredGrowth(k) = requiredGrowthIn(base, media(k), biomass, ...
        opt.minGrowth, opt.growthFraction);
end

unknown = keepRxns(~ismember(keepRxns, originalRxns));
if ~isempty(unknown)
    error('minimalReactome:unknownKeepRxn', ...
        'keepRxns names %d reaction(s) the model does not have, for example ''%s''.', ...
        numel(unknown), unknown{1});
end

tasks = opt.tasks;
if ~any(requiredGrowth > 0) && isempty(products) && isempty(tasks)
    hint = '';
    if isempty(biomass)
        hint = [' This model has no unambiguous objective reaction, so ' ...
                'growth was not required of it; pass biomassReaction to ' ...
                'say which reaction growth is.'];
    end
    error('minimalReactome:nothingRequired', ...
        ['Nothing was asked of the minimal reactome, so the smallest ' ...
         'network meeting the requirement is the empty one. Give a growth ' ...
         'requirement, a product or a task list. A medium on its own is ' ...
         'not a requirement - it says what may be taken up, not what has ' ...
         'to happen.%s'], hint);
end

%% every requirement has to be reachable before anything is removed
checkTasksArePossible(base, tasks);
for k = 1:numel(media)
    checkGrowthIsPossible(base, media(k), biomass, requiredGrowth(k), products);
end

%% the conditions, each a copy of the network under its own bounds
[conditions, keptButBlocked, forceable] = buildConditions(base, media, biomass, ...
    requiredGrowth, products, tasks, keepRxns, tol);

%% what can be decided before the MILP starts
blocked = {};
essential = {};
if opt.preprocess
    [blocked, essential] = preprocessConditions(conditions, originalRxns, ...
        keepRxns, tol, opt.seed);
end

protectedRxns = setdiff(union(keepRxns, essential), blocked);
counted = setdiff(originalRxns, union(protectedRxns, blocked));
% setdiff sorts; the binaries are cheaper to read in model order.
counted = originalRxns(ismember(originalRxns, counted));

if isempty(counted)
    fprintf(['minimalReactome: preprocessing settled every reaction, so ' ...
             'there is nothing left for the MILP to decide.\n']);
end

%% one copy of the network per condition, tied together by one binary each
[joint, groups, copyOf] = replicateConditions(conditions, originalRxns, ...
    blocked, counted);

costs = zeros(numel(counted), 1);
for g = 1:numel(counted)
    costs(g) = weights(strcmp(originalRxns, counted{g}));
end

%% the reactions that must carry flux, named in the joint model
core = [];
for k = 1:numel(conditions)
    usable = forceable{k};
    for i = 1:numel(usable)
        at = copyOf{k}(strcmp(originalRxns, usable{i}));
        if at > 0
            core(end+1, 1) = at;                                  %#ok<AGROW>
        end
    end
end

fprintf(['minimalReactome: minimising %d reactions over %d condition(s) ' ...
         '(%s): %d reactions and %d metabolites in the joint problem, %d ' ...
         'forced.\n'], numel(counted), numel(conditions), ...
    strjoin(cellfun(@(c) c.label, conditions, 'UniformOutput', false), ', '), ...
    numel(joint.rxns), numel(joint.mets), numel(core));

%% solve, and again for each alternative asked for
direction = settleDirections(joint, core, tol);
prevSols = {};
solutions = cell(0, 1);
for s = 1:opt.nSol
    [reacInd, ~, stat, selected] = minNet(joint, direction, costs, tol, 1, ...
        'minNetMILP', opt.solveTime, [], prevSols, groups, inclusionCutoff);
    if stat ~= 1 && stat ~= 3
        if s == 1
            error('minimalReactome:infeasible', ...
                ['The joint problem could not be solved, so no minimal ' ...
                 'reactome was found. Every requirement was reachable on ' ...
                 'its own, so the usual cause is that they cannot all be ' ...
                 'met at once.']);
        end
        fprintf(['minimalReactome: no further alternative solution after ' ...
                 '%d of %d.\n'], s - 1, opt.nSol);
        break
    end
    solutions{end+1, 1} = reacInd;                                %#ok<AGROW>
    prevSols = [prevSols; {find(selected)}];                      %#ok<AGROW>
end

%% map each solution back to the model's own reactions, and report
result = repmat(emptyResult(), numel(solutions), 1);
for s = 1:numel(solutions)
    kept = keptReactions(solutions{s}, joint, originalRxns, copyOf, protectedRxns);
    result(s) = reportReactome(base, originalRxns, kept, media, biomass, ...
        requiredGrowth, products, tasks, conditions, essential, ...
        keptButBlocked, tol, opt.check);
    fprintf('minimalReactome: solution %d: %s\n', s, firstLine(result(s).summary));
end
if isscalar(result)
    return
end
end


function biomass = resolveBiomass(model, given, required)
% Returns the reaction whose flux counts as growth, or '' if there is none.
if ~isempty(given)
    if ~ismember(given, model.rxns)
        error('minimalReactome:noSuchReaction', ...
            'No reaction ''%s'' in the model.', given);
    end
    biomass = given;
    return
end
inObjective = model.rxns(model.c ~= 0);
if isscalar(inObjective)
    biomass = inObjective{1};
    return
end
biomass = '';
if required
    error('minimalReactome:ambiguousObjective', ...
        ['The model has %d reactions in its objective, so which one is ' ...
         'growth is ambiguous. Pass biomassReaction, or growthFraction 0 ' ...
         'if growth is not what you are requiring.'], numel(inObjective));
end
end


function media = normaliseMedia(media)
% Returns the media as a labelled struct array, or one unnamed condition.
if isempty(media)
    % No medium at all is still one condition: the model's own bounds.
    media = struct('label', {'growth'}, 'rxns', {{}}, 'lb', {[]});
    return
end
if ~all(isfield(media, {'rxns', 'lb'}))
    error('minimalReactome:badMedia', ...
        ['media must be a struct array with fields rxns and lb, and ' ...
         'optionally label.']);
end
for k = 1:numel(media)
    if numel(media(k).rxns) ~= numel(media(k).lb)
        error('minimalReactome:mediumLengthMismatch', ...
            'medium %d names %d reactions but gives %d bounds.', ...
            k, numel(media(k).rxns), numel(media(k).lb));
    end
    if ~isfield(media, 'label') || isempty(media(k).label)
        media(k).label = sprintf('medium %d', k);
    end
end
labels = {media.label};
if numel(unique(labels)) ~= numel(labels)
    error('minimalReactome:duplicateMediumLabel', ...
        ['Two media share a label, which would make them one copy of the ' ...
         'network rather than two. Media labels must be distinct.']);
end
end


function checkMedium(model, medium)
% Stops a medium naming an exchange the model does not have.
unknown = medium.rxns(~ismember(medium.rxns, model.rxns));
if ~isempty(unknown)
    error('minimalReactome:unknownExchange', ...
        ['Medium ''%s'' names %d reaction(s) the model does not have, for ' ...
         'example ''%s''.'], medium.label, numel(unknown), unknown{1});
end
end


function constrained = applyMedium(model, medium)
% Closes every exchange, then opens the ones the medium names.
%
% A medium is checked against the model it was given for, but it is applied
% again to the reduced result, which no longer has every exchange. One the
% result dropped is one it does not need, so it is skipped rather than
% treated as an error.
constrained = model;
if isempty(medium.rxns)
    return
end
isExc = findExcRxns(constrained);
constrained.lb(isExc) = 0;
at = findRxnIDs(constrained, medium.rxns);
bounds = medium.lb(:);
constrained.lb(at(at > 0)) = bounds(at > 0);
end


function products = resolveProducts(model, given)
% Turns requested products into reactions with a rate to demand.
%
% A metabolite is resolved to the boundary reaction the model already has
% for it, and nothing is added. Adding a demand reaction for an unexchanged
% metabolite is the obvious convenience and it is wrong: the result would be
% the smallest subnetwork of a model the caller does not have, feasible only
% because of the drain that was added to measure it.
products = struct('rxn', {}, 'rate', {}, 'asked', {});
if isempty(given)
    return
end
if ~all(isfield(given, {'id', 'rate'}))
    error('minimalReactome:badProducts', ...
        'products must be a struct array with fields id and rate.');
end

isExc = findExcRxns(model);
for i = 1:numel(given)
    rate = given(i).rate;
    if rate <= 0
        error('minimalReactome:badProductRate', ...
            ['The production rate demanded of ''%s'' is %g, which requires ' ...
             'nothing: a zero flux satisfies it. Give a positive rate.'], ...
            given(i).id, rate);
    end
    if ismember(given(i).id, model.rxns)
        products(end+1) = struct('rxn', given(i).id, 'rate', rate, ...
            'asked', given(i).id);                                %#ok<AGROW>
        continue
    end
    at = findMetIDs(model, given(i).id);
    if at == 0
        error('minimalReactome:unknownProduct', ...
            ['''%s'' is neither a reaction nor a metabolite of the model, ' ...
             'so there is no way to require its production.'], given(i).id);
    end

    touching = find(model.S(at, :) ~= 0);
    exchanges = touching(isExc(touching));
    if isempty(exchanges)
        error('minimalReactome:unexchangedProduct', ...
            ['''%s'' has no exchange reaction, so its production cannot be ' ...
             'required without adding one - and adding one would change the ' ...
             'network being minimised. The answer would be the smallest ' ...
             'subnetwork of a model you do not have, feasible only because ' ...
             'of the drain added to measure it. Name a metabolite the model ' ...
             'exchanges, or add the reaction to your own model and name ' ...
             'it.'], given(i).id);
    end
    if ~isscalar(exchanges)
        error('minimalReactome:ambiguousProduct', ...
            ['''%s'' has %d boundary reactions, so which one its production ' ...
             'should be measured through is ambiguous. Name the reaction ' ...
             'rather than the metabolite.'], given(i).id, numel(exchanges));
    end
    coefficient = model.S(at, exchanges);
    if coefficient > 0
        error('minimalReactome:productIsSource', ...
            ['%s is written as a source of ''%s'' (coefficient %+g), so ' ...
             'flux through it is uptake and requiring production would mean ' ...
             'a negative bound rather than a positive one. Reverse the ' ...
             'reaction, or name a reaction whose forward direction is the ' ...
             'production meant.'], model.rxns{exchanges}, given(i).id, ...
            coefficient);
    end
    products(end+1) = struct('rxn', model.rxns{exchanges}, 'rate', rate, ...
        'asked', given(i).id);                                    %#ok<AGROW>
    fprintf('minimalReactome: the production of %s is required through %s.\n', ...
        given(i).id, model.rxns{exchanges});
end
end


function required = requiredGrowthIn(model, medium, biomass, minGrowth, fraction)
% Returns the growth rate to impose as a lower bound in one medium.
if isempty(biomass)
    required = 0;
    return
end
if ~isempty(minGrowth)
    required = minGrowth;
    return
end
if isempty(fraction)
    fraction = 1;
end
if fraction <= 0
    required = 0;
    return
end
constrained = applyMedium(model, medium);
constrained = changeObjective(constrained, biomass);
solution = optimizeCbModel(constrained);
if solution.stat ~= 1
    error('minimalReactome:mediumInfeasible', ...
        ['The model cannot be solved at all in medium ''%s'', so there is ' ...
         'no growth rate to take a fraction of.'], medium.label);
end
required = fraction * solution.f;
end


function checkTasksArePossible(model, tasks)
% Stops a task the model cannot perform from being required of a subset.
for i = 1:numel(tasks)
    if isfield(tasks(i), 'shouldFail') && tasks(i).shouldFail
        continue
    end
    outcome = checkTask(model, tasks(i));
    if ~outcome.feasible
        error('minimalReactome:taskImpossible', ...
            ['The model cannot perform task ''%s'' as it stands, so no ' ...
             'subset of it can. Gap-fill towards the task first, with ' ...
             'gapfillForTasks.'], tasks(i).id);
    end
end
end


function checkGrowthIsPossible(model, medium, biomass, required, products)
% Stops an unreachable requirement from being blamed on the minimisation.
if required <= 0 && isempty(products)
    return
end
constrained = growthCondition(model, medium, biomass, required, products);
if ~isFeasible(constrained)
    error('minimalReactome:growthImpossible', ...
        ['The model itself cannot meet the requirement in medium ''%s'', ' ...
         'so no subset of it can. Lower growthFraction or minGrowth, check ' ...
         'the medium, or ask for less production.'], medium.label);
end
end


function constrained = growthCondition(model, medium, biomass, required, products)
% Applies the medium, the growth floor and the product floors.
constrained = applyMedium(model, medium);
if ~isempty(biomass) && required > 0
    at = findRxnIDs(constrained, biomass);
    if at > 0
        constrained.lb(at) = required;
        constrained.ub(at) = max(constrained.ub(at), required);
    end
end
for i = 1:numel(products)
    at = findRxnIDs(constrained, products(i).rxn);
    if at > 0
        constrained.lb(at) = products(i).rate;
        constrained.ub(at) = max(constrained.ub(at), products(i).rate);
    end
end
end


function [conditions, keptButBlocked, forceable] = buildConditions(base, media, ...
    biomass, requiredGrowth, products, tasks, keepRxns, tol)
% Returns the conditions to impose, each as a label and a constrained model.
%
% Growth, production and the medium are one condition, solved as one flux
% vector. A task is not: it closes the model's own exchanges and opens its
% own, so each task is a different environment, and a network that performs
% two tasks does so with two different flux distributions.
conditions = cell(0, 1);
keptButBlocked = struct('label', {}, 'rxns', {});
forceable = cell(0, 1);

for k = 1:numel(media)
    % A growth condition with no floor and nothing to produce is satisfied
    % by the all-zero flux vector, so it constrains nothing and only costs a
    % copy of the network.
    if requiredGrowth(k) <= 0 && isempty(products)
        if ~isempty(media(k).rxns)
            fprintf(['minimalReactome: medium ''%s'' requires nothing of the ' ...
                     'network - no growth rate and no product - so it ' ...
                     'constrains nothing and is dropped. Each task defines ' ...
                     'its own medium.\n'], media(k).label);
        end
        continue
    end
    constrained = growthCondition(base, media(k), biomass, requiredGrowth(k), ...
        products);
    conditions{end+1, 1} = struct('label', media(k).label, ...
        'model', constrained, 'isGrowth', true, 'at', k);         %#ok<AGROW>

    % A reaction the caller insists on keeping is forced to carry flux only
    % where it can: in a condition that blocks it, nothing can make it work
    % and forcing it would make the whole problem infeasible.
    [usable, stuck] = splitForceable(constrained, keepRxns, tol);
    forceable{end+1, 1} = usable;                                 %#ok<AGROW>
    if ~isempty(stuck)
        keptButBlocked(end+1) = struct('label', media(k).label, ...
            'rxns', {stuck});                                     %#ok<AGROW>
        fprintf(['minimalReactome: %d of the reactions in keepRxns cannot ' ...
                 'carry %g of flux in ''%s''. They are kept, but nothing ' ...
                 'forces them to do anything there.\n'], ...
            numel(stuck), tol, media(k).label);
    end
end

for i = 1:numel(tasks)
    if isfield(tasks(i), 'shouldFail') && tasks(i).shouldFail
        fprintf(['minimalReactome: task ''%s'' is a negative control and ' ...
                 'cannot be required of a minimal reactome, so it is ' ...
                 'skipped.\n'], tasks(i).id);
        continue
    end
    if ~taskIsForced(tasks(i))
        fprintf(['minimalReactome: task ''%s'' forces no flux, so the ' ...
                 'all-zero solution performs it and requiring it ' ...
                 'constrains nothing.\n'], tasks(i).id);
    end
    constrained = taskConstraints(base, tasks(i));
    conditions{end+1, 1} = struct('label', tasks(i).id, ...
        'model', constrained, 'isGrowth', false, 'at', i);         %#ok<AGROW>
    forceable{end+1, 1} = cell(0, 1);                              %#ok<AGROW>
end

if isempty(conditions)
    error('minimalReactome:noConditions', ...
        ['No condition survived, so nothing is required of the result. ' ...
         'Give a growth requirement, a product or a task that forces flux.']);
end
labels = cellfun(@(c) c.label, conditions, 'UniformOutput', false);
if numel(unique(labels)) ~= numel(labels)
    error('minimalReactome:duplicateCondition', ...
        ['Two conditions share a label, which would make them one copy of ' ...
         'the network rather than two. Media labels and task identifiers ' ...
         'must all be distinct.']);
end
end


function [usable, stuck] = splitForceable(constrained, keepRxns, tol)
% Splits keepRxns into those that can carry tol here and those that cannot.
usable = cell(0, 1);
stuck = cell(0, 1);
for i = 1:numel(keepRxns)
    if canCarry(constrained, keepRxns{i}, tol)
        usable{end+1, 1} = keepRxns{i};                           %#ok<AGROW>
    else
        stuck{end+1, 1} = keepRxns{i};                            %#ok<AGROW>
    end
end
end


function able = canCarry(model, rxnID, tol)
% Returns whether a reaction can reach tol in either direction here.
at = findRxnIDs(model, rxnID);
if at == 0
    able = false;
    return
end
able = false;
for sense = {'max', 'min'}
    probe = changeObjective(model, rxnID);
    probe.osenseStr = sense{1};
    solution = optimizeCbModel(probe, sense{1});
    if solution.stat == 1 && abs(solution.f) >= tol
        able = true;
        return
    end
end
end


function [blocked, essential] = preprocessConditions(conditions, originalRxns, ...
    keepRxns, tol, seed)
% Decides before the MILP what the MILP would have had to decide.
%
% Both passes are exact, so the optimum is unchanged - only the number of
% binaries the solver has to branch on. Blocked is an intersection over
% conditions and essential a union, and the asymmetry is the point: a
% reaction useless under the growth medium may be the only route a task has,
% while a reaction some task cannot do without has to be in the answer
% whatever the growth condition thinks of it.
% Any nonzero flux counts, which is what fastcc does; see below.
detectionCutoff = getCobraSolverParams('LP', 'feasTol');

blocked = [];
for k = 1:numel(conditions)
    % Count a reaction as able to carry flux if it can carry *any*, not if
    % it can reach tol. The difference decides whether a reaction carrying
    % a load-bearing trace is deleted here, and deleting one breaks the
    % chain it belonged to - the same trap inclusionCutoff exists for.
    % Erring towards keeping costs one binary; erring the other way costs
    % correctness. spectraCC reports indices into the model it was given,
    % which holds the task's own scaffolding as well as the network.
    consistent = spectraCC(conditions{k}.model, tol, 'stoichiometry', seed, ...
        detectionCutoff);
    here = setdiff(originalRxns, conditions{k}.model.rxns(consistent));
    fprintf('minimalReactome: %s: %d of %d reactions cannot carry flux.\n', ...
        conditions{k}.label, numel(here), numel(originalRxns));
    if isempty(blocked)
        blocked = here;
    else
        blocked = intersect(blocked, here);
    end
end
blocked = setdiff(blocked, keepRxns);

candidates = setdiff(originalRxns, union(blocked, keepRxns));
essential = {};
for k = 1:numel(conditions)
    remaining = setdiff(candidates, essential);
    if isempty(remaining)
        break
    end
    found = essentialIn(conditions{k}.model, remaining, tol);
    fprintf('minimalReactome: %s: %d reactions are essential to it.\n', ...
        conditions{k}.label, numel(found));
    essential = union(essential, found);
end

fprintf(['minimalReactome: preprocessing found %d blocked everywhere ' ...
         '(deleted) and %d essential somewhere (no binary), leaving %d ' ...
         'binaries to decide.\n'], numel(blocked), numel(essential), ...
    numel(originalRxns) - numel(blocked) - numel(essential) - numel(keepRxns));
blocked = blocked(:);
essential = essential(:);
end


function essential = essentialIn(constrained, candidates, tol)
% Returns the reactions this condition becomes infeasible without.
%
% The test is feasibility, not growth: with the requirements already
% expressed as bounds, "the condition is infeasible without this reaction"
% is "the required rate cannot be reached without it". Only a reaction
% carrying flux in some feasible solution can be essential, so the
% candidates are first shrunk by minimising total flux: anything a
% parsimonious solution leaves at zero has an alternative and is dropped,
% many at a time, for one LP. Only the survivors are knocked out one by one.
essential = {};
if ~isFeasible(constrained)
    return
end

pool = candidates(ismember(candidates, constrained.rxns));
at = findRxnIDs(constrained, pool);

% A reaction whose own bounds exclude zero cannot be switched off at all, so
% it is in every solution by construction. Knocking one out would not even
% be the right question: the biomass reaction carries the growth requirement
% in its lower bound, so deleting it deletes the requirement rather than
% violating it, and the condition comes back feasible.
forced = constrained.lb(at) > 0 | constrained.ub(at) < 0;
essential = pool(forced);
pool = pool(~forced);

previous = -1;
while ~isempty(pool) && previous ~= numel(pool)
    previous = numel(pool);
    n = numel(constrained.rxns);
    [~, x] = minNet(constrained, zeros(n, 1), ones(n, 1), tol, 1, 'minNetLP');
    if isempty(x)
        break
    end
    flux = x(1:n);
    pool = pool(abs(flux(findRxnIDs(constrained, pool))) > 1e-9);
end

for i = 1:numel(pool)
    probe = changeRxnBounds(constrained, pool{i}, 0, 'b');
    if ~isFeasible(probe)
        essential{end+1, 1} = pool{i};                            %#ok<AGROW>
    end
end
essential = essential(:);
end


function feasible = isFeasible(model)
% Returns whether the constraints can be satisfied at all.
probe = model;
probe.c(:) = 0;
solution = optimizeCbModel(probe);
feasible = solution.stat == 1;
end


function [joint, groups, copyOf] = replicateConditions(conditions, originalRxns, ...
    blocked, counted)
% Builds one copy of the network per condition, side by side in one model.
%
% No metabolite is shared between the copies, so each condition is satisfied
% by its own flux vector; the binaries are what ties them together. With a
% single condition no suffix is added and the copy is the model itself, so
% the common case pays nothing for the generality.
nConditions = numel(conditions);
single = nConditions == 1;
keepMask = ~ismember(originalRxns, blocked);

joint = struct();
joint.rxns = cell(0, 1);
joint.mets = cell(0, 1);
joint.S = sparse(0, 0);
joint.lb = [];
joint.ub = [];
joint.c = [];
joint.b = [];
joint.csense = char([]);
joint.osenseStr = 'max';
joint.modelID = 'minimalReactome';
joint.description = 'minimalReactome';
joint.rev = [];

copyOf = cell(nConditions, 1);
countedAt = cell(nConditions, 1);
for k = 1:nConditions
    constrained = conditions{k}.model;
    suffix = '';
    if ~single
        suffix = ['__' conditions{k}.label];
    end

    % A task adds its own exchange reactions to the model. They are
    % scaffolding rather than network, so they are carried into the joint
    % problem without a binary and never reported.
    [isOriginal, at] = ismember(constrained.rxns, originalRxns);
    keep = true(numel(constrained.rxns), 1);
    keep(isOriginal) = keepMask(at(isOriginal));

    rows = numel(joint.mets);
    cols = numel(joint.rxns);
    joint.mets = [joint.mets; strcat(constrained.mets(:), suffix)];
    joint.rxns = [joint.rxns; strcat(constrained.rxns(keep), suffix)];
    joint.S = blkdiag(joint.S, constrained.S(:, keep));
    joint.lb = [joint.lb; constrained.lb(keep)];
    joint.ub = [joint.ub; constrained.ub(keep)];
    joint.c = [joint.c; zeros(sum(keep), 1)];
    joint.rev = [joint.rev; double(constrained.lb(keep) < 0)];
    joint.b = [joint.b; zeros(numel(constrained.mets), 1)];
    joint.csense = [joint.csense; repmat('E', numel(constrained.mets), 1)];
    assert(rows + numel(constrained.mets) == numel(joint.mets));

    % Where every reaction of the model ended up in this condition's copy.
    place = zeros(numel(originalRxns), 1);
    kept = constrained.rxns(keep);
    [found, at] = ismember(originalRxns, kept);
    place(found) = at(found) + cols;
    copyOf{k} = place;
    countedAt{k} = place(ismember(originalRxns, counted));
end

% One binary per counted reaction, governing its copy in every condition.
groups = cell(numel(counted), 1);
countedRows = find(ismember(originalRxns, counted));
for g = 1:numel(counted)
    members = zeros(nConditions, 1);
    for k = 1:nConditions
        members(k) = copyOf{k}(countedRows(g));
    end
    groups{g} = members(members > 0);
end

% Every condition is built from the same model and none of them removes a
% reaction, so today no reaction ends up without a copy. If one ever did, a
% binary governing nothing would be free to sit at one and inflate the
% count, so it is dropped instead.
empty = cellfun(@isempty, groups);
if any(empty)
    fprintf(['minimalReactome: %d reaction(s) appear in no condition, so ' ...
             'nothing requires them and they are dropped.\n'], sum(empty));
    groups = groups(~empty);
end
end


function kept = keptReactions(reacInd, joint, originalRxns, copyOf, protectedRxns)
% Maps a solution over the joint model back to the model's own reactions.
%
% A reaction is kept if any of its copies is, which is what the shared binary
% already enforces; reading it back this way also catches a copy kept for
% carrying a trace of flux with its binary off.
included = false(numel(originalRxns), 1);
for k = 1:numel(copyOf)
    at = copyOf{k};
    has = at > 0;
    included(has) = included(has) | reacInd(at(has));
end
kept = originalRxns(included);
kept = union(kept, protectedRxns);
kept = originalRxns(ismember(originalRxns, kept));
assert(numel(joint.rxns) >= numel(kept));
end


function blank = emptyResult()
% Returns one result entry with every field present and empty.
blank = struct('model', [], 'kept', {cell(0,1)}, 'removed', {cell(0,1)}, ...
    'growth', [], 'requiredGrowth', [], 'production', [], ...
    'requiredProduction', [], 'conditions', {cell(0,1)}, ...
    'essential', {cell(0,1)}, 'keptButBlocked', [], ...
    'unsatisfied', {cell(0,1)}, 'summary', '');
end


function out = firstLine(text)
% Returns the first line of a multi-line string.
lines = strsplit(text, newline);
out = lines{1};
end


function direction = settleDirections(model, core, tol)
% Settles which way each forced reaction should carry its flux.
%
% A reaction that must carry at least tol has to be told which sign that
% flux takes, and for a reversible one either would do. The direction is
% read off a flux distribution that gets every forced reaction working at
% once, built by blending a forward-biased solution with a reverse-biased
% one until nothing is left idle - the same phase spectraME runs, and the
% reason a core set is settled rather than guessed one reaction at a time.
n = numel(model.rxns);
direction = zeros(n, 1);
if isempty(core)
    return
end

isCore = false(n, 1);
isCore(core) = true;
remaining = isCore;
flux = zeros(n, 1);

while any(remaining)
    flux1 = forward(model, remaining, tol, 1);
    assertForcedFeasible(flux1, sum(isCore), tol);
    if sum(abs(flux)) == 0
        flux = flux1;
    else
        c1 = unifrnd(0.45, 0.55, 1);
        flux = (c1 * flux) + ((1 - c1) * flux1);
    end
    remaining(isCore & abs(flux) >= tol * 1e-1) = 0;
    if ~any(remaining)
        break
    end
    flux2 = reverse(model, remaining, tol, 1);
    assertForcedFeasible(flux2, sum(isCore), tol);
    c1 = unifrnd(0.45, 0.55, 1);
    flux = (c1 * flux) + ((1 - c1) * flux2);
    remaining(isCore & abs(flux) >= tol * 1e-1) = 0;
end

direction(isCore & flux > 0) = 1;
direction(isCore & flux < 0) = -1;
end


function assertForcedFeasible(flux, nCore, tol)
% Stops a direction LP that came back infeasible from becoming a hang.
if any(isnan(flux))
    error('minimalReactome:infeasibleKeeps', ...
        ['No flux distribution lets all %d forced reactions carry at least ' ...
         'tol=%g at once, although each can on its own. The reactions in ' ...
         'keepRxns are mutually inconsistent in some condition.'], nCore, tol);
end
end


function report = reportReactome(base, originalRxns, kept, media, biomass, ...
    requiredGrowth, products, tasks, conditions, essential, keptButBlocked, ...
    tol, check)
% Builds the minimal reactome as a model and measures what it still does.
removed = setdiff(originalRxns, kept);
minimal = removeRxns(base, removed, 'metFlag', false);

report = emptyResult();
report.model = minimal;
report.kept = kept(:);
report.removed = originalRxns(ismember(originalRxns, removed));
report.conditions = cellfun(@(c) c.label, conditions, 'UniformOutput', false);
report.essential = essential(:);
report.keptButBlocked = keptButBlocked;
report.requiredProduction = products;

growth = struct('label', {}, 'reached', {}, 'required', {});
production = struct('label', {}, 'rxn', {}, 'carried', {});
unsatisfied = cell(0, 1);

if check
    for k = 1:numel(media)
        if requiredGrowth(k) <= 0 && isempty(products)
            continue
        end
        constrained = growthCondition(minimal, media(k), biomass, 0, ...
            struct('rxn', {}, 'rate', {}, 'asked', {}));
        reached = 0;
        if ~isempty(biomass) && ismember(biomass, minimal.rxns)
            constrained = changeObjective(constrained, biomass);
            solution = optimizeCbModel(constrained);
            if solution.stat == 1
                reached = solution.f;
            end
        end
        growth(end+1) = struct('label', media(k).label, ...
            'reached', reached, 'required', requiredGrowth(k));   %#ok<AGROW>

        % Read the products off the same solve the requirement is met in,
        % rather than maximising each on its own: they have to hold in one
        % flux distribution, which is how they were required.
        met = true;
        if requiredGrowth(k) > 0 && reached < requiredGrowth(k) * (1 - 1e-6)
            met = false;
        end
        if ~isempty(products)
            together = growthCondition(minimal, media(k), biomass, ...
                requiredGrowth(k), products);
            if isFeasible(together)
                probe = together;
                probe.c(:) = 0;
                solution = optimizeCbModel(probe);
                for i = 1:numel(products)
                    at = findRxnIDs(probe, products(i).rxn);
                    carried = 0;
                    if at > 0 && solution.stat == 1
                        flux = solution.v;
                        carried = flux(at);
                    end
                    production(end+1) = struct('label', media(k).label, ...
                        'rxn', products(i).rxn, 'carried', carried); %#ok<AGROW>
                end
            else
                met = false;
            end
        end
        if ~met
            unsatisfied{end+1, 1} = media(k).label;              %#ok<AGROW>
        end
    end

    for i = 1:numel(tasks)
        if isfield(tasks(i), 'shouldFail') && tasks(i).shouldFail
            continue
        end
        outcome = checkTask(minimal, tasks(i));
        if ~outcome.ok
            unsatisfied{end+1, 1} = tasks(i).id;                  %#ok<AGROW>
        end
    end
end

report.growth = growth;
report.production = production;
report.requiredGrowth = requiredGrowth(:);
report.unsatisfied = unsatisfied;
report.summary = summarise(report, numel(originalRxns), tol);
end


function text = summarise(report, nTotal, tol)
% Returns a printable account of what the minimal reactome is and does.
lines = cell(0, 1);
lines{end+1} = sprintf('minimal reactome: %d of %d reactions (tol %g)', ...
    numel(report.kept), nTotal, tol);
for i = 1:numel(report.growth)
    lines{end+1} = sprintf('  growth on %s: %.6g (needed %.6g)', ...
        report.growth(i).label, report.growth(i).reached, ...
        report.growth(i).required);                               %#ok<AGROW>
end
for i = 1:numel(report.production)
    lines{end+1} = sprintf('  %s carries %.6g in %s', ...
        report.production(i).rxn, report.production(i).carried, ...
        report.production(i).label);                              %#ok<AGROW>
end
if numel(report.conditions) > 1
    lines{end+1} = sprintf('  %d conditions: %s', numel(report.conditions), ...
        strjoin(report.conditions', ', '));
end
if ~isempty(report.essential)
    lines{end+1} = sprintf('  %d reactions proved essential before the solve', ...
        numel(report.essential));
end
if ~isempty(report.keptButBlocked)
    stuck = unique([report.keptButBlocked.rxns]);
    lines{end+1} = sprintf('  %d kept but unable to carry flux', numel(stuck));
end
if isempty(report.unsatisfied)
    lines{end+1} = '  every condition met';
else
    lines{end+1} = sprintf('  UNMET: %s', strjoin(report.unsatisfied', ', '));
end
text = strjoin(lines, newline);
end
