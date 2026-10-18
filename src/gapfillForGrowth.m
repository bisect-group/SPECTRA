function result = gapfillForGrowth(draft,universal,media,minGrowth,objectiveRxn,tol,probType,freeRxns,solveTime,seed)
% Adds reactions until a draft grows in every medium of a panel
%
% USAGE:
%
%    result = gapfillForGrowth(draft,universal,media,minGrowth,objectiveRxn,tol,probType,freeRxns,solveTime,seed)
%
% INPUTS:
%    draft:             COBRA model structure that does not yet grow
%                       everywhere.
%    universal:         COBRA model structure to borrow reactions from. The
%                       draft's reactions are expected to be among them;
%                       any that are not cannot be kept and are reported.
%    media:             Struct array, one entry per medium, with fields:
%
%                         * rxns - exchange reaction IDs to open, as a cell
%                           array of char
%                         * lb - lower bound for each, as a numeric vector
%                         * label - name for the medium (optional; defaults
%                           to its position, as 'medium 1')
%
%                       Every exchange not named is closed. A medium is
%                       defined as much by what it withholds as by what it
%                       supplies.
%
% OPTIONAL INPUTS:
%    minGrowth:         Growth rate to demand, as a lower bound on the
%                       objective reaction (Default: 0.1). Pass [] to demand
%                       no particular rate, in which case anything above
%                       1e-6 counts as growth. A plausible rate rather than
%                       a numerical floor: the point is a model that grows,
%                       and one scraping along at 1e-6 does not. Replace it
%                       with a measured rate where you have one.
%    objectiveRxn:      Reaction whose flux counts as growth (Default: the
%                       single reaction in the universal model's objective;
%                       an error if that is ambiguous).
%    tol:               Flux threshold handed to spectraME (Default: 1e-4).
%    probType:          'minNetMILP' (Default) or 'minNetLP'. 'tradeOff' is
%                       refused: it requires every included reaction to
%                       carry at least tol, while a growth solution needs
%                       trace fluxes far below that, so it is infeasible at
%                       every tolerance.
%    freeRxns:          Further reaction IDs that cost nothing to include,
%                       beyond the draft's own and the universal model's
%                       exchanges, which are always free. (Default: none)
%    solveTime:         Maximum runtime per MILP, in seconds (Default: 300).
%    seed:              Seed for spectraME's random coefficients. Leave it
%                       out to draw them unseeded. (Default: none)
%
% OUTPUTS:
%    result:            Structure with fields:
%
%                         * model - the gap-filled model, the draft plus
%                           what was added
%                         * added - reaction IDs taken from the universal
%                           model, sorted
%                         * satisfied - labels of the media the result
%                           grows in
%                         * unsatisfied - labels of the media it still
%                           does not
%                         * failed - struct mapping a medium's label to why
%                           its solve failed, where one did
%                         * missingFromUniversal - draft reactions absent
%                           from the universal model, which could not be
%                           carried over
%                         * summary - a readable multi-line summary, as char
%
% EXAMPLE:
%
%    ::
%
%        M9 = {'EX_pi(e)','EX_h2o(e)','EX_h(e)','EX_nh4(e)','EX_o2(e)','EX_co2(e)'};
%        carbon = {'EX_glc(e)','EX_ac(e)','EX_succ(e)'};
%        media = struct('label',{},'rxns',{},'lb',{});
%        for i = 1:numel(carbon)
%            media(i).label = carbon{i};
%            media(i).rxns  = [M9, carbon(i)];
%            media(i).lb    = [-1000*ones(1,numel(M9)), -10];
%        end
%        result = gapfillForGrowth(draft, universal, media);
%
% NOTE:
%    The growth requirement is a bound on the objective reaction, not a
%    large tol. The two are easy to confuse, because making biomass the
%    core reaction and raising tol to a growth rate does appear to work. It
%    also silently demands that same rate of every other core reaction,
%    which breaks as soon as you supply one. Keeping them apart
%    additionally makes the answer stop moving with tol: on iJO1366 the
%    same 15 reactions come back at every tolerance from 1e-4 to 1e-7.
%
%    Media are worked through in order, and what one medium adds is free
%    for the next, so the panel is filled jointly rather than
%    independently.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('minGrowth', 'var')
    minGrowth = 0.1;
end
if ~exist('objectiveRxn', 'var')
    objectiveRxn = '';
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end
if ~exist('probType', 'var') || isempty(probType)
    probType = 'minNetMILP';
end
if ~exist('freeRxns', 'var') || isempty(freeRxns)
    freeRxns = {};
end
if ~exist('solveTime', 'var') || isempty(solveTime)
    solveTime = 300;
end
if ~exist('seed', 'var')
    seed = [];
end

if ~any(strcmp(probType,{'minNetMILP','minNetLP'}))
    error('gapfillForGrowth:badProbType', ...
        ['probType must be ''minNetMILP'' or ''minNetLP'' to gap-fill ' ...
         'towards growth, not ''%s''. tradeOff requires every included ' ...
         'reaction to carry at least tol, which growth cannot satisfy.'], ...
        probType);
end

media = normaliseMedia(media);
biomass = objectiveReaction(universal, objectiveRxn);

draftIDs = draft.rxns(:);
missingFromUniversal = draftIDs(~ismember(draftIDs, universal.rxns));
if ~isempty(missingFromUniversal)
    warning('gapfillForGrowth:draftNotInUniversal', ...
        ['%d of the draft''s reactions are not in the universal model, so ' ...
         'they cannot be carried into the result. The first is ''%s''.'], ...
        numel(missingFromUniversal), missingFromUniversal{1});
    draftIDs = draftIDs(ismember(draftIDs, universal.rxns));
end

isExc = findExcRxns(universal);
freeIDs = unique([draftIDs; universal.rxns(isExc); freeRxns(:)]);

added = {};
failedLabels = {};
failedWhy = {};

for k = 1:numel(media)
    sofar = rebuildModel(universal, [draftIDs; added]);
    if grows(sofar, media(k), growthFloor(minGrowth))
        continue
    end

    constrained = applyMedium(universal, media(k));
    if ~isempty(minGrowth)
        constrained.lb(findRxnIDs(constrained,biomass)) = minGrowth;
    end

    free = unique([freeIDs; added]);
    weights = ones(numel(constrained.rxns),1);
    weights(ismember(constrained.rxns, free)) = 0;

    try
        filled = spectraME(constrained, findRxnIDs(constrained,biomass), ...
            tol, 'stoichiometry', weights, 1, 'coreDirection', probType, ...
            solveTime, 0, {}, seed);
    catch err
        failedLabels{end+1,1} = media(k).label; %#ok<AGROW>
        failedWhy{end+1,1} = sprintf('%s: %s', err.identifier, err.message); %#ok<AGROW>
        warning('gapfillForGrowth:mediumFailed', ...
            '%s could not be gap-filled (%s)', media(k).label, err.message);
        continue
    end

    newlyAdded = setdiff(filled.rxns, draftIDs);
    added = unique([added; newlyAdded(:)]);
end

model = rebuildModel(universal, [draftIDs; added]);

threshold = growthFloor(minGrowth);
isSatisfied = false(numel(media),1);
for k = 1:numel(media)
    isSatisfied(k) = grows(model, media(k), threshold);
end
labels = {media.label}';

result.model                = model;
result.added                = sort(added);
result.satisfied            = labels(isSatisfied);
result.unsatisfied          = labels(~isSatisfied);
result.missingFromUniversal = missingFromUniversal;
if isempty(failedLabels)
    result.failed = struct();
else
    result.failed = cell2struct(failedWhy, matlab.lang.makeValidName(failedLabels), 1);
end

lines = {sprintf('%d reactions added to a %d-reaction draft', ...
    numel(result.added), numel(draftIDs))};
lines{end+1} = sprintf('  grows in %d of %d media', sum(isSatisfied), numel(media));
if any(~isSatisfied)
    still = labels(~isSatisfied);
    lines{end+1} = sprintf('  still does not grow in: %s', ...
        strjoin(still(1:min(5,end))', ', '));
end
if ~isempty(failedLabels)
    lines{end+1} = sprintf('  %d medium/media could not be solved', numel(failedLabels));
end
result.summary = strjoin(lines, newline);
end


function media = normaliseMedia(media)
% Checks the media struct array and fills in any missing labels
%
% USAGE:
%
%    media = normaliseMedia(media)
%
% INPUTS:
%    media:             Struct array with fields rxns and lb, and
%                       optionally label.
%
% OUTPUTS:
%    media:             The same array, with every entry labelled.

if ~isstruct(media) || ~all(isfield(media,{'rxns','lb'}))
    error('gapfillForGrowth:badMedia', ...
        ['media must be a struct array with fields rxns and lb, and ' ...
         'optionally label.']);
end
for k = 1:numel(media)
    if numel(media(k).rxns) ~= numel(media(k).lb)
        error('gapfillForGrowth:mediumLengthMismatch', ...
            'medium %d names %d reactions but gives %d bounds.', ...
            k, numel(media(k).rxns), numel(media(k).lb));
    end
    if ~isfield(media,'label') || isempty(media(k).label)
        media(k).label = sprintf('medium %d', k);
    end
end
end


function constrained = applyMedium(model,medium)
% Closes every exchange, then opens the ones the medium names
%
% USAGE:
%
%    constrained = applyMedium(model,medium)
%
% INPUTS:
%    model:             COBRA model structure.
%    medium:            One entry of the media struct array.
%
% OUTPUTS:
%    constrained:       The model with that medium applied.

unknown = medium.rxns(~ismember(medium.rxns, model.rxns));
if ~isempty(unknown)
    error('gapfillForGrowth:unknownExchange', ...
        'The medium names %d reaction(s) the model does not have, for example ''%s''.', ...
        numel(unknown), unknown{1});
end
constrained = model;
isExc = findExcRxns(constrained);
constrained.lb(isExc) = 0;
constrained.lb(findRxnIDs(constrained, medium.rxns)) = medium.lb(:);
end


function rxnID = objectiveReaction(model,given)
% Returns the reaction whose flux counts as growth
%
% USAGE:
%
%    rxnID = objectiveReaction(model,given)
%
% INPUTS:
%    model:             COBRA model structure.
%    given:             The caller's choice, or '' to read it off the
%                       model's objective.
%
% OUTPUTS:
%    rxnID:             The reaction identifier.

if ~isempty(given)
    if ~ismember(given, model.rxns)
        error('gapfillForGrowth:noSuchReaction', ...
            'No reaction ''%s'' in the model.', given);
    end
    rxnID = given;
    return
end
inObjective = model.rxns(model.c ~= 0);
if numel(inObjective) ~= 1
    error('gapfillForGrowth:ambiguousObjective', ...
        ['The model has %d reactions in its objective, so which one is ' ...
         'growth is ambiguous. Pass objectiveRxn.'], numel(inObjective));
end
rxnID = inObjective{1};
end


function yes = grows(model,medium,threshold)
% Returns whether the model reaches the growth threshold in a medium
%
% USAGE:
%
%    yes = grows(model,medium,threshold)
%
% INPUTS:
%    model:             COBRA model structure.
%    medium:            One entry of the media struct array.
%    threshold:         The growth rate to reach.
%
% OUTPUTS:
%    yes:               Whether it grows that fast.
%
% NOTE:
%    A medium naming reactions this cut-down model does not have cannot be
%    applied to it, so it certainly does not grow in it.

try
    constrained = applyMedium(model,medium);
catch
    yes = false;
    return
end
solution = optimizeCbModel(constrained);
yes = (solution.stat == 1) && ~isnan(solution.f) && (solution.f >= threshold);
end


function floorValue = growthFloor(minGrowth)
% Returns the growth rate that counts as growing
%
% USAGE:
%
%    floorValue = growthFloor(minGrowth)
%
% INPUTS:
%    minGrowth:         The demanded rate, or [].
%
% OUTPUTS:
%    floorValue:        The rate to test against.

if isempty(minGrowth)
    floorValue = 1e-6;
else
    floorValue = minGrowth;
end
end


function reduced = rebuildModel(universal,keepIDs)
% Returns the universal model cut down to the given reactions
%
% USAGE:
%
%    reduced = rebuildModel(universal,keepIDs)
%
% INPUTS:
%    universal:         COBRA model structure to cut down.
%    keepIDs:           Reaction IDs to retain.
%
% OUTPUTS:
%    reduced:           The reduced model.
%
% NOTE:
%    Orphaned metabolites stay: the draft may still need them, and they
%    cost nothing. COBRA's removeRxns drops them by default, hence metFlag
%    false.

drop = setdiff(universal.rxns, keepIDs);
if isempty(drop)
    reduced = universal;
else
    reduced = removeRxns(universal, drop, 'metFlag', false);
end
end
