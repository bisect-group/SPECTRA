function result = gapfillForTasks(draft, universal, tasks, tol, probType, freeRxns, solveTime, seed, essentialRxns, coreSolve, repair, weights)
% Adds reactions until a draft can perform a list of metabolic tasks
%
% USAGE:
%
%    result = gapfillForTasks(draft, universal, tasks, tol, probType, freeRxns, solveTime, seed, essentialRxns, coreSolve, repair, weights)
%
% INPUTS:
%    draft:             COBRA model structure that cannot yet perform every
%                       task.
%    universal:         COBRA model structure to borrow reactions from.
%    tasks:             Task structures from metabolicTask, as a struct
%                       array.
%
% OPTIONAL INPUTS:
%    tol:               Flux threshold handed to extraction/repair (Default: 1e-4).
%    probType:          Formulation for the core solve: 'minNetMILP' (Default),
%                       'minNetLP', 'tradeOff', or 'optimBiomass'.
%    freeRxns:          Further reaction IDs that cost nothing to include,
%                       beyond the draft's own (Default: {}).
%    solveTime:         Maximum runtime per MILP, in seconds (Default: 300).
%    seed:              Seed for random direction coefficients (Default: []).
%    essentialRxns:     Reaction IDs to use as the core set instead of
%                       computing them via essentialReactionsForTasks (Default: {}).
%    coreSolve:         Whether to run the Step-1 core extraction (Default: true).
%    repair:            Whether to follow up task by task with 'minNetMILP' on
%                       any task still failing after Step 1 (Default: true).
%    weights:           Struct array with fields rxn and weight overriding
%                       reaction weights. For 'minNetMILP', 'minNetLP', and
%                       'optimBiomass', weights must be non-negative (>= 0;
%                       draft/free default to 0 except objective reactions in
%                       'optimBiomass' which default to 1, others default to 1).
%                       For 'tradeOff', any signed real weights are accepted
%                       (draft/free default to 0, candidate additions to -1)
%                       and converted monotonically into non-negative repair
%                       penalties in the Step-2 'minNetMILP' pass.
%
% OUTPUTS:
%    result:            Structure with fields:
%                         * model       - the gap-filled model
%                         * added       - reaction IDs taken from universal, sorted
%                         * satisfied   - IDs of usable tasks the result performs
%                         * unsatisfied - IDs of usable tasks it still does not
%                         * failed      - struct mapping failed solve IDs to reasons
%                         * skipped     - IDs of tasks universal itself cannot perform
%                         * summary     - readable multi-line summary (char)
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end
if ~isnumeric(tol) || ~isscalar(tol) || ~isfinite(tol) || tol <= 0
    error('gapfillForTasks:badTol', 'tol must be a positive finite scalar.');
end
tol = double(tol);

if ~exist('probType', 'var') || isempty(probType)
    probType = 'minNetMILP';
end
if ~ismember(probType, {'minNetLP', 'minNetMILP', 'tradeOff', 'optimBiomass'})
    error('gapfillForTasks:badProbType', ...
        'probType must be ''minNetLP'', ''minNetMILP'', ''tradeOff'', or ''optimBiomass'', not ''%s''.', probType);
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
if ~exist('essentialRxns', 'var')
    essentialRxns = {};
end
if ~exist('coreSolve', 'var') || isempty(coreSolve)
    coreSolve = true;
end
if ~exist('repair', 'var') || isempty(repair)
    repair = true;
end
if ~exist('weights', 'var') || isempty(weights)
    weights = struct('rxn', {}, 'weight', {});
end
if ~isstruct(weights) || (~isempty(weights) && ~all(isfield(weights, {'rxn', 'weight'})))
    error('gapfillForTasks:badWeights', ...
        'weights must be a struct array with fields rxn and weight.');
end

% Validate user weight overrides up front (0 LPs)
localValidateWeightStruct(weights, universal.rxns, probType);

draftIDs = draft.rxns(:);
draftIDs = draftIDs(ismember(draftIDs, universal.rxns));
freeIDs  = unique([draftIDs; freeRxns(:)]);

% Only work on tasks the universal model can actually perform
outcomes = checkTasks(universal, tasks, 0);
isUsable = [outcomes.ok]' & ~[tasks.shouldFail]';
usable   = tasks(isUsable);
skipped  = {tasks(~isUsable & ~[tasks.shouldFail]').id}';
if ~isempty(skipped)
    warning('gapfillForTasks:unusableTasks', ...
        ['%d task(s) the universal model cannot perform are skipped; ' ...
         'there is nothing to add that would help. The first is ''%s''.'], ...
        numel(skipped), skipped{1});
end

failedIDs = {};
failedWhy = {};
added     = {};

% =========================================================================
% Step 1: Core Solve — Extract Around Essential Task Reactions
% =========================================================================
if ~isempty(essentialRxns)
    core = essentialRxns(:);
elseif coreSolve && ~isempty(usable)
    core = essentialReactionsForTasks(universal, usable, {}, 0);
else
    core = {};
end

if ~isempty(core) && coreSolve
    coreWeights = localBuildCoreWeights(universal, freeIDs, weights, probType);
    try
        % Use spectraCCME so blocked reactions or task-scaffold-only core
        % reactions in an unpruned universal model never abort Step 1
        [filled, ~] = spectraCCME(universal, findRxnIDs(universal, core), tol, ...
            'stoichiometry', coreWeights, 1, [], probType, solveTime, 0, seed);
        added = unique([added; setdiff(filled.rxns, draftIDs)]);
    catch err
        failedIDs{end+1, 1} = 'essentialReactions';
        failedWhy{end+1, 1} = sprintf('%s: %s', err.identifier, err.message);
        warning('gapfillForTasks:coreSolveFailed', ...
            'the core solve failed (%s)', err.message);
    end
end

model = rebuildModel(universal, [draftIDs; added]);

% =========================================================================
% Step 2: Task-by-Task Repair Pass via minNetMILP (Non-Negative Penalties)
% =========================================================================
if repair && ~isempty(usable)
    outcomes = checkTasks(model, usable, 0);
    for k = 1:numel(usable)
        if outcomes(k).ok
            continue;
        end
        present = unique([draftIDs; added]);
        try
            [constrained, taskAddedRxns] = taskConstraints(universal, usable(k));
        catch err
            failedIDs{end+1, 1} = usable(k).id; %#ok<AGROW>
            failedWhy{end+1, 1} = sprintf('%s: %s', err.identifier, err.message); %#ok<AGROW>
            continue;
        end

        n = numel(constrained.rxns);
        repairCost = localBuildRepairCosts( ...
            constrained, universal, present, freeRxns, taskAddedRxns, weights, probType);

        [reacInd, ~, stat] = minNet(constrained, zeros(n, 1), repairCost, tol, 1, ...
            'minNetMILP', solveTime);
        if isempty(reacInd) || (~isempty(stat) && stat ~= 1 && stat ~= 3)
            failedIDs{end+1, 1} = usable(k).id; %#ok<AGROW>
            failedWhy{end+1, 1} = sprintf('repair solve returned status %d', stat); %#ok<AGROW>
            continue;
        end

        kept  = constrained.rxns(reacInd);
        kept  = kept(ismember(kept, universal.rxns));
        added = unique([added; setdiff(kept, draftIDs)]);
    end
    model = rebuildModel(universal, [draftIDs; added]);
end

% =========================================================================
% Step 3: Verify & Summarize
% =========================================================================
outcomes  = checkTasks(model, usable, 0);
isOK      = [outcomes.ok]';
usableIDs = {usable.id}';

result.model       = model;
result.added       = sort(added);
result.satisfied   = usableIDs(isOK);
result.unsatisfied = usableIDs(~isOK);
result.skipped     = skipped;
if isempty(failedIDs)
    result.failed = struct();
else
    result.failed = cell2struct(failedWhy, matlab.lang.makeValidName(failedIDs), 1);
end

lines = {sprintf('%d reactions added to a %d-reaction draft', ...
    numel(result.added), numel(draftIDs))};
lines{end+1} = sprintf('  performs %d of %d usable tasks', sum(isOK), numel(usable));
if any(~isOK)
    still = usableIDs(~isOK);
    lines{end+1} = sprintf('  still fails: %s', strjoin(still(1:min(5, end))', ', '));
end
if ~isempty(skipped)
    lines{end+1} = sprintf('  %d task(s) skipped: the universal model cannot do them either', ...
        numel(skipped));
end
result.summary = strjoin(lines, newline);
end

% =========================================================================
% Weight Validation & Formulation-Aware Weight Construction
% =========================================================================
function localValidateWeightStruct(weights, univRxns, probType)
for i = 1:numel(weights)
    if ~ischar(weights(i).rxn) && ~(isstring(weights(i).rxn) && isscalar(weights(i).rxn))
        error('gapfillForTasks:badWeights', ...
            'Each weights entry must specify a reaction ID in field rxn.');
    end
    if ~ismember(char(weights(i).rxn), univRxns)
        error('gapfillForTasks:unknownWeight', ...
            'No reaction ''%s'' in the universal model, so it cannot be weighted.', ...
            char(weights(i).rxn));
    end
    wVal = weights(i).weight;
    if ~isnumeric(wVal) || ~isscalar(wVal) || ~isfinite(wVal)
        error('gapfillForTasks:badWeights', ...
            'Weight for reaction ''%s'' must be a finite numeric scalar.', ...
            char(weights(i).rxn));
    end
    if ~strcmp(probType, 'tradeOff') && (wVal < 0)
        error('gapfillForTasks:negativeWeights', ...
            'Reaction weights must be non-negative for probType ''%s''; use ''tradeOff'' for signed weights.', ...
            probType);
    end
end
end

% -------------------------------------------------------------------------
function coreWeights = localBuildCoreWeights(universal, freeIDs, weights, probType)
nUniv  = numel(universal.rxns);
isFree = ismember(universal.rxns, freeIDs);

if strcmp(probType, 'tradeOff')
    % In tradeOff (max sum w_i z_i), candidate additions default to -1 penalty
    % and draft/free reactions default to 0 (neutral/free).
    coreWeights         = -ones(nUniv, 1);
    coreWeights(isFree) = 0.0;
elseif strcmp(probType, 'optimBiomass')
    % Keep positive default reward (1) on objective reactions even if in draft
    isObj = (universal.c(:) ~= 0);
    if ~any(isObj)
        warning('gapfillForTasks:noObjectiveForOptimBiomass', ...
            'universal.c has no non-zero objective entries; ''optimBiomass'' reduces to L1 flux minimization.');
    end
    coreWeights                  = ones(nUniv, 1);
    coreWeights(isFree & ~isObj) = 0.0;
else
    % minNetMILP / minNetLP: non-negative inclusion costs
    coreWeights         = ones(nUniv, 1);
    coreWeights(isFree) = 0.0;
end

coreWeights = applyWeights(coreWeights, universal.rxns, weights, true);
end

% -------------------------------------------------------------------------
function repairCost = localBuildRepairCosts(constrained, universal, present, freeRxns, taskAddedRxns, weights, probType)
n              = numel(constrained.rxns);
isAlreadyThere = ismember(constrained.rxns, [present(:); taskAddedRxns(:)]);
isFreeNow      = isAlreadyThere | ismember(constrained.rxns, freeRxns(:));

if strcmp(probType, 'tradeOff')
    % Build signed tradeOff evidence vector first, then map monotonically
    % onto non-negative minNetMILP repair penalties:
    %   - Negative weight w < 0  -> penalty -w > 0
    %   - Already-present / task scaffolding with w >= 0 -> 0
    %   - Non-present candidate with w >= 0 -> tiny positive tie-breaker
    %     epsBase / (1 + w) > 0 so higher positive evidence is preferred
    %     without ever switching on idle (zero-flux) binaries.
    signedW            = -ones(n, 1);
    signedW(isFreeNow) = 0.0;
    signedW            = applyWeights(signedW, constrained.rxns, weights, false);

    negMask = (signedW < 0);
    if any(negMask)
        minNeg = min([1.0; -signedW(negMask)]);
    else
        minNeg = 1.0;
    end
    epsBase = 1e-4 * minNeg;

    repairCost          = zeros(n, 1);
    repairCost(negMask) = -signedW(negMask);

    candNonNeg = (~negMask) & (~isAlreadyThere);
    if any(candNonNeg)
        repairCost(candNonNeg) = epsBase ./ (1.0 + signedW(candNonNeg));
    end

elseif strcmp(probType, 'optimBiomass')
    % Do not turn a Step-1 biomass reward weight into a Step-2 repair penalty
    repairCost            = ones(n, 1);
    repairCost(isFreeNow) = 0.0;
    if ~isempty(weights)
        objRxns       = universal.rxns(universal.c(:) ~= 0);
        nonObjWeights = weights(~ismember({weights.rxn}, objRxns));
        repairCost    = applyWeights(repairCost, constrained.rxns, nonObjWeights, false);
    end

else
    % minNetMILP / minNetLP
    repairCost            = ones(n, 1);
    repairCost(isFreeNow) = 0.0;
    repairCost            = applyWeights(repairCost, constrained.rxns, weights, false);
end
end

% -------------------------------------------------------------------------
function cost = applyWeights(cost, rxns, weights, strict)
for i = 1:numel(weights)
    at = find(strcmp(rxns, char(weights(i).rxn)), 1);
    if isempty(at)
        if strict
            error('gapfillForTasks:unknownWeight', ...
                'No reaction ''%s'' in the universal model, so it cannot be weighted.', ...
                char(weights(i).rxn));
        end
        continue;
    end
    cost(at) = double(weights(i).weight);
end
end

% -------------------------------------------------------------------------
function reduced = rebuildModel(universal, keepIDs)
drop = setdiff(universal.rxns, keepIDs);
if isempty(drop)
    reduced = universal;
else
    reduced = removeRxns(universal, drop, 'metFlag', false);
end
end
