function [Model, LPS, selection] = spectraME(model, coreRxns, tol, consType, weights, nSol, altSolMethod, probType, solveTime, remGene, prevSols, seed, indicatorGroups, inclusionCutoff, solverTol, useCoupling)
% USAGE:
%   [Model, LPS, selection] = spectraME(model, coreRxns, tol, consType, weights, nSol, altSolMethod, probType, solveTime, remGene, prevSols, seed, indicatorGroups, inclusionCutoff, solverTol, useCoupling)
%
% INPUTS:
%   model:           Consistent COBRA model structure (run spectraCC first, or
%                    use spectraCCME on inconsistent models).
%   coreRxns:        Core reactions required in the extracted model (numeric
%                    indices, logical mask, or reaction ID cell/string array).
%
% OPTIONAL INPUTS:
%   tol:             Minimum absolute flux required for core reactions (Default: 1e-4).
%   consType:        'stoichiometry' (S*v = 0, default) or 'topology' (S*v >= 0).
%   weights:         Reaction weights (Default: ones(n, 1)). Non-negative costs
%                    for 'minNetLP', 'minNetMILP', and 'optimBiomass'; any
%                    signed real numbers (positive rewards, negative penalties)
%                    for 'tradeOff'.
%   nSol:            Number of alternative solutions (Default: 1).
%   altSolMethod:    Method to generate alternative solutions:
%                    'coreDirection' (default for LP formulations 'minNetLP'
%                    and 'optimBiomass') or 'pathwayExclusion' (default for
%                    MILP formulations 'minNetMILP' and 'tradeOff'; MILP only).
%   probType:        Network inference formulation:
%                    'minNetLP' (default), 'minNetMILP', 'tradeOff', or
%                    'optimBiomass'.
%   solveTime:       Time limit in seconds for MILP solves (Default: 7200).
%   remGene:         Remove unused genes from extracted model (Default: 0).
%   prevSols:        Previously excluded binary solution sets for MILP (Default: {}).
%   seed:            Optional RNG seed for reproducibility (Default: []).
%   indicatorGroups: Optional MILP binary indicator groups (Default: {}).
%   inclusionCutoff: Minimum absolute flux to retain a non-core reaction (Default: derived).
%   solverTol:       LP/MILP feasibility & optimality (and integrality) tolerance,
%                    or 'default' (Default: 1e-9).
%   useCoupling:     Honor coupling constraints model.C, d, dsense (Default: false).
%
% OUTPUTS:
%   Model:           Extracted consistent model (or cell array when nSol > 1).
%   LPS:             Number of LPs solved (or cell array when nSol > 1).
%   selection:       Group-selection struct when MILP binaries are used.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end
if ~isnumeric(tol) || ~isscalar(tol) || ~isfinite(tol) || tol <= 0
    error('spectraME:badTol', 'tol must be a positive finite scalar.');
end
tol = double(tol);

[~, n] = size(model.S);
if ~exist('consType', 'var') || isempty(consType)
    consType = 'stoichiometry';
end
if strcmp(consType, 'stoichiometry')
    steadystate = 1;
elseif strcmp(consType, 'topology')
    steadystate = 0;
else
    error('spectraME:badConsType', ...
        'consType must be ''stoichiometry'' or ''topology'', not ''%s''.', consType);
end

if ~exist('probType', 'var') || isempty(probType)
    probType = 'minNetLP';
end
if ~ismember(probType, {'minNetLP', 'minNetMILP', 'tradeOff', 'optimBiomass'})
    error('spectraME:badProbType', ...
        'probType must be ''minNetLP'', ''minNetMILP'', ''tradeOff'', or ''optimBiomass'', not ''%s''.', probType);
end

isMilpProb = ismember(probType, {'minNetMILP', 'tradeOff'});
if ~exist('nSol', 'var') || isempty(nSol)
    nSol = 1;
end
if ~exist('altSolMethod', 'var') || isempty(altSolMethod)
    if isMilpProb
        altSolMethod = 'pathwayExclusion';
    else
        altSolMethod = 'coreDirection';
    end
end
if nSol > 1
    if ~ismember(altSolMethod, {'coreDirection', 'pathwayExclusion'})
        error('spectraME:badAltSolMethod', ...
            'altSolMethod must be ''coreDirection'' or ''pathwayExclusion'', not ''%s''.', altSolMethod);
    end
    if strcmp(altSolMethod, 'pathwayExclusion') && ~isMilpProb
        error('spectraME:pathwayExclusionNeedsMilp', ...
            'altSolMethod ''pathwayExclusion'' requires ''minNetMILP'' or ''tradeOff'', not ''%s''.', probType);
    end
end

if ~exist('weights', 'var') || isempty(weights)
    weights       = ones(n, 1);
    userGaveWeights = false;
else
    if ~isnumeric(weights) || any(~isfinite(weights(:)))
        error('spectraME:badWeights', 'weights must be finite numeric values.');
    end
    weights       = double(weights(:));
    userGaveWeights = true;
end
if ~exist('solveTime', 'var') || isempty(solveTime)
    solveTime = 7200;
end
if ~exist('remGene', 'var') || isempty(remGene)
    remGene = 0;
end
if ~exist('prevSols', 'var') || isempty(prevSols)
    prevSols = {};
end
if exist('seed', 'var') && ~isempty(seed)
    rng(seed);
end
if ~exist('indicatorGroups', 'var')
    indicatorGroups = {};
end
if ~exist('inclusionCutoff', 'var')
    inclusionCutoff = [];
end
if ~exist('solverTol', 'var')
    solverTol = [];
end
if ~exist('useCoupling', 'var') || isempty(useCoupling)
    useCoupling = false;
else
    useCoupling = logical(useCoupling);
end

% Validate indicatorGroups and weights
[indicatorGroups, milpWeights] = localResolveGroupArgument( ...
    indicatorGroups, weights, userGaveWeights, model.rxns, probType);
if isempty(indicatorGroups) && numel(weights) ~= n
    error('spectraME:badWeights', ...
        'weights must have one entry per reaction (%d), got %d.', n, numel(weights));
elseif ~isempty(indicatorGroups) && numel(milpWeights) ~= numel(indicatorGroups)
    error('spectraME:groupWeightCount', ...
        ['With indicatorGroups, weights must hold one entry per reaction (%d) ' ...
         'or one entry per group (%d), got %d.'], n, numel(indicatorGroups), numel(weights));
end
if ~strcmp(probType, 'tradeOff') && (any(weights < 0) || any(milpWeights < 0))
    error('spectraME:negativeWeights', ...
        'weights must be non-negative for probType ''%s''; use ''tradeOff'' for signed weights.', probType);
end

% Validate solverTol and inclusionCutoff up front
lockTol = localResolveLockTol(tol, solverTol, inclusionCutoff);

% Canonicalize reverse-only reactions (ub <= 0) on internal working copy
origModel = model;
coreMask  = localParseCoreMask(coreRxns, model.rxns, n);

IrR = model.ub <= 0;
tempUb = model.ub(IrR);
model.S(:, IrR) = -model.S(:, IrR);
if isfield(model, 'C') && ~isempty(model.C)
    model.C(:, IrR) = -model.C(:, IrR);
end
model.ub(IrR) = -model.lb(IrR);
model.lb(IrR) = -tempUb;
model.lb  = full(double(model.lb(:)));
model.ub  = full(double(model.ub(:)));
model.rev = double(model.lb < 0);

% =========================================================================
% Stage 1: Core Direction Discovery (Skipped When Core Is Empty or All-Irrev)
% =========================================================================
[direction, witnessFlux, pass0FluxReady, dirLPs] = localDiscoverCoreDirections( ...
    model, coreMask, weights, probType, tol, steadystate, useCoupling, ...
    inclusionCutoff, solverTol, lockTol, false);

LPS = dirLPs;

% =========================================================================
% Stage 2: Primary Network Extraction (j = 1)
% =========================================================================
[reacInd, extLPs, selected, stat] = localRunSingleExtraction( ...
    model, direction, coreMask, weights, milpWeights, probType, tol, steadystate, ...
    solveTime, witnessFlux, pass0FluxReady, prevSols, indicatorGroups, ...
    inclusionCutoff, solverTol, useCoupling);
LPS = LPS + extLPs;

if isempty(reacInd) && any(coreMask) && isempty(prevSols)
    localRaiseInfeasibleCore(sum(coreMask), tol);
end

Model = localSubsetModel(origModel, reacInd, remGene);
selection = localReportSelection(indicatorGroups, selected, model.rxns, direction);

% =========================================================================
% Alternative Solutions (nSol > 1)
% =========================================================================
if nSol > 1
    Model = {Model};
    if strcmp(altSolMethod, 'coreDirection')
        % Skip duplicate solves when no core reaction can carry >= tol in both directions
        canBoth = coreMask & (model.lb <= -lockTol) & (model.ub >= lockTol);
        if ~any(canBoth)
            Model = repmat(Model, nSol, 1);
            LPS   = num2cell([LPS; zeros(nSol - 1, 1)]);
            return;
        end
        LPS = {LPS};
        for j = 2:nSol
            [directionJ, witnessJ, ~, lpsDirJ] = localDiscoverCoreDirections( ...
                model, coreMask, weights, probType, tol, steadystate, useCoupling, ...
                inclusionCutoff, solverTol, lockTol, true);
            [reacIndJ, lpsExtJ, selected] = localRunSingleExtraction( ...
                model, directionJ, coreMask, weights, milpWeights, probType, tol, ...
                steadystate, solveTime, witnessJ, [], {}, indicatorGroups, ...
                inclusionCutoff, solverTol, useCoupling);
            if isempty(reacIndJ) && any(coreMask)
                localRaiseInfeasibleCore(sum(coreMask), tol);
            end
            ModJ  = localSubsetModel(origModel, reacIndJ, remGene);
            Model = [Model; {ModJ}]; %#ok<AGROW>
            LPS   = [LPS; {lpsDirJ + lpsExtJ}]; %#ok<AGROW>
        end
    elseif strcmp(altSolMethod, 'pathwayExclusion')
        for j = 2:nSol
            if isempty(selected) || ~any(selected) || (stat ~= 1 && stat ~= 3)
                break;
            end
            if isempty(indicatorGroups)
                freeIdx  = find(direction == 0);
                prevSols = [prevSols; {freeIdx(selected)}]; %#ok<AGROW>
            else
                prevSols = [prevSols; {find(selected)}]; %#ok<AGROW>
            end
            [reacIndJ, ~, selected, stat] = localRunSingleExtraction( ...
                model, direction, coreMask, weights, milpWeights, probType, tol, ...
                steadystate, solveTime, [], [], prevSols, indicatorGroups, ...
                inclusionCutoff, solverTol, useCoupling);
            if stat ~= 1
                break;
            end
            ModJ  = localSubsetModel(origModel, reacIndJ, remGene);
            Model = [Model; {ModJ}]; %#ok<AGROW>
        end
    end
    if ~isempty(selected)
        selection = localReportSelection(indicatorGroups, selected, model.rxns, direction);
    end
end
end

% =========================================================================
% Stage 1: Core Direction Discovery (Zero Convex Combination, Smart Skips)
% =========================================================================
function [direction, lastFlux, pass0FluxReady, nLPs] = localDiscoverCoreDirections( ...
    model, coreMask, weights, probType, tol, steadystate, useCoupling, ...
    inclusionCutoff, solverTol, lockTol, isAltSol)

n = numel(model.rxns);
direction      = zeros(n, 1);
lastFlux       = [];
pass0FluxReady = [];
nLPs           = 0;

% Skip 1: Empty core set (even if positive lower bounds or coupling exist)
if ~any(coreMask)
    return;
end

% Classify core reactions by bound reachability at tol (0 LPs)
canFwd = (model.ub >= lockTol);
canRev = (model.lb <= -lockTol);

% Fail immediately in 0 LPs if any core reaction's box bounds forbid +/-tol
if any(coreMask & (~canFwd) & (~canRev))
    localRaiseInfeasibleCore(sum(coreMask), tol);
end

onlyFwdCore = coreMask & canFwd & (~canRev);
onlyRevCore = coreMask & (~canFwd) & canRev;
twoWayCore  = coreMask & canFwd & canRev;

direction(onlyFwdCore) = 1;
direction(onlyRevCore) = -1;

% Skip 2: Every core reaction has a single bound-feasible sign — 0 Stage-1 LPs
if ~any(twoWayCore)
    return;
end

% Build non-negative Stage-1 Anchor penalty vector appropriate for probType:
%   - minNetLP / minNetMILP: non-negative cost weights directly
%   - tradeOff (signed weights): penalize negatively weighted reactions (-w_i > 0)
%     while treating positively weighted reactions (w_i >= 0) as zero-penalty
wAnchor   = [];
useAnchor = false;
hasFixed  = any(onlyFwdCore) || any(onlyRevCore) || ~localZeroFluxIsFeasible(model, useCoupling);
if (~isAltSol) && (numel(weights) == n) && hasFixed
    if ismember(probType, {'minNetLP', 'minNetMILP'})
        wAnchor          = weights;
        wAnchor(coreMask) = 0.0;
        useAnchor        = true;
    elseif strcmp(probType, 'tradeOff') && any(weights < 0)
        wAnchor          = zeros(n, 1);
        negW             = (weights < 0);
        wAnchor(negW)    = -weights(negW);
        wAnchor(coreMask) = 0.0;
        useAnchor        = true;
    end
end

if useAnchor
    [~, vAnchor, statAnc, ~, ancLPs] = minNet(model, direction, wAnchor, tol, ...
        steadystate, 'minNetLP', [], [], {}, {}, inclusionCutoff, solverTol, useCoupling, false);
    nLPs = nLPs + ancLPs;
    if statAnc ~= 1 || isempty(vAnchor)
        localRaiseInfeasibleCore(sum(coreMask), tol);
    end
    lastFlux = vAnchor;

    settledRev = twoWayCore & (abs(vAnchor) >= 0.1 * tol);
    direction(settledRev) = sign(vAnchor(settledRev));
    unsettledMask = twoWayCore & (direction == 0);

    if ~any(unsettledMask)
        if all(abs(vAnchor(coreMask)) >= tol)
            if strcmp(probType, 'minNetLP')
                pass0FluxReady = vAnchor;
            end
            return;
        elseif all(abs(vAnchor(coreMask)) >= lockTol)
            return;
        end
        direction(twoWayCore & (abs(vAnchor) < lockTol)) = 0;
        unsettledMask = twoWayCore & (direction == 0);
    end
else
    vAnchor       = zeros(n, 1);
    unsettledMask = twoWayCore;
end

% Iterative mPinned forward/reverse on unsettled two-way core reactions
stalls   = 0;
roundIdx = 0;
while any(unsettledMask)
    roundIdx = roundIdx + 1;
    before   = sum(unsettledMask);

    if isAltSol
        firstIsFwd = (randi(2) == 1);
    else
        firstIsFwd = true;
    end

    % Reuse iterW across forward and reverse within the round
    iterW = zeros(n, 1);
    iterW(unsettledMask) = 1.0 + 0.1 * rand(before, 1);

    passOrder = [firstIsFwd, ~firstIsFwd];
    for p = 1:2
        if ~any(unsettledMask)
            break;
        end
        doFwd = passOrder(p);

        mPinned = model;
        mPinned.lb(direction == 1)  = max(tol, mPinned.lb(direction == 1));
        mPinned.ub(direction == -1) = min(-tol, mPinned.ub(direction == -1));

        nLPs = nLPs + 1;
        if doFwd
            [f, ~] = forward(mPinned, unsettledMask, tol, steadystate, useCoupling, solverTol, iterW);
        else
            [f, ~] = reverse(mPinned, unsettledMask, tol, steadystate, useCoupling, solverTol, iterW);
        end

        if any(isnan(f)) && useAnchor && roundIdx == 1 && p == 1
            direction(twoWayCore & (abs(vAnchor) < lockTol)) = 0;
            unsettledMask = twoWayCore & (direction == 0);
            mPinned = model;
            mPinned.lb(direction == 1)  = max(tol, mPinned.lb(direction == 1));
            mPinned.ub(direction == -1) = min(-tol, mPinned.ub(direction == -1));
            iterW = zeros(n, 1);
            iterW(unsettledMask) = 1.0 + 0.1 * rand(sum(unsettledMask), 1);
            nLPs = nLPs + 1;
            if doFwd
                [f, ~] = forward(mPinned, unsettledMask, tol, steadystate, useCoupling, solverTol, iterW);
            else
                [f, ~] = reverse(mPinned, unsettledMask, tol, steadystate, useCoupling, solverTol, iterW);
            end
        end
        if any(isnan(f))
            localRaiseInfeasibleCore(sum(coreMask), tol);
        end

        lastFlux = f;
        newly = unsettledMask & (abs(f) >= lockTol);
        direction(newly)     = sign(f(newly));
        unsettledMask(newly) = false;
    end

    stalls = assertCoreProgress(unsettledMask, before, stalls, tol);
end
end

% =========================================================================
% Stage 2: Dispatch to Selected Formulation (minNet / tradeOff / optimBiomass)
% =========================================================================
function [reacInd, extLPs, selected, stat] = localRunSingleExtraction( ...
    model, direction, coreMask, weights, milpWeights, probType, tol, steadystate, ...
    solveTime, witnessFlux, pass0FluxReady, prevSols, indicatorGroups, ...
    inclusionCutoff, solverTol, useCoupling)

n        = numel(model.rxns);
extLPs   = 0;
selected = [];
zeroFeas = ~any(direction ~= 0) && localZeroFluxIsFeasible(model, useCoupling);

if strcmp(probType, 'minNetLP')
    [reacInd, ~, stat, ~, extLPs] = minNet(model, direction, weights, tol, ...
        steadystate, 'minNetLP', [], pass0FluxReady, {}, {}, ...
        inclusionCutoff, solverTol, useCoupling, true);

elseif strcmp(probType, 'minNetMILP')
    if zeroFeas && isempty(prevSols) && isempty(indicatorGroups) && all(milpWeights > 0)
        reacInd  = false(n, 1);
        selected = false(n, 1);
        stat     = 1;
        return;
    end
    if ~any(coreMask) && isempty(witnessFlux) && isempty(prevSols) && ~zeroFeas
        witnessFlux = localOptionalWarmStart(model, useCoupling);
    end
    [reacInd, ~, stat, selected] = minNet(model, direction, milpWeights, tol, ...
        steadystate, 'minNetMILP', solveTime, witnessFlux, prevSols, ...
        indicatorGroups, inclusionCutoff, solverTol, useCoupling);

elseif strcmp(probType, 'optimBiomass')
    if zeroFeas && (~any(model.c ~= 0) || all(weights(model.c ~= 0) == 0))
        reacInd = false(n, 1);
        stat    = 1;
        return;
    end
    [reacInd, ~, stat, extLPs] = optimBiomass(model, direction, weights, tol, ...
        steadystate, useCoupling, solverTol, inclusionCutoff);

elseif strcmp(probType, 'tradeOff')
    if zeroFeas && isempty(prevSols) && isempty(indicatorGroups) && all(milpWeights < 0)
        reacInd  = false(n, 1);
        selected = false(n, 1);
        stat     = 1;
        return;
    end
    if ~any(coreMask) && isempty(witnessFlux) && isempty(prevSols) && ~zeroFeas
        witnessFlux = localOptionalWarmStart(model, useCoupling);
    end
    [reacInd, ~, stat, selected] = tradeOff(model, direction, milpWeights, tol, ...
        steadystate, solveTime, witnessFlux, prevSols, indicatorGroups, ...
        inclusionCutoff, solverTol, useCoupling);
end
end

% =========================================================================
% Validation & Helper Utilities
% =========================================================================
function lockTol = localResolveLockTol(tol, solverTol, inclusionCutoff)
FLOOR = 1e-9;
lpTolerances(solverTol); % Validates solverTol immediately
if isempty(solverTol)
    effSolverTol = FLOOR;
elseif (ischar(solverTol) || (isstring(solverTol) && isscalar(solverTol))) && ...
        strcmpi(solverTol, 'default')
    effSolverTol = double(getCobraSolverParams('LP', 'feasTol'));
else
    effSolverTol = max(FLOOR, double(solverTol));
end
if ~isempty(inclusionCutoff)
    if ~isnumeric(inclusionCutoff) || ~isscalar(inclusionCutoff) || ...
            ~isfinite(inclusionCutoff) || inclusionCutoff <= 0
        error('spectraME:badInclusionCutoff', ...
            'inclusionCutoff must be a positive finite scalar.');
    end
end
lockTol = max(0.99 * tol, tol - 10.0 * effSolverTol);
end

% -------------------------------------------------------------------------
function tf = localZeroFluxIsFeasible(model, useCoupling)
tf = all(model.lb <= 0) && all(model.ub >= 0);
if tf && useCoupling && isfield(model, 'C') && ~isempty(model.C)
    d  = reshape(model.d, [], 1);
    ds = reshape(model.dsense, [], 1);
    tf = ~any(((ds == 'G') & (d > 0)) | ((ds == 'L') & (d < 0)) | ((ds == 'E') & (d ~= 0)));
end
end

% -------------------------------------------------------------------------
function x0 = localOptionalWarmStart(model, useCoupling)
x0 = [];
if ~any(model.c ~= 0)
    return;
end
mWarm = model;
if ~useCoupling && isfield(mWarm, 'C')
    mWarm = rmfield(mWarm, intersect(fieldnames(mWarm), {'C', 'd', 'dsense', 'ctrs'}));
end
sol = optimizeCbModel(mWarm);
if isfield(sol, 'stat') && sol.stat == 1 && isfield(sol, 'x') && ~isempty(sol.x)
    x0 = sol.x;
end
end

% -------------------------------------------------------------------------
function coreMask = localParseCoreMask(coreRxns, rxns, n)
if isempty(coreRxns)
    coreMask = false(n, 1);
elseif iscell(coreRxns) || isstring(coreRxns) || ischar(coreRxns)
    if ischar(coreRxns), coreRxns = {coreRxns}; end
    coreMask = ismember(rxns, cellstr(coreRxns(:)));
elseif islogical(coreRxns) && numel(coreRxns) == n
    coreMask = coreRxns(:);
else
    coreMask = ismember((1:n)', coreRxns(:));
end
end

% -------------------------------------------------------------------------
function localRaiseInfeasibleCore(nCore, tol)
error('spectraME:infeasibleCore', ...
    ['No flux distribution lets all %d core reactions carry at least ' ...
     'tol=%g at once. Either the core set is mutually inconsistent, or ' ...
     'some of its reactions are blocked in this model. spectraCCME ' ...
     'handles the latter by dropping them and returning them as ' ...
     'BlockedCoreRxns.'], nCore, tol);
end

% -------------------------------------------------------------------------
function subModel = localSubsetModel(origModel, reacInd, remGene)
n = numel(origModel.rxns);
keep = false(n, 1);
if islogical(reacInd)
    keep = reacInd(:);
else
    keep(reacInd) = true;
end
if all(keep)
    subModel = origModel;
else
    subModel = removeRxns(origModel, origModel.rxns(~keep));
end
if remGene
    subModel = removeUnusedGenes(subModel);
end
end

% -------------------------------------------------------------------------
function [groups, milpWeights] = localResolveGroupArgument(given, weights, userGaveWeights, rxns, probType)
groups      = {};
milpWeights = weights;
if isempty(given)
    return;
end
if ~ismember(probType, {'minNetMILP', 'tradeOff'})
    error('spectraME:groupsWithoutBinaries', ...
        ['indicatorGroups groups the MILP''s binaries, and ''%s'' has ' ...
         'none. Use ''minNetMILP'' or ''tradeOff''.'], probType);
end
if ~iscell(given)
    error('spectraME:badIndicatorGroups', ...
        'indicatorGroups must be a cell array, of reaction index vectors or of reaction IDs.');
end
named = all(cellfun(@(g) ischar(g) || (isstring(g) && isscalar(g)), given(:)'));
if named
    ids = cellstr(given(:));
    [found, at] = ismember(ids, rxns);
    if ~all(found)
        missing = ids(~found);
        error('spectraME:unknownIndicatorReaction', ...
            'No reaction %s in the model, so it cannot be given a binary.', missing{1});
    end
    groups = num2cell(at(:));
else
    groups = cellfun(@(g) g(:), given(:), 'UniformOutput', false);
end
if ~userGaveWeights
    milpWeights = ones(numel(groups), 1);
elseif numel(weights) == numel(rxns) && numel(groups) ~= numel(weights)
    milpWeights = cellfun(@(g) sum(weights(g)), groups);
    milpWeights = milpWeights(:);
end
end

% -------------------------------------------------------------------------
function selection = localReportSelection(indicatorGroups, selected, rxns, direction)
if isempty(selected)
    selection = [];
    return;
end
[groupOf, nGroups, ungoverned] = resolveIndicatorGroups(indicatorGroups, ...
    direction == 0, numel(rxns));
groups = cell(nGroups, 1);
for g = 1:nGroups
    groups{g} = find(groupOf == g);
end
on = false(numel(rxns), 1);
for g = reshape(find(selected(:)), 1, [])
    on(groups{g}) = true;
end
selection = struct('groups', {groups}, 'selected', {logical(selected(:))}, ...
    'selectedRxns', {rxns(on)}, 'ungoverned', {ungoverned});
end
