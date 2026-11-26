function [Model, BlockedCoreRxns, LPs, selection] = spectraCCME(model, coreRxns, tol, consType, weights, nSol, altSolMethod, probType, solveTime, remGene, seed, useCoupling, inclusionCutoff, solverTol, indicatorGroups)
% USAGE:
%   [Model, BlockedCoreRxns, LPs, selection] = spectraCCME(model, coreRxns, tol, consType, weights, nSol, altSolMethod, probType, solveTime, remGene, seed, useCoupling, inclusionCutoff, solverTol, indicatorGroups)
%
% INPUTS:
%   model:           COBRA model structure (may contain blocked/inconsistent reactions).
%   coreRxns:        Core reactions requested in the extracted model (numeric
%                    indices, logical mask, or reaction ID cell/string array).
%                    Blocked core reactions are omitted and returned in BlockedCoreRxns.
%
% OPTIONAL INPUTS:
%   tol:             Minimum absolute flux required for consistency (Default: 1e-4).
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
%   seed:            Optional RNG seed for reproducibility (Default: []).
%   useCoupling:     Honor coupling constraints model.C, d, dsense (Default: false).
%   inclusionCutoff: Minimum absolute flux to retain a non-core reaction (Default: derived).
%   solverTol:       LP/MILP feasibility & optimality (and integrality) tolerance,
%                    or 'default' (Default: 1e-9).
%   indicatorGroups: Optional MILP binary indicator groups (Default: {}).
%
% OUTPUTS:
%   Model:           Extracted consistent model (or cell array when nSol > 1).
%   BlockedCoreRxns: Core reactions blocked at tol in the input model.
%   LPs:             Number of LPs solved (or cell array when nSol > 1).
%   selection:       Group-selection struct when MILP binaries are used.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end
if ~isnumeric(tol) || ~isscalar(tol) || ~isfinite(tol) || tol <= 0
    error('spectraCCME:badTol', 'tol must be a positive finite scalar.');
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
    error('spectraCCME:badConsType', ...
        'consType must be ''stoichiometry'' or ''topology'', not ''%s''.', consType);
end

if ~exist('probType', 'var') || isempty(probType)
    probType = 'minNetLP';
end
if ~ismember(probType, {'minNetLP', 'minNetMILP', 'tradeOff', 'optimBiomass'})
    error('spectraCCME:badProbType', ...
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
        error('spectraCCME:badAltSolMethod', ...
            'altSolMethod must be ''coreDirection'' or ''pathwayExclusion'', not ''%s''.', altSolMethod);
    end
    if strcmp(altSolMethod, 'pathwayExclusion') && ~isMilpProb
        error('spectraCCME:pathwayExclusionNeedsMilp', ...
            'altSolMethod ''pathwayExclusion'' requires ''minNetMILP'' or ''tradeOff'', not ''%s''.', probType);
    end
end

if ~exist('weights', 'var') || isempty(weights)
    weights         = ones(n, 1);
    userGaveWeights = false;
else
    if ~isnumeric(weights) || any(~isfinite(weights(:)))
        error('spectraCCME:badWeights', 'weights must be finite numeric values.');
    end
    weights         = double(weights(:));
    userGaveWeights = true;
end
if ~exist('solveTime', 'var') || isempty(solveTime)
    solveTime = 7200;
end
if ~exist('remGene', 'var') || isempty(remGene)
    remGene = 0;
end
if exist('seed', 'var') && ~isempty(seed)
    rng(seed);
end
if ~exist('useCoupling', 'var') || isempty(useCoupling)
    useCoupling = false;
else
    useCoupling = logical(useCoupling);
end
if ~exist('inclusionCutoff', 'var')
    inclusionCutoff = [];
end
if ~exist('solverTol', 'var')
    solverTol = [];
end
if ~exist('indicatorGroups', 'var')
    indicatorGroups = {};
end

% Validate indicatorGroups and weights
[indicatorGroups, milpWeights] = localResolveGroupArgument( ...
    indicatorGroups, weights, userGaveWeights, model.rxns, probType);
if isempty(indicatorGroups) && numel(weights) ~= n
    error('spectraCCME:badWeights', ...
        'weights must have one entry per reaction (%d), got %d.', n, numel(weights));
elseif ~isempty(indicatorGroups) && numel(milpWeights) ~= numel(indicatorGroups)
    error('spectraCCME:groupWeightCount', ...
        ['With indicatorGroups, weights must hold one entry per reaction (%d) ' ...
         'or one entry per group (%d), got %d.'], n, numel(indicatorGroups), numel(weights));
end
if ~strcmp(probType, 'tradeOff') && (any(weights < 0) || any(milpWeights < 0))
    error('spectraCCME:negativeWeights', ...
        'weights must be non-negative for probType ''%s''; use ''tradeOff'' for signed weights.', probType);
end

tolCfg = localResolveTolerances(tol, solverTol, inclusionCutoff);

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
% Stage 1a: 0-LP Dead-End Pruning + Core-Scoped forwardcc/reverse on mPinned
%           (Skipped Automatically When coreRxns Is Empty!)
% =========================================================================
[direction, witnessFlux, trulyBlockedCore, survCols, dirLPs] = localDiscoverConsistentCoreDirections( ...
    model, origModel.rxns, coreMask, tol, steadystate, useCoupling, solverTol, tolCfg, false);

BlockedCoreRxns = origModel.rxns(trulyBlockedCore);
LPs             = dirLPs;
selection       = [];

workModel = model;
shutCols  = (~survCols) | trulyBlockedCore;
if any(shutCols)
    workModel.lb(shutCols)  = 0.0;
    workModel.ub(shutCols)  = 0.0;
    workModel.rev(shutCols) = 0;
end

% Return empty model immediately if core was non-empty and 100% of core is blocked
if any(coreMask) && ~any(direction ~= 0)
    Model = removeRxns(origModel, origModel.rxns);
    if remGene
        Model = removeUnusedGenes(Model);
    end
    if nSol > 1
        Model = repmat({Model}, nSol, 1);
        if strcmp(altSolMethod, 'coreDirection')
            LPs = num2cell([LPs; zeros(nSol - 1, 1)]);
        end
    end
    return;
end

% =========================================================================
% Stage 1b: Anchor-Refined Reversible Active Core Sign Selection
%           (Supports non-negative cost weights & tradeOff negative penalties)
% =========================================================================
activeCorePrimary = (direction ~= 0);
pass0FluxReady    = [];
canFwdWork        = (workModel.ub >= tolCfg.lockTol);
canRevWork        = (workModel.lb <= -tolCfg.lockTol);
twoWayAct         = activeCorePrimary & canFwdWork & canRevWork;
oneWayAct         = activeCorePrimary & (~twoWayAct);

wAnchorRef   = [];
useAnchorRef = false;
if (numel(weights) == n) && any(twoWayAct) && ...
        (any(oneWayAct) || ~localZeroFluxIsFeasible(workModel, useCoupling))
    if ismember(probType, {'minNetLP', 'minNetMILP'})
        wAnchorRef                     = weights;
        wAnchorRef(activeCorePrimary)  = 0.0;
        useAnchorRef                   = true;
    elseif strcmp(probType, 'tradeOff') && any(weights < 0)
        wAnchorRef                     = zeros(n, 1);
        negW                           = (weights < 0);
        wAnchorRef(negW)               = -weights(negW);
        wAnchorRef(activeCorePrimary)  = 0.0;
        useAnchorRef                   = true;
    end
end

if useAnchorRef
    [directionRef, witnessRef, lpsRef, okRef, exactPass0] = localRefineRevCoreViaAnchorPinned( ...
        workModel, activeCorePrimary, twoWayAct, direction, wAnchorRef, tol, ...
        steadystate, useCoupling, inclusionCutoff, solverTol, tolCfg);
    LPs = LPs + lpsRef;
    if okRef
        direction   = directionRef;
        witnessFlux = witnessRef;
        if exactPass0 && strcmp(probType, 'minNetLP')
            pass0FluxReady = witnessRef;
        end
    end
end

% =========================================================================
% Stage 2: Primary Network Extraction (j = 1)
% =========================================================================
[reacInd, extLPs, selected, stat] = localRunSingleExtraction( ...
    workModel, direction, weights, milpWeights, probType, tol, steadystate, ...
    solveTime, witnessFlux, pass0FluxReady, {}, indicatorGroups, ...
    inclusionCutoff, solverTol, useCoupling, tolCfg);
LPs = LPs + extLPs;

Model     = localSubsetModel(origModel, reacInd, remGene);
selection = localReportSelection(indicatorGroups, selected, workModel.rxns, direction);

% =========================================================================
% Alternative Solutions (nSol > 1)
% =========================================================================
if nSol > 1
    Model          = {Model};
    activeCoreMask = coreMask & (~trulyBlockedCore);
    if strcmp(altSolMethod, 'coreDirection')
        % Skip redundant solves if no active core reaction can carry >= tol both ways
        canBothAct = activeCoreMask & (workModel.lb <= -tolCfg.lockTol) & (workModel.ub >= tolCfg.lockTol);
        if ~any(canBothAct)
            Model = repmat(Model, nSol, 1);
            LPs   = num2cell([LPs; zeros(nSol - 1, 1)]);
            return;
        end
        LPs = {LPs};
        for j = 2:nSol
            [directionJ, witnessJ, ~, ~, lpsDirJ] = localDiscoverConsistentCoreDirections( ...
                workModel, origModel.rxns, activeCoreMask, tol, steadystate, ...
                useCoupling, solverTol, tolCfg, true);
            [reacIndJ, lpsExtJ, selected] = localRunSingleExtraction( ...
                workModel, directionJ, weights, milpWeights, probType, tol, ...
                steadystate, solveTime, witnessJ, [], {}, indicatorGroups, ...
                inclusionCutoff, solverTol, useCoupling, tolCfg);
            ModJ  = localSubsetModel(origModel, reacIndJ, remGene);
            Model = [Model; {ModJ}]; %#ok<AGROW>
            LPs   = [LPs; {lpsDirJ + lpsExtJ}]; %#ok<AGROW>
        end
    elseif strcmp(altSolMethod, 'pathwayExclusion')
        prevSols = {};
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
            [reacIndJ, lpsExtJ, selected, stat] = localRunSingleExtraction( ...
                workModel, direction, weights, milpWeights, probType, tol, ...
                steadystate, solveTime, [], [], prevSols, indicatorGroups, ...
                inclusionCutoff, solverTol, useCoupling, tolCfg);
            LPs = LPs + lpsExtJ;
            if stat ~= 1
                break;
            end
            ModJ  = localSubsetModel(origModel, reacIndJ, remGene);
            Model = [Model; {ModJ}]; %#ok<AGROW>
        end
    end
    if ~isempty(selected)
        selection = localReportSelection(indicatorGroups, selected, workModel.rxns, direction);
    end
end
end

% =========================================================================
% Stage 1a: Core-Scoped Iterative forwardcc + reverse on mPinned
% =========================================================================
function [direction, lastFlux, trulyBlockedCore, survCols, nLPs] = localDiscoverConsistentCoreDirections( ...
    model, origRxns, coreMask, tol, steadystate, useCoupling, solverTol, tolCfg, isAltSol)

n                = numel(model.rxns);
direction        = zeros(n, 1);
lastFlux         = [];
trulyBlockedCore = false(n, 1);
nLPs             = 0;

openCols = (abs(model.lb) >= tolCfg.deadEndZeroTol) | (abs(model.ub) >= tolCfg.deadEndZeroTol);
if ~useCoupling
    survCols = localStripDeadEnds(model.S, model.lb, model.ub, openCols, ...
        false(n, 1), steadystate, tolCfg.deadEndZeroTol);
else
    survCols = openCols;
end

% Skip all Stage-1 LPs when coreRxns is empty (e.g. only lb > 0 constraints)
if ~any(coreMask)
    return;
end

% Core reactions whose box bounds cannot reach +/-lockTol are blocked in 0 LPs
canReachTol = (model.ub >= tolCfg.lockTol) | (model.lb <= -tolCfg.lockTol);
if ~isAltSol
    trulyBlockedCore(coreMask & (~survCols | ~canReachTol)) = true;
end

workMod  = model;
deadCols = ~survCols;
if any(deadCols)
    workMod.lb(deadCols)  = 0.0;
    workMod.ub(deadCols)  = 0.0;
    workMod.rev(deadCols) = 0;
end

candCore = coreMask & survCols & canReachTol;
if ~any(candCore)
    return;
end

[direction, lastFlux, remU, loopLPs] = localRunPinnedCCLoop( ...
    workMod, zeros(n, 1), candCore, tol, steadystate, useCoupling, solverTol, tolCfg, isAltSol);
nLPs = nLPs + loopLPs;

if ~any(remU) || isAltSol
    return;
end

% Unpinned BlockedCoreRxns Verification Guard on remaining unsettled core U
[stillBlockedU, unpinnedSigns, unpLPs] = localVerifyBlockedCoreUnpinned( ...
    workMod, remU, tol, steadystate, useCoupling, solverTol, tolCfg);
nLPs = nLPs + unpLPs;

trulyBlockedCore(stillBlockedU) = true;
conflictOrGreedy = remU & (~stillBlockedU);

if any(conflictOrGreedy)
    seedDir  = zeros(n, 1);
    seedDir(conflictOrGreedy) = unpinnedSigns(conflictOrGreedy);
    restCore = (candCore & (~stillBlockedU)) & (~conflictOrGreedy);
    [dirRescue, fluxRescue, ~, rescLPs] = localRunPinnedCCLoop( ...
        workMod, seedDir, restCore, tol, steadystate, useCoupling, solverTol, tolCfg, false);
    nLPs = nLPs + rescLPs;

    if ~isempty(fluxRescue) && sum(dirRescue ~= 0) > sum(direction ~= 0)
        direction = dirRescue;
        lastFlux  = fluxRescue;
    elseif isempty(fluxRescue) && sum(direction ~= 0) == 0
        [dirRescue2, fluxRescue2, ~, rescLPs2] = localRunPinnedCCLoop( ...
            workMod, zeros(n, 1), conflictOrGreedy, tol, steadystate, useCoupling, solverTol, tolCfg, false);
        nLPs = nLPs + rescLPs2;
        if sum(dirRescue2 ~= 0) > sum(direction ~= 0)
            direction = dirRescue2;
            lastFlux  = fluxRescue2;
        end
    end

    droppedMutual = (candCore & (~stillBlockedU)) & (direction == 0);
    if any(droppedMutual)
        sampleIDs = strjoin(origRxns(find(droppedMutual, min(5, sum(droppedMutual)))), ', ');
        warning('spectraCCME:mutualCoreConflict', ...
            ['%d core reaction(s) are individually consistent at tol=%g (so not ' ...
             'returned in BlockedCoreRxns), but could not carry >= tol simultaneously ' ...
             'with the other %d active core reactions under finite model bounds ' ...
             '(e.g., %s). Keeping maximal simultaneously consistent core subset.'], ...
            sum(droppedMutual), tol, sum(direction ~= 0), sampleIDs);
    end
end
end

% -------------------------------------------------------------------------
function [direction, lastFlux, U, nLPs] = localRunPinnedCCLoop( ...
    model, initDir, initU, tol, steadystate, useCoupling, solverTol, tolCfg, randOrder)

n         = numel(model.rxns);
direction = initDir(:);
U         = logical(initU(:));
lastFlux  = [];
nLPs      = 0;
lockTol   = tolCfg.lockTol;
canFwd    = (model.ub >= lockTol);
canRev    = (model.lb <= -lockTol);

prevBatchMask = (direction ~= 0);
prevBatchFlux = [];

while any(U)
    before = sum(U);
    if randOrder
        firstIsFwd = (randi(2) == 1);
    else
        firstIsFwd = true;
    end

    % Reuse iterW across forwardcc and reverse so anti-correlated pairs never cancel
    iterW    = zeros(n, 1);
    iterW(U) = 1.0 + 0.1 * rand(before, 1);

    passOrder = [firstIsFwd, ~firstIsFwd];
    for p = 1:2
        doFwd = passOrder(p);
        if ~any(U)
            break;
        end
        if doFwd && ~any(U & canFwd)
            continue;
        end
        if ~doFwd && ~any(U & canRev)
            continue;
        end

        mPinned = model;
        mPinned.lb(direction == 1)  = max(tol, mPinned.lb(direction == 1));
        mPinned.ub(direction == -1) = min(-tol, mPinned.ub(direction == -1));

        if doFwd
            [f, ~] = forwardcc(mPinned, U & canFwd, tol, steadystate, useCoupling, solverTol, iterW);
        else
            [f, ~] = reverse(mPinned, U & canRev, tol, steadystate, useCoupling, solverTol, iterW);
        end
        nLPs = nLPs + 1;

        if any(isnan(f))
            if ~isempty(prevBatchFlux) && any(prevBatchMask)
                borderline = prevBatchMask & (abs(prevBatchFlux) < tol);
                if any(borderline)
                    direction(borderline) = 0;
                    mPinned = model;
                    mPinned.lb(direction == 1)  = max(tol, mPinned.lb(direction == 1));
                    mPinned.ub(direction == -1) = min(-tol, mPinned.ub(direction == -1));
                    if doFwd
                        [f, ~] = forwardcc(mPinned, U & canFwd, tol, steadystate, useCoupling, solverTol, iterW);
                    else
                        [f, ~] = reverse(mPinned, U & canRev, tol, steadystate, useCoupling, solverTol, iterW);
                    end
                    nLPs = nLPs + 1;
                end
            end
            if any(isnan(f))
                lastFlux = [];
                return;
            end
        end

        lastFlux = f;
        newly    = U & (abs(f) >= lockTol);
        if any(newly)
            direction(newly) = sign(f(newly));
            U(newly)         = false;
            prevBatchMask    = newly;
            prevBatchFlux    = f;
        end
    end

    if sum(U) == before
        break;
    end
end
end

% -------------------------------------------------------------------------
function [stillBlockedU, signsFound, nLPs] = localVerifyBlockedCoreUnpinned( ...
    model, candU, tol, steadystate, useCoupling, solverTol, tolCfg)

n          = numel(model.rxns);
U          = logical(candU(:));
signsFound = zeros(n, 1);
nLPs       = 0;
lockTol    = tolCfg.lockTol;
canFwd     = (model.ub >= lockTol);
canRev     = (model.lb <= -lockTol);

while any(U)
    before   = sum(U);
    iterW    = zeros(n, 1);
    iterW(U) = 1.0 + 0.1 * rand(before, 1);

    if any(U & canFwd)
        [f1, ~] = forwardcc(model, U & canFwd, tol, steadystate, useCoupling, solverTol, iterW);
        nLPs = nLPs + 1;
        if ~any(isnan(f1))
            hit1 = U & (abs(f1) >= lockTol);
            signsFound(hit1) = sign(f1(hit1));
            U(hit1)          = false;
        end
    end

    if any(U & canRev)
        [f2, ~] = reverse(model, U & canRev, tol, steadystate, useCoupling, solverTol, iterW);
        nLPs = nLPs + 1;
        if ~any(isnan(f2))
            hit2 = U & (abs(f2) >= lockTol);
            signsFound(hit2) = sign(f2(hit2));
            U(hit2)          = false;
        end
    end

    if sum(U) == before
        break;
    end
end
stillBlockedU = U;
end

% =========================================================================
% Stage 1b: Anchor-Refined Reversible Active Core Sign Selection
% =========================================================================
function [direction, lastFlux, nLPs, okAll, exactPass0] = localRefineRevCoreViaAnchorPinned( ...
    workModel, activeCoreMask, twoWayAct, prevDirection, wAnchor, tol, ...
    steadystate, useCoupling, inclusionCutoff, solverTol, tolCfg)

n          = numel(workModel.rxns);
exactPass0 = false;

direction = prevDirection(:);
direction(twoWayAct) = 0;

[~, vAnchor, statAnc, ~, nLPs] = minNet(workModel, direction, wAnchor, tol, ...
    steadystate, 'minNetLP', [], [], {}, {}, inclusionCutoff, solverTol, useCoupling, false);
if statAnc ~= 1 || isempty(vAnchor)
    lastFlux = [];
    okAll    = false;
    return;
end
lastFlux = vAnchor;

settledRev = twoWayAct & (abs(vAnchor) >= 0.1 * tol);
direction(settledRev) = sign(vAnchor(settledRev));
unsettledRev = twoWayAct & (direction == 0);

if ~any(unsettledRev)
    if all(abs(vAnchor(activeCoreMask)) >= tol)
        okAll      = true;
        exactPass0 = true;
        return;
    elseif all(abs(vAnchor(activeCoreMask)) >= tolCfg.lockTol)
        okAll = true;
        return;
    end
    direction(twoWayAct & (abs(vAnchor) < tolCfg.lockTol)) = 0;
    unsettledRev = twoWayAct & (direction == 0);
end

[direction, fluxPin, remU, pinLPs] = localRunPinnedCCLoop( ...
    workModel, direction, unsettledRev, tol, steadystate, useCoupling, solverTol, tolCfg, false);
nLPs = nLPs + pinLPs;
if ~isempty(fluxPin)
    lastFlux = fluxPin;
end
okAll = ~any(remU) && (sum(direction ~= 0) == sum(activeCoreMask));
end

% =========================================================================
% Stage 2: Dispatch to Selected Formulation (minNet / tradeOff / optimBiomass)
% =========================================================================
function [reacInd, extraLPs, selected, stat] = localRunSingleExtraction( ...
    model, direction, weights, milpWeights, probType, tol, steadystate, ...
    solveTime, witnessFlux, pass0FluxReady, prevSols, indicatorGroups, ...
    inclusionCutoff, solverTol, useCoupling, tolCfg)

n              = numel(model.rxns);
extraLPs       = 0;
selected       = [];
corePinnedMask = (direction ~= 0);
zeroFeas       = ~any(corePinnedMask) && localZeroFluxIsFeasible(model, useCoupling);

if strcmp(probType, 'minNetLP')
    [reacInd, ~, stat, ~, extraLPs] = minNet(model, direction, weights, tol, ...
        steadystate, 'minNetLP', [], pass0FluxReady, {}, {}, ...
        inclusionCutoff, solverTol, useCoupling, true);
    if stat ~= 1 && ~isempty(witnessFlux)
        reacInd = corePinnedMask | (abs(witnessFlux(:)) >= tolCfg.suppTol);
        [reacInd, cLPs] = localCertifyMilpSupport(model, reacInd, witnessFlux, ...
            corePinnedMask, false, steadystate, useCoupling, solverTol, tolCfg);
        extraLPs = extraLPs + cLPs;
        stat = 1;
    elseif isempty(reacInd)
        reacInd = false(n, 1);
    end

elseif strcmp(probType, 'minNetMILP')
    if zeroFeas && isempty(prevSols) && isempty(indicatorGroups) && all(milpWeights > 0)
        reacInd  = false(n, 1);
        selected = false(n, 1);
        stat     = 1;
        return;
    end
    [reacInd, x, stat, selected, warmLPs] = minNet(model, direction, milpWeights, tol, ...
        steadystate, 'minNetMILP', solveTime, [], prevSols, ...
        indicatorGroups, inclusionCutoff, solverTol, useCoupling);
    extraLPs = extraLPs + warmLPs;
    if stat == 1 || stat == 3
        [reacInd, cLPs] = localCertifyMilpSupport(model, reacInd, x, corePinnedMask, ...
            ~isempty(indicatorGroups), steadystate, useCoupling, solverTol, tolCfg);
        extraLPs = extraLPs + cLPs;
    else
        reacInd = false(n, 1);
    end

elseif strcmp(probType, 'optimBiomass')
    if zeroFeas && (~any(model.c ~= 0) || all(weights(model.c ~= 0) == 0))
        reacInd = false(n, 1);
        stat    = 1;
        return;
    end
    [reacInd, ~, stat, extraLPs] = optimBiomass(model, direction, weights, tol, ...
        steadystate, useCoupling, solverTol, inclusionCutoff);
    if isempty(reacInd)
        reacInd = false(n, 1);
    end

elseif strcmp(probType, 'tradeOff')
    if zeroFeas && isempty(prevSols) && isempty(indicatorGroups) && all(milpWeights < 0)
        reacInd  = false(n, 1);
        selected = false(n, 1);
        stat     = 1;
        return;
    end
    [reacInd, x, stat, selected, warmLPs] = tradeOff(model, direction, milpWeights, tol, ...
        steadystate, solveTime, [], prevSols, indicatorGroups, ...
        inclusionCutoff, solverTol, useCoupling);
    extraLPs = extraLPs + warmLPs;
    if stat == 1 || stat == 3
        [reacInd, cLPs] = localCertifyMilpSupport(model, reacInd, x, corePinnedMask, ...
            ~isempty(indicatorGroups), steadystate, useCoupling, solverTol, tolCfg);
        extraLPs = extraLPs + cLPs;
    else
        reacInd = false(n, 1);
    end
end
end

% =========================================================================
% Support Certification & Validation Utilities for Inconsistent Models
% =========================================================================
function tolCfg = localResolveTolerances(tol, solverTol, inclusionCutoff)
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
        error('spectraCCME:badInclusionCutoff', ...
            'inclusionCutoff must be a positive finite scalar.');
    end
    suppTol = double(inclusionCutoff);
    userSetCutoff = true;
else
    suppTol = min(0.1 * tol, max(effSolverTol, 1e-5 * tol));
    userSetCutoff = false;
end

tolCfg = struct();
tolCfg.tol            = tol;
tolCfg.effSolverTol   = effSolverTol;
tolCfg.suppTol        = suppTol;
tolCfg.deadEndZeroTol = min(effSolverTol, suppTol);
tolCfg.rescueTol      = 1e-2 * min(effSolverTol, suppTol);
tolCfg.primalResTol   = min(0.1 * tol, 1e2 * effSolverTol);
if userSetCutoff
    tolCfg.certUnblockedTol = min(tol, max(suppTol, effSolverTol));
else
    tolCfg.certUnblockedTol = min(tol, max(1e2 * effSolverTol, suppTol));
end
tolCfg.lockTol = max(0.99 * tol, tol - 10.0 * effSolverTol);
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
function active = localStripDeadEnds(S, lb, ub, keepMask, coreMask, steadystate, deadEndZeroTol)
active   = logical(keepMask(:));
coreMask = logical(coreMask(:));
while true
    cols = find(active);
    if isempty(cols), break; end
    Ssub    = S(:, cols);
    rowNnz  = full(sum(Ssub ~= 0, 2));
    actRows = find(rowNnz > 0);
    if isempty(actRows), break; end
    Ssub = Ssub(actRows, :);
    rNnz = rowNnz(actRows);

    [rIdx, cIdx, vals] = find(Ssub);
    lbc = lb(cols);
    ubc = ub(cols);

    canPos = ((vals > 0) & (ubc(cIdx) >= deadEndZeroTol)) | ...
             ((vals < 0) & (lbc(cIdx) <= -deadEndZeroTol));
    canNeg = ((vals > 0) & (lbc(cIdx) <= -deadEndZeroTol)) | ...
             ((vals < 0) & (ubc(cIdx) >= deadEndZeroTol));

    nAct   = numel(actRows);
    posCnt = accumarray(rIdx, double(canPos), [nAct, 1]);
    negCnt = accumarray(rIdx, double(canNeg), [nAct, 1]);

    if steadystate == 1
        deadRows = (rNnz == 1) | (posCnt == 0) | (negCnt == 0);
    else
        deadRows = (posCnt == 0) & (negCnt > 0);
    end
    if ~any(deadRows), break; end
    [~, deadLocalCols] = find(Ssub(deadRows, :));
    deadGlobal  = cols(unique(deadLocalCols));
    deadNonCore = deadGlobal(~coreMask(deadGlobal));
    if isempty(deadNonCore), break; end
    active(deadNonCore) = false;
end
end

% -------------------------------------------------------------------------
function [keepMask, extraLPs] = localCertifyMilpSupport( ...
    model, rawReacInd, x, corePinnedMask, hasGroups, steadystate, useCoupling, solverTol, tolCfg)
extraLPs = 0;
n        = numel(model.rxns);
keepMask = false(n, 1);
if islogical(rawReacInd)
    keepMask = rawReacInd(:);
else
    keepMask(rawReacInd) = true;
end

% Fast 0-LP conic ray-scaling check when no ungoverned group reactions exist
if ~hasGroups && ~isempty(x) && numel(x) >= n
    vSol = reshape(double(x(1:n)), n, 1);
    keepCand = corePinnedMask | (abs(vSol) >= tolCfg.suppTol);
    if ~useCoupling
        keepCand = localStripDeadEnds(model.S, model.lb, model.ub, keepCand, ...
            corePinnedMask, steadystate, tolCfg.deadEndZeroTol);
    end
    cols = find(keepCand);
    if isempty(cols)
        keepMask = keepCand;
        return;
    end
    u  = vSol(cols);
    sv = model.S(:, cols) * u;
    if steadystate == 1
        resInf = full(max(abs(sv)));
    else
        resInf = full(max(max(0.0, -sv)));
    end
    if isempty(resInf), resInf = 0.0; end
    resLimit = 0.5 * tolCfg.primalResTol;
    if resInf <= resLimit
        pos = u > 0;
        neg = u < 0;
        maxScale = max(50.0, tolCfg.tol / max(tolCfg.suppTol, eps));
        alphaUb  = maxScale;
        alphaLb  = maxScale;
        if any(pos), alphaUb = min(model.ub(cols(pos)) ./ u(pos)); end
        if any(neg), alphaLb = min(model.lb(cols(neg)) ./ u(neg)); end
        if resInf > 0
            alphaRes = resLimit / resInf;
        else
            alphaRes = maxScale;
        end
        alpha = max(1.0, min([maxScale, alphaUb, alphaLb, alphaRes]));
        cplOK = true;
        if useCoupling && isfield(model, 'C') && ~isempty(model.C)
            cu    = model.C(:, cols) * (alpha * u);
            d     = reshape(model.d, [], 1);
            ds    = reshape(model.dsense, [], 1);
            violL = any((ds == 'L') & (cu - d > tolCfg.primalResTol));
            violG = any((ds == 'G') & (d - cu > tolCfg.primalResTol));
            violE = any((ds == 'E') & (abs(cu - d) > tolCfg.primalResTol));
            cplOK = ~(violL || violG || violE);
        end
        if cplOK
            uScaled = abs(u) * alpha;
            coreLoc = corePinnedMask(cols);
            if all(uScaled(coreLoc) >= 0.99 * tolCfg.tol) && ...
                    all(uScaled(~coreLoc) >= tolCfg.certUnblockedTol)
                keepMask = keepCand;
                return;
            end
        end
    end
    keepMask = keepCand;
end

cols = find(keepMask);
if isempty(cols), return; end
subMod  = struct();
Ssub    = model.S(:, cols);
actRows = full(sum(Ssub ~= 0, 2)) > 0;
subMod.S    = Ssub(actRows, :);
subMod.lb   = model.lb(cols);
subMod.ub   = model.ub(cols);
subMod.rxns = model.rxns(cols);
subMod.mets = model.mets(actRows);
if useCoupling && isfield(model, 'C') && ~isempty(model.C)
    subMod.C      = model.C(:, cols);
    subMod.d      = model.d;
    subMod.dsense = model.dsense;
end
cType = 'stoichiometry';
if steadystate == 0, cType = 'topology'; end
[consIdx, ccLPs] = spectraCC(subMod, tolCfg.certUnblockedTol, cType, [], [], useCoupling, solverTol);
extraLPs  = extraLPs + ccLPs;
isConsSub = false(numel(cols), 1);
if iscell(consIdx)
    isConsSub = ismember(subMod.rxns, consIdx);
else
    isConsSub(consIdx) = true;
end
keepMask(cols(~isConsSub & ~corePinnedMask(cols))) = false;
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
    error('spectraCCME:groupsWithoutBinaries', ...
        ['indicatorGroups groups the MILP''s binaries, and ''%s'' has ' ...
         'none. Use ''minNetMILP'' or ''tradeOff''.'], probType);
end
if ~iscell(given)
    error('spectraCCME:badIndicatorGroups', ...
        'indicatorGroups must be a cell array, of reaction index vectors or of reaction IDs.');
end
named = all(cellfun(@(g) ischar(g) || (isstring(g) && isscalar(g)), given(:)'));
if named
    ids = cellstr(given(:));
    [found, at] = ismember(ids, rxns);
    if ~all(found)
        missing = ids(~found);
        error('spectraCCME:unknownIndicatorReaction', ...
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
