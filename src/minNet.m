function [reacInd, x, stat, selected, nLPs] = minNet(model, direction, weights, tol, steadyState, probType, solveTime, x0, prevSols, indicatorGroups, inclusionCutoff, solverTol, useCoupling, runIRLS)
% USAGE:
%   [reacInd, x, stat, selected, nLPs] = minNet(model, direction, weights, tol, steadyState, probType, solveTime, x0, prevSols, indicatorGroups, inclusionCutoff, solverTol, useCoupling, runIRLS)
%
% INPUTS:
%   model:           COBRA model structure.
%   direction:       Vector of length n with entries in {-1, 0, 1} specifying
%                    pinned core reaction directions (0 for free reactions).
%   weights:         Non-negative reaction (or group) cost weights. Higher
%                    weights penalize inclusion more strongly; zero is free.
%   tol:             Minimum absolute flux required for pinned core reactions.
%   steadyState:     1 for stoichiometric steady state (S*v = 0),
%                    0 for topological accumulation (S*v >= 0).
%
% OPTIONAL INPUTS:
%   probType:        'minNetLP' (default: Split-Variable L1 + warm-started IRLS
%                    with "With-Chain" median-normalized weights) or
%                    'minNetMILP' (exact binary reaction-count minimization).
%   solveTime:       Time limit in seconds for 'minNetMILP' (Default: 7200).
%   x0:              Initial feasible flux vector (optional MIP incumbent
%                    warm start for 'minNetMILP' — when empty and prevSols
%                    is empty, a sparse Split-Variable minNetLP + unsplit
%                    polishing warm start is computed automatically; optional
%                    Pass-0 flux vector when already solved at >= tol for
%                    'minNetLP').
%   prevSols:        Cell array of excluded binary solution sets ('minNetMILP').
%   indicatorGroups: Cell array of reaction index vectors sharing one binary
%                    in 'minNetMILP' (Default: one binary per free reaction).
%   inclusionCutoff: Minimum absolute flux to retain a non-core reaction
%                    (Default: derived from (tol, solverTol) for 'minNetLP',
%                    tol*1e-7 for 'minNetMILP').
%   solverTol:       Solver feasibility/optimality (and integrality) tolerance,
%                    or 'default' (Default: 1e-9).
%   useCoupling:     Honor coupling constraints model.C, d, dsense (Default: false).
%   runIRLS:         Logical flag for 'minNetLP': true (default) runs IRLS
%                    refinement + support certification; false solves a single
%                    Pass-0 Split-Variable minNetLP (used as Stage-1 Anchor).
%
% OUTPUTS:
%   reacInd:         Logical mask of reactions retained in the extracted model.
%   x:               Flux solution vector (n x 1).
%   stat:            Solver optimization status (1 = optimal, 3 = time-limit incumbent).
%   selected:        Logical vector over MILP binary groups ([] for 'minNetLP').
%   nLPs:            Number of LPs solved (including MILP warm-start LPs).
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('probType', 'var') || isempty(probType)
    probType = 'minNetLP';
end
if strcmp(probType, 'minNetMILP') && (~exist('solveTime', 'var') || isempty(solveTime))
    solveTime = 7200;
end
if ~exist('x0', 'var')
    x0 = [];
end
if ~exist('prevSols', 'var') || isempty(prevSols)
    prevSols = {};
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
if ~exist('runIRLS', 'var') || isempty(runIRLS)
    runIRLS = true;
else
    runIRLS = logical(runIRLS);
end

selected = [];
nLPs     = 0;
[m, n]   = size(model.S);
direction = reshape(double(direction), n, 1);
dir0      = (direction == 0);

if ~isnumeric(weights) || isempty(weights) || any(~isfinite(weights(:)))
    error('minNet:badWeights', 'weights must be finite numeric values.');
end
weights = double(weights(:));
if any(weights < 0)
    error('minNet:negativeWeights', ...
        'weights must be non-negative for ''%s''; use ''tradeOff'' for signed weights.', probType);
end

% =========================================================================
% Branch 1: minNetLP (Split-Variable LP + "With-Chain" Weights + IRLS)
% =========================================================================
if strcmp(probType, 'minNetLP')
    if numel(weights) ~= n
        error('minNet:badWeights', ...
            'weights must have one entry per reaction (%d), got %d.', n, numel(weights));
    end

    tolCfg = localResolveTolerances(tol, solverTol, inclusionCutoff);
    corePinnedMask = ~dir0;

    % Skip all LPs when no core reaction is pinned and v = 0 is feasible
    if runIRLS && ~any(corePinnedMask) && localZeroFluxIsFeasible(model, useCoupling)
        reacInd = false(n, 1);
        x       = zeros(n, 1);
        stat    = 1;
        return;
    end

    % Compute "With-Chain" median-normalized non-negative weights (0 LPs)
    effWeights = localPrepareMergedWeights( ...
        model.S, model.lb, model.ub, corePinnedMask, weights, steadyState);

    if runIRLS
        deltas = [1e-3, 1e-5, 1e-6, 1e-8];
    else
        deltas = [];
    end

    [vOri, stat, irlsLPs] = localRunSplitMinNetIRLS( ...
        model, direction, effWeights, steadyState, useCoupling, tolCfg, deltas, x0);
    nLPs = nLPs + irlsLPs;

    if stat ~= 1
        x       = [];
        reacInd = [];
        return;
    end

    x = vOri;
    if ~runIRLS
        reacInd = corePinnedMask | (abs(vOri) >= tolCfg.suppTol);
        return;
    end

    [reacInd, certLPs] = localCertifySupport( ...
        model, vOri, corePinnedMask, steadyState, useCoupling, solverTol, tolCfg);
    nLPs = nLPs + certLPs;

% =========================================================================
% Branch 2: minNetMILP (Exact Binary Reaction-Count Minimization)
% =========================================================================
elseif strcmp(probType, 'minNetMILP')
    if isempty(inclusionCutoff)
        milpCutoff    = tol * 1e-7;
        userSetCutoff = false;
    else
        if ~isnumeric(inclusionCutoff) || ~isscalar(inclusionCutoff) || ...
                ~isfinite(inclusionCutoff) || inclusionCutoff <= 0
            error('minNet:badInclusionCutoff', ...
                'inclusionCutoff must be a positive finite scalar.');
        end
        milpCutoff    = double(inclusionCutoff);
        userSetCutoff = true;
    end
    tolParams = milpTolerances(solverTol);

    [groupOf, nGroups, ungoverned] = resolveIndicatorGroups(indicatorGroups, dir0, n);
    governed = find(groupOf > 0);
    nGov     = numel(governed);

    if isempty(indicatorGroups)
        if numel(weights) ~= n
            error('minNet:badWeights', ...
                'weights must have one entry per reaction (%d), got %d.', n, numel(weights));
        end
        groupWeights = weights(dir0);
    else
        if numel(weights) ~= nGroups
            error('minNet:groupWeightCount', ...
                ['With indicatorGroups, weights must hold one entry per ' ...
                 'group, but got %d weights for %d groups. What a group ' ...
                 'costs to keep is looked up by the group, not by its ' ...
                 'members.'], numel(weights), nGroups);
        end
        groupWeights = weights(:);
    end

    corePinnedMask = ~dir0;
    if ~any(corePinnedMask) && isempty(prevSols) && isempty(indicatorGroups) && ...
            all(groupWeights > 0) && localZeroFluxIsFeasible(model, useCoupling)
        reacInd  = false(n, 1);
        x        = zeros(n, 1);
        stat     = 1;
        selected = false(nGroups, 1);
        return;
    end

    warmTolCfg = localResolveMilpWarmConfig(tol, solverTol);
    if isempty(x0) && isempty(prevSols)
        wWarm = zeros(n, 1);
        if nGov > 0
            grpSizes = accumarray(groupOf(governed), 1, [nGroups, 1]);
            wWarm(governed) = groupWeights(groupOf(governed)) ./ max(1, grpSizes(groupOf(governed)));
        end
        effWarmW = localPrepareMergedWeights( ...
            model.S, model.lb, model.ub, corePinnedMask, wWarm, steadyState);
        [vWarm, statWarm, warmLPs] = localRunSplitMinNetIRLS( ...
            model, direction, effWarmW, steadyState, useCoupling, warmTolCfg, ...
            [1e-3, 1e-5, 1e-6, 1e-8], []);
        nLPs = nLPs + warmLPs;
        if statWarm == 1 && ~isempty(vWarm)
            x0 = vWarm;
        end
    end

    f   = [zeros(n, 1); groupWeights];
    Aeq = [model.S, sparse(m, nGroups)];
    beq = zeros(m, 1);
    if steadyState
        csenseeq = repmat('E', m, 1);
    else
        csenseeq = repmat('G', m, 1);
    end

    pick  = sparse(1:nGov, governed, 1, nGov, n);
    share = sparse(1:nGov, groupOf(governed), 1, nGov, nGroups);

    Aineq1      = [pick, -spdiags(model.lb(governed), 0, nGov, nGov) * share];
    bineq1      = zeros(nGov, 1);
    csenseineq1 = repmat('G', nGov, 1);

    Aineq2      = [pick, -spdiags(model.ub(governed), 0, nGov, nGov) * share];
    bineq2      = zeros(nGov, 1);
    csenseineq2 = repmat('L', nGov, 1);

    Aineq3      = sparse(0, n + nGroups);
    bineq3      = zeros(0, 1);
    csenseineq3 = repmat('L', 0, 1);
    if ~isempty(prevSols)
        nExcl       = numel(prevSols);
        Aineq3      = sparse(nExcl, n + nGroups);
        bineq3      = zeros(nExcl, 1);
        csenseineq3 = repmat('L', nExcl, 1);
        freeList    = find(dir0);
        for id = 1:nExcl
            if isempty(indicatorGroups)
                maskExcl = ismember(freeList, prevSols{id});
            else
                maskExcl = ismember((1:nGroups)', prevSols{id});
            end
            maskExcl = reshape(double(maskExcl), 1, nGroups);
            Aineq3(id, n+1:end) = maskExcl;
            bineq3(id)          = sum(maskExcl) - 1;
        end
    end

    lb = full(double(model.lb(:)));
    ub = full(double(model.ub(:)));
    lb(direction == 1)  = max(tol, lb(direction == 1));
    ub(direction == -1) = min(-tol, ub(direction == -1));
    lb = [lb; zeros(nGroups, 1)];
    ub = [ub; ones(nGroups, 1)];

    [Acpl, bcpl, csensecpl] = couplingRows(model, nGroups, useCoupling);

    MILPproblem = struct( ...
        'A',       [Aeq; Aineq1; Aineq2; Aineq3; Acpl], ...
        'b',       [beq; bineq1; bineq2; bineq3; bcpl], ...
        'lb',      lb, ...
        'ub',      ub, ...
        'c',       f, ...
        'osense',  1, ...
        'vartype', [repmat('C', n, 1); repmat('B', nGroups, 1)], ...
        'csense',  [csenseeq; csenseineq1; csenseineq2; csenseineq3; csensecpl]);

    if ~isempty(x0) && numel(x0) >= n && ~any(isnan(x0(1:n)))
        [x0Full, polLPs] = localBuildPolishedMilpStart( ...
            MILPproblem, Aeq(:, 1:n), beq, csenseeq, Acpl(:, 1:n), bcpl, csensecpl, ...
            x0(1:n), groupOf, governed, nGroups, warmTolCfg);
        nLPs = nLPs + polLPs;
        if ~isempty(x0Full)
            MILPproblem.x0 = x0Full;
        end
    end

    solution = solveCobraMILP(MILPproblem, 'timeLimit', solveTime, tolParams{:});
    stat      = solution.stat;
    if stat == 1 || stat == 3
        x          = solution.cont;
        z          = solution.int;
        selected   = z >= 0.5;
        selectedOn = ismember(groupOf, find(selected));
        binMask    = (direction ~= 0) | selectedOn | ungoverned;
        if userSetCutoff
            reacInd = binMask | (abs(x(1:n)) >= milpCutoff);
        else
            [reacInd, x, certLPs] = localResolveMilpSupportOptionC( ...
                x(1:n), binMask, milpCutoff, Aeq(:, 1:n), beq, csenseeq, ...
                Acpl(:, 1:n), bcpl, csensecpl, lb(1:n), ub(1:n), warmTolCfg);
            nLPs = nLPs + certLPs;
        end
    else
        x        = [];
        reacInd  = [];
        selected = [];
    end
else
    error('minNet:badProbType', ...
        'probType must be ''minNetLP'' or ''minNetMILP'', not ''%s''.', probType);
end
end

% =========================================================================
% Derived Numerical Tolerances
% =========================================================================
function tolCfg = localResolveTolerances(tol, solverTol, inclusionCutoff)
FLOOR    = 1e-9;
lpParams = lpTolerances(solverTol);

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
        error('minNet:badInclusionCutoff', ...
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
tolCfg.lpParams       = lpParams;
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
tolCfg.weightCap = max(1e4, 1.0 / max(effSolverTol, 1e-4 * tol));
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

% =========================================================================
% "With-Chain" Median-Normalized Non-Negative Weights (0 LPs)
% =========================================================================
function effWeights = localPrepareMergedWeights(S, lb, ub, corePinnedMask, rawWeights, steadyState)
[~, n] = size(S);
corePinnedMask = logical(corePinnedMask(:));
openMask    = (abs(lb(:)) > 0) | (abs(ub(:)) > 0);
userW       = double(rawWeights(:));
paidNonCore = openMask & (~corePinnedMask) & (userW > 0);
uPos        = userW(paidNonCore);

if isempty(uPos)
    uMed = 1.0;
else
    uMed = median(uPos);
end

uNorm  = zeros(n, 1);
posAll = openMask & (userW > 0);
uNorm(posAll) = min(1e4, max(1e-4, userW(posAll) / uMed));

if steadyState == 1
    mChain = localDegree2ChainMultipliers(S, openMask, paidNonCore);
else
    mChain = ones(n, 1);
end

effWeights = uNorm .* mChain;
end

% -------------------------------------------------------------------------
function multipliers = localDegree2ChainMultipliers(S, openMask, paidNonCore)
[~, n] = size(S);
Sact = S;
if any(~openMask)
    Sact(:, ~openMask) = 0;
end
rowNnz   = full(sum(Sact ~= 0, 2));
deg2Rows = find(rowNnz == 2);
if isempty(deg2Rows)
    multipliers = ones(n, 1);
    return;
end
St     = Sact';
parent = (1:n)';
for idx = 1:numel(deg2Rows)
    r    = deg2Rows(idx);
    cols = find(St(:, r));
    if numel(cols) == 2
        ri = cols(1);
        while parent(ri) ~= ri
            parent(ri) = parent(parent(ri));
            ri = parent(ri);
        end
        rj = cols(2);
        while parent(rj) ~= rj
            parent(rj) = parent(parent(rj));
            rj = parent(rj);
        end
        if ri ~= rj
            parent(ri) = rj;
        end
    end
end
roots = zeros(n, 1);
for i = 1:n
    r = i;
    while parent(r) ~= r
        r = parent(r);
    end
    curr = i;
    while parent(curr) ~= r
        nxt = parent(curr);
        parent(curr) = r;
        curr = nxt;
    end
    roots(i) = r;
end

compPaidCounts = accumarray(roots(paidNonCore), 1, [n, 1]);
multipliers    = ones(n, 1);
multipliers(paidNonCore) = max(1, compPaidCounts(roots(paidNonCore)));
end

% =========================================================================
% Single-Loop Warm-Started Split-Variable minNetLP IRLS
% =========================================================================
function [vOri, statOut, nLPs] = localRunSplitMinNetIRLS( ...
    model, direction, effWeights, steadyState, useCoupling, tolCfg, deltas, pass0FluxReady)

[m, n]   = size(model.S);
tol      = tolCfg.tol;
Mcap     = tolCfg.weightCap;
suppTol  = tolCfg.suppTol;
resTol   = tolCfg.primalResTol;
freeMask = (direction == 0);

lb = full(double(model.lb(:)));
ub = full(double(model.ub(:)));
w        = effWeights(:);
freePos  = find(freeMask & (lb >= 0) & (ub > 0));
freeNeg  = find(freeMask & (ub <= 0) & (lb < 0));
freeRev  = find(freeMask & (lb < 0) & (ub > 0) & (w > 0));
nSplit   = numel(freeRev);
openFree = freeMask & ((lb < 0) | (ub > 0));
paidMask = openFree & (w > 0);

colLb = lb;
colUb = ub;
colLb(direction == 1)  = max(tol, colLb(direction == 1));
colUb(direction == -1) = min(-tol, colUb(direction == -1));
colLb(freeRev) = 0.0;

splitLb = zeros(nSplit, 1);
splitUb = -lb(freeRev);

if nSplit > 0
    Asplit = [model.S, -model.S(:, freeRev)];
else
    Asplit = model.S;
end
if steadyState == 1
    rowSense = repmat('E', m, 1);
else
    rowSense = repmat('G', m, 1);
end
rhs = zeros(m, 1);

if useCoupling && isfield(model, 'C') && ~isempty(model.C)
    if nSplit > 0
        Acpl = [sparse(model.C), -sparse(model.C(:, freeRev))];
    else
        Acpl = sparse(model.C);
    end
    Asplit   = [Asplit; Acpl];
    rhs      = [rhs; reshape(model.d, [], 1)];
    rowSense = [rowSense; reshape(model.dsense, [], 1)];
end

c0 = zeros(n + nSplit, 1);
c0(freePos) = w(freePos);
c0(freeNeg) = -w(freeNeg);
c0(freeRev) = w(freeRev);
if nSplit > 0
    c0(n+1:end) = w(freeRev);
end

LP = struct( ...
    'A',      Asplit, ...
    'b',      rhs, ...
    'csense', rowSense, ...
    'lb',     [colLb; splitLb], ...
    'ub',     [colUb; splitUb], ...
    'c',      c0, ...
    'osense', 1);

nLPs      = 0;
statOut   = 1;
lastBasis = [];

if isempty(pass0FluxReady) || numel(pass0FluxReady) < n || any(isnan(pass0FluxReady(1:n)))
    sol0 = solveCobraLP(LP, tolCfg.lpParams{:});
    nLPs = nLPs + 1;
    if sol0.stat ~= 1 || isempty(sol0.full)
        statOut = sol0.stat;
        if isempty(statOut) || statOut == 1
            statOut = 0;
        end
        vOri = zeros(n, 1);
        return;
    end
    vOri = sol0.full(1:n);
    if nSplit > 0
        vOri(freeRev) = vOri(freeRev) - sol0.full(n+1:end);
    end
    if isfield(sol0, 'basis') && ~isempty(sol0.basis)
        lastBasis = sol0.basis;
    end
else
    vOri = reshape(double(pass0FluxReady(1:n)), n, 1);
end

if isempty(deltas) || ~any(paidMask)
    return;
end

prevSupport = sum(abs(vOri(paidMask)) >= suppTol);
for t = 1:numel(deltas)
    delta = deltas(t);
    cStep = tol * w .* min(Mcap, 1.0 ./ (abs(vOri) + delta));

    cVec = zeros(n + nSplit, 1);
    cVec(freePos) = cStep(freePos);
    cVec(freeNeg) = -cStep(freeNeg);
    cVec(freeRev) = cStep(freeRev);
    if nSplit > 0
        cVec(n+1:end) = cStep(freeRev);
    end
    LP.c = cVec;
    if ~isempty(lastBasis)
        LP.basis = lastBasis;
    end

    solT = solveCobraLP(LP, tolCfg.lpParams{:});
    nLPs = nLPs + 1;
    if solT.stat == 1 && ~isempty(solT.full)
        vCand = solT.full(1:n);
        if nSplit > 0
            vCand(freeRev) = vCand(freeRev) - solT.full(n+1:end);
        end
        if localCheckRowResidual(model.S, vCand, steadyState) <= resTol
            vOri = vCand;
            if isfield(solT, 'basis') && ~isempty(solT.basis)
                lastBasis = solT.basis;
            end
        else
            break;
        end
    else
        break;
    end

    curSupport = sum(abs(vOri(paidMask)) >= suppTol);
    if curSupport == prevSupport && delta <= 0.1 * tol
        break;
    end
    prevSupport = curSupport;
end
end

% -------------------------------------------------------------------------
function resInf = localCheckRowResidual(S, v, steadyState)
sv = S * v;
if steadyState == 1
    resInf = full(max(abs(sv)));
else
    resInf = full(max(max(0.0, -sv)));
end
if isempty(resInf)
    resInf = 0.0;
end
end

% =========================================================================
% Tolerance-Derived Support Certification
% =========================================================================
function active = localStripDeadEnds(S, lb, ub, keepMask, coreMask, steadyState, deadEndZeroTol)
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

    if steadyState == 1
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
function [keepMask, extraLPs] = localCertifySupport( ...
    model, vSol, coreMask, steadyState, useCoupling, solverTol, tolCfg)

S  = model.S;
lb = full(double(model.lb(:)));
ub = full(double(model.ub(:)));
extraLPs = 0;

suppTol    = tolCfg.suppTol;
rescueTol  = tolCfg.rescueTol;
resLimit   = 0.5 * tolCfg.primalResTol;
unblockTol = tolCfg.certUnblockedTol;
zeroTol    = tolCfg.deadEndZeroTol;
tol        = tolCfg.tol;

keepMask = coreMask | (abs(vSol(:)) >= suppTol);
if ~useCoupling
    keepMask = localStripDeadEnds(S, lb, ub, keepMask, coreMask, steadyState, zeroTol);
end
cols = find(keepMask);
if isempty(cols), return; end

u      = vSol(cols);
resInf = localCheckRowResidual(S(:, cols), u, steadyState);
if resInf > resLimit
    keepMask = coreMask | (abs(vSol(:)) >= rescueTol);
    if ~useCoupling
        keepMask = localStripDeadEnds(S, lb, ub, keepMask, coreMask, steadyState, zeroTol);
    end
    cols   = find(keepMask);
    u      = vSol(cols);
    resInf = localCheckRowResidual(S(:, cols), u, steadyState);
end

% Zero-LP Conic Ray-Scaling Certificate (pos = u > 0, neg = u < 0)
if resInf <= resLimit
    pos = u > 0;
    neg = u < 0;
    maxScale = max(50.0, tol / max(suppTol, eps));
    alphaUb  = maxScale;
    alphaLb  = maxScale;
    if any(pos), alphaUb = min(ub(cols(pos)) ./ u(pos)); end
    if any(neg), alphaLb = min(lb(cols(neg)) ./ u(neg)); end
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
        coreLoc = coreMask(cols);
        if all(uScaled(coreLoc) >= 0.99 * tol) && all(uScaled(~coreLoc) >= unblockTol)
            return;
        end
    end
end

subMod  = struct();
Ssub    = S(:, cols);
actRows = full(sum(Ssub ~= 0, 2)) > 0;
subMod.S    = Ssub(actRows, :);
subMod.lb   = lb(cols);
subMod.ub   = ub(cols);
subMod.rxns = model.rxns(cols);
subMod.mets = model.mets(actRows);
if useCoupling && isfield(model, 'C') && ~isempty(model.C)
    subMod.C      = model.C(:, cols);
    subMod.d      = model.d;
    subMod.dsense = model.dsense;
end

if steadyState == 1
    cType = 'stoichiometry';
else
    cType = 'topology';
end
[consIdx, ccLPs] = spectraCC(subMod, unblockTol, cType, [], [], useCoupling, solverTol);
extraLPs  = extraLPs + ccLPs;
isConsSub = false(numel(cols), 1);
if iscell(consIdx)
    isConsSub = ismember(subMod.rxns, consIdx);
else
    isConsSub(consIdx) = true;
end
deadGlobal = cols(~isConsSub & ~coreMask(cols));
keepMask(deadGlobal) = false;
end

% =========================================================================
% MILP Warm-Start & Unsplit Feasibility Polisher
% =========================================================================
function warmCfg = localResolveMilpWarmConfig(tol, solverTol)
FLOOR      = 1e-9;
useDefault = false;
if isempty(solverTol)
    effMilpTol = 1e-8;
    effLpTol   = FLOOR;
elseif (ischar(solverTol) || (isstring(solverTol) && isscalar(solverTol))) && ...
        strcmpi(solverTol, 'default')
    milpFeas   = getCobraSolverParams('MILP', {'feasTol'}, struct());
    effMilpTol = max(FLOOR, double(milpFeas));
    effLpTol   = FLOOR;
    useDefault = true;
else
    effMilpTol = max(FLOOR, double(solverTol));
    effLpTol   = effMilpTol;
end

suppTol = min(0.1 * tol, max(effLpTol, 1e-5 * tol));
warmCfg = struct();
warmCfg.tol          = tol;
warmCfg.lpParams     = {'feasTol', effLpTol, 'optTol', effLpTol};
warmCfg.effSolverTol = effLpTol;
warmCfg.effMilpTol   = effMilpTol;
warmCfg.useDefault   = useDefault;
warmCfg.suppTol      = suppTol;
warmCfg.primalResTol = min(0.1 * tol, 1e2 * effLpTol);
warmCfg.weightCap    = max(1e4, 1.0 / max(effLpTol, 1e-4 * tol));
end

% -------------------------------------------------------------------------
function [x0Full, polLPs] = localBuildPolishedMilpStart( ...
    MILPproblem, AeqFlux, beq, csenseeq, AcplFlux, bcpl, csensecpl, ...
    rawFlux, groupOf, governed, nGroups, warmCfg)

x0Full = [];
polLPs = 0;
n      = size(AeqFlux, 2);
colLb  = MILPproblem.lb(1:n);
colUb  = MILPproblem.ub(1:n);
vIn    = reshape(double(rawFlux), n, 1);
nGov   = numel(governed);

forcedGov = false(nGov, 1);
if nGov > 0
    forcedGov = (colLb(governed) > 0) | (colUb(governed) < 0);
end

Apol     = [AeqFlux; AcplFlux];
bpol     = [beq; bcpl];
sensepol = [csenseeq; csensecpl];
cuts     = [max(10.0 * warmCfg.effSolverTol, 1e-11), 0.1 * warmCfg.effSolverTol];

for cIdx = 1:numel(cuts)
    warm = false(nGroups, 1);
    if nGov > 0
        onGov = forcedGov | (abs(vIn(governed)) >= cuts(cIdx));
        if any(onGov)
            warm(groupOf(governed(onGov))) = true;
        end
    end

    offMask = false(n, 1);
    if nGov > 0
        offMask(governed) = ~warm(groupOf(governed));
    end

    % Fast 0-LP check when rawFlux already satisfies MILP tolerances cleanly
    vTry          = vIn;
    vTry(offMask) = 0.0;
    vTry          = min(colUb, max(colLb, vTry));
    cand0         = [vTry; double(warm)];
    if localMaxMilpStartViolation(MILPproblem, cand0) <= 0.5 * warmCfg.effMilpTol
        x0Full = cand0;
        return;
    end

    % 1-LP unsplit polish (no v^+ - v^- cancellation; off binaries fixed to 0)
    lbPol          = colLb;
    ubPol          = colUb;
    lbPol(offMask) = 0.0;
    ubPol(offMask) = 0.0;
    LPpol = struct( ...
        'A',      Apol, ...
        'b',      bpol, ...
        'csense', sensepol, ...
        'lb',     lbPol, ...
        'ub',     ubPol, ...
        'c',      zeros(n, 1), ...
        'osense', 1);

    solPol = solveCobraLP(LPpol, warmCfg.lpParams{:});
    polLPs = polLPs + 1;
    if (solPol.stat ~= 1 || isempty(solPol.full)) && warmCfg.useDefault
        solPol = solveCobraLP(LPpol);
        polLPs = polLPs + 1;
    end

    if solPol.stat == 1 && ~isempty(solPol.full)
        vPol          = solPol.full(1:n);
        vPol(offMask) = 0.0;
        vPol          = min(colUb, max(colLb, vPol));
        candPol       = [vPol; double(warm)];
        if localMaxMilpStartViolation(MILPproblem, candPol) <= 0.95 * warmCfg.effMilpTol
            x0Full = candPol;
        elseif strcmpi(getCobraSolver('MILP'), 'gurobi')
            x0Full = [nan(n, 1); double(warm)];
        else
            x0Full = candPol;
        end
        return;
    end
end
end

% -------------------------------------------------------------------------
function maxViol = localMaxMilpStartViolation(MILPproblem, xVec)
bndViol = max([0.0; MILPproblem.lb - xVec; xVec - MILPproblem.ub]);
Ax      = MILPproblem.A * xVec;
rhs     = MILPproblem.b;
cs      = MILPproblem.csense;
isE     = (cs == 'E');
isG     = (cs == 'G');
isL     = (cs == 'L');
eqViol  = 0.0;
gViol   = 0.0;
lViol   = 0.0;
if any(isE), eqViol = full(max(abs(Ax(isE) - rhs(isE)))); end
if any(isG), gViol  = full(max(max(0.0, rhs(isG) - Ax(isG)))); end
if any(isL), lViol  = full(max(max(0.0, Ax(isL) - rhs(isL)))); end
maxViol = max([bndViol, eqViol, gViol, lViol]);
end

% -------------------------------------------------------------------------
function [reacInd, xOut, certLPs] = localResolveMilpSupportOptionC( ...
    xRaw, binMask, milpCutoff, AeqFlux, beq, csenseeq, ...
    AcplFlux, bcpl, csensecpl, colLb, colUb, warmCfg)

n          = numel(binMask);
xOut       = reshape(double(xRaw), n, 1);
certLPs    = 0;
extraTrace = (~binMask) & (abs(xOut) >= milpCutoff);

% Step 0: Fast exit when no z_g = 0 reaction carries flux >= milpCutoff
if ~any(extraTrace)
    reacInd = binMask;
    return;
end

Aall     = [AeqFlux; AcplFlux];
ball     = [beq; bcpl];
senseall = [csenseeq; csensecpl];
milpTolP = {'feasTol', warmCfg.effMilpTol, 'optTol', warmCfg.effMilpTol};

% Step 1: 1-LP feasibility check on binMask alone
lbChk           = colLb;
ubChk           = colUb;
lbChk(~binMask) = 0.0;
ubChk(~binMask) = 0.0;
if all(lbChk <= ubChk)
    lpChk = struct( ...
        'A',      Aall, ...
        'b',      ball, ...
        'csense', senseall, ...
        'lb',     lbChk, ...
        'ub',     ubChk, ...
        'c',      zeros(n, 1), ...
        'osense', 1);
    solChk  = solveCobraLP(lpChk, warmCfg.lpParams{:});
    certLPs = certLPs + 1;
    if (solChk.stat ~= 1 || isempty(solChk.full)) && (warmCfg.effMilpTol > warmCfg.effSolverTol)
        solChk  = solveCobraLP(lpChk, milpTolP{:});
        certLPs = certLPs + 1;
    end
    if solChk.stat == 1 && ~isempty(solChk.full)
        vChk           = solChk.full(1:n);
        vChk(~binMask) = 0.0;
        reacInd        = binMask;
        xOut           = vChk;
        return;
    end
end

% Step 2: 1-LP minimal split-variable L1 rescue over extraTrace (e.g. trace cofactors)
allowMask           = binMask | extraTrace;
lbRes               = colLb;
ubRes               = colUb;
lbRes(~allowMask)   = 0.0;
ubRes(~allowMask)   = 0.0;
if any(lbRes > ubRes)
    reacInd = allowMask;
    return;
end

tracePos = find(extraTrace & (lbRes >= 0) & (ubRes > 0));
traceNeg = find(extraTrace & (ubRes <= 0) & (lbRes < 0));
traceRev = find(extraTrace & (lbRes < 0) & (ubRes > 0));
nSplit   = numel(traceRev);

splitLb         = zeros(nSplit, 1);
splitUb         = -colLb(traceRev);
lbRes(traceRev) = 0.0;

if nSplit > 0
    Ares = [Aall, -Aall(:, traceRev)];
else
    Ares = Aall;
end

cRes           = zeros(n + nSplit, 1);
cRes(tracePos) = 1.0;
cRes(traceNeg) = -1.0;
cRes(traceRev) = 1.0;
if nSplit > 0
    cRes(n+1:end) = 1.0;
end

lpRes = struct( ...
    'A',      Ares, ...
    'b',      ball, ...
    'csense', senseall, ...
    'lb',     [lbRes; splitLb], ...
    'ub',     [ubRes; splitUb], ...
    'c',      cRes, ...
    'osense', 1);

solRes  = solveCobraLP(lpRes, milpTolP{:});
certLPs = certLPs + 1;
if solRes.stat == 1 && ~isempty(solRes.full)
    vRes = solRes.full(1:n);
    if nSplit > 0
        vRes(traceRev) = vRes(traceRev) - solRes.full(n+1:end);
    end
    vRes(~allowMask) = 0.0;
    reacInd          = binMask | (extraTrace & (abs(vRes) >= milpCutoff));
    xOut             = vRes;
else
    reacInd = allowMask;
end
end
