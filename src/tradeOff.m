function [reacInd, x, stat, selected, nLPs] = tradeOff(model, direction, weights, tol, steadyState, solveTime, x0, prevSols, indicatorGroups, inclusionCutoff, solverTol, useCoupling)
% USAGE:
%   [reacInd, x, stat, selected, nLPs] = tradeOff(model, direction, weights, tol, steadyState, solveTime, x0, prevSols, indicatorGroups, inclusionCutoff, solverTol, useCoupling)
%
% INPUTS:
%   model:           COBRA model structure.
%   direction:       Vector of length n with entries in {-1, 0, 1} specifying
%                    pinned core reaction directions (0 for free reactions).
%   weights:         Signed reaction (or group) evidence weights: any finite
%                    real numbers. Positive weights reward inclusion; negative
%                    weights penalize inclusion; zero is neutral.
%   tol:             Minimum absolute flux required for pinned core reactions
%                    and switched-on tradeOff binaries.
%   steadyState:     1 for stoichiometric steady state (S*v = 0),
%                    0 for topological accumulation (S*v >= 0).
%
% OPTIONAL INPUTS:
%   solveTime:       Time limit in seconds for the MILP solve (Default: 7200).
%   x0:              Initial feasible flux vector for MIP warm start (Default: [];
%                    when empty and prevSols is empty, a sparse LP + unsplit
%                    polishing warm start is computed automatically).
%   prevSols:        Cell array of previously excluded binary solution sets (Default: {}).
%   indicatorGroups: Cell array of reaction index vectors sharing one binary.
%                    When provided, weights must hold one entry per group;
%                    free reactions outside all groups get no binary and are
%                    retained unconditionally (Default: {}).
%   inclusionCutoff: Minimum absolute flux to retain a non-core reaction (Default: tol*1e-7).
%   solverTol:       Solver feasibility, optimality, and integrality tolerance,
%                    or 'default' (Default: 1e-9).
%   useCoupling:     Honor coupling constraints model.C, d, dsense (Default: false).
%
% OUTPUTS:
%   reacInd:         Logical mask of reactions retained in the extracted model.
%   x:               Continuous flux vector from the MILP solution.
%   stat:            MILP solver status (1 = optimal, 3 = time-limit incumbent).
%   selected:        Logical vector over MILP binary groups (z >= 0.5).
%   nLPs:            Number of LPs solved during warm-start generation/polishing.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: minNet, optimBiomass, resolveIndicatorGroups, milpTolerances

if ~exist('solveTime', 'var') || isempty(solveTime)
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
if ~exist('inclusionCutoff', 'var') || isempty(inclusionCutoff)
    inclusionCutoff = tol * 1e-7;
    userSetCutoff   = false;
else
    if ~isnumeric(inclusionCutoff) || ~isscalar(inclusionCutoff) || ...
            ~isfinite(inclusionCutoff) || inclusionCutoff <= 0
        error('tradeOff:badInclusionCutoff', ...
            'inclusionCutoff must be a positive finite scalar.');
    end
    inclusionCutoff = double(inclusionCutoff);
    userSetCutoff   = true;
end
if ~exist('solverTol', 'var')
    solverTol = [];
end
if ~exist('useCoupling', 'var') || isempty(useCoupling)
    useCoupling = false;
else
    useCoupling = logical(useCoupling);
end

nLPs      = 0;
[m, n]    = size(model.S);
direction = reshape(double(direction), n, 1);
dir0      = (direction == 0);
tolParams = milpTolerances(solverTol);

if ~isnumeric(weights) || isempty(weights) || any(~isfinite(weights(:)))
    error('tradeOff:badWeights', 'weights must be finite numeric values (positive, zero, or negative).');
end
weights = double(weights(:));

lb    = full(double(model.lb(:)));
ub    = full(double(model.ub(:)));
isRev = (lb < 0);

[groupOf, nGroups, ungoverned] = resolveIndicatorGroups(indicatorGroups, dir0, n);
governed = (groupOf > 0);
govIdx   = find(governed);
govIrr   = find(governed & (~isRev));
govRev   = find(governed & isRev);
nIrr     = numel(govIrr);
n_rev    = numel(govRev);

if isempty(indicatorGroups)
    if numel(weights) ~= n
        error('tradeOff:badWeights', ...
            'weights must have one entry per reaction (%d), got %d.', n, numel(weights));
    end
    groupWeights = weights(dir0);
else
    if numel(weights) ~= nGroups
        error('tradeOff:groupWeightCount', ...
            ['With indicatorGroups, weights must hold one entry per group, ' ...
             'but got %d weights for %d groups. What a group scores is ' ...
             'looked up by the group, not by its members.'], ...
            numel(weights), nGroups);
    end
    groupWeights = weights(:);
end

corePinnedMask = ~dir0;
if ~any(corePinnedMask) && isempty(prevSols) && isempty(indicatorGroups) && ...
        all(groupWeights < 0) && localZeroFluxIsFeasible(model, useCoupling)
    reacInd  = false(n, 1);
    x        = zeros(n, 1);
    stat     = 1;
    selected = false(nGroups, 1);
    return;
end

warmCfg = localResolveTradeOffWarmConfig(solverTol);
if isempty(x0) && isempty(prevSols)
    negGroupW    = max(0.0, -groupWeights);
    wWarm        = zeros(n, 1);
    runWarmIRLS  = false;
    if any(negGroupW > 0) && ~isempty(govIdx)
        grpSizes      = accumarray(groupOf(govIdx), 1, [nGroups, 1]);
        wWarm(govIdx) = negGroupW(groupOf(govIdx)) ./ max(1, grpSizes(groupOf(govIdx)));
        runWarmIRLS   = true;
    end
    [~, vWarm, statWarm, ~, warmLPs] = minNet(model, direction, wWarm, tol, ...
        steadyState, 'minNetLP', [], [], {}, {}, [], warmCfg.warmSolverArg, ...
        useCoupling, runWarmIRLS);
    nLPs = nLPs + warmLPs;
    if statWarm == 1 && ~isempty(vWarm)
        x0 = vWarm;
    end
end

% Decision variables: [v (n); z (nGroups); a_fwd (n_rev); b_rev (n_rev)]
f = [zeros(n, 1); groupWeights; zeros(2 * n_rev, 1)];

% Stoichiometry or topology constraints
Aeq1 = [model.S, sparse(m, nGroups + 2 * n_rev)];
beq1 = zeros(m, 1);
if steadyState
    csenseeq1 = repmat('E', m, 1);
else
    csenseeq1 = repmat('G', m, 1);
end

pickIrr  = sparse(1:nIrr, govIrr, 1, nIrr, n);
shareIrr = sparse(1:nIrr, groupOf(govIrr), 1, nIrr, nGroups);
pickRev  = sparse(1:n_rev, govRev, 1, n_rev, n);
shareRev = sparse(1:n_rev, groupOf(govRev), 1, n_rev, nGroups);

% Irreversible governed reactions: tol * z <= v <= ub * z
Aineq1      = [-pickIrr, tol * shareIrr, sparse(nIrr, 2 * n_rev)];
bineq1      = zeros(nIrr, 1);
csenseineq1 = repmat('L', nIrr, 1);

Aineq2      = [pickIrr, -spdiags(ub(govIrr), 0, nIrr, nIrr) * shareIrr, sparse(nIrr, 2 * n_rev)];
bineq2      = zeros(nIrr, 1);
csenseineq2 = repmat('L', nIrr, 1);

% Exclusion cuts on previously selected binary indicator sets
Aineq3      = sparse(0, n + nGroups + 2 * n_rev);
bineq3      = zeros(0, 1);
csenseineq3 = repmat('L', 0, 1);
if ~isempty(prevSols)
    nExcl       = numel(prevSols);
    Aineq3      = sparse(nExcl, n + nGroups + 2 * n_rev);
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
        Aineq3(id, n+1:n+nGroups) = maskExcl;
        bineq3(id)                = sum(maskExcl) - 1;
    end
end

% Reversible governed reactions: a_fwd + b_rev = z,
%   tol * a_fwd + lb * b_rev <= v <= ub * a_fwd - tol * b_rev
Aeq2      = [sparse(n_rev, n), -shareRev, speye(n_rev), speye(n_rev)];
beq2      = zeros(n_rev, 1);
csenseeq2 = repmat('E', n_rev, 1);

Aineq4      = [pickRev, sparse(n_rev, nGroups), ...
               -spdiags(tol * ones(n_rev, 1), 0, n_rev, n_rev), ...
               -spdiags(lb(govRev), 0, n_rev, n_rev)];
bineq4      = zeros(n_rev, 1);
csenseineq4 = repmat('G', n_rev, 1);

Aineq5      = [pickRev, sparse(n_rev, nGroups), ...
               -spdiags(ub(govRev), 0, n_rev, n_rev), ...
               spdiags(tol * ones(n_rev, 1), 0, n_rev, n_rev)];
bineq5      = zeros(n_rev, 1);
csenseineq5 = repmat('L', n_rev, 1);

% Bounds (preserving tighter box bounds on pinned core reactions)
lb(direction == 1)  = max(tol, lb(direction == 1));
ub(direction == -1) = min(-tol, ub(direction == -1));
lb = [lb; zeros(nGroups + 2 * n_rev, 1)];
ub = [ub; ones(nGroups + 2 * n_rev, 1)];

[Acpl, bcpl, csensecpl] = couplingRows(model, nGroups + 2 * n_rev, useCoupling);

MILPproblem = struct( ...
    'A',       [Aeq1; Aeq2; Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Acpl], ...
    'b',       [beq1; beq2; bineq1; bineq2; bineq3; bineq4; bineq5; bcpl], ...
    'lb',      lb, ...
    'ub',      ub, ...
    'c',       f, ...
    'osense',  -1, ...
    'vartype', [repmat('C', n, 1); repmat('B', nGroups + 2 * n_rev, 1)], ...
    'csense',  [csenseeq1; csenseeq2; csenseineq1; csenseineq2; csenseineq3; csenseineq4; csenseineq5; csensecpl]);

if ~isempty(x0) && numel(x0) >= n && ~any(isnan(x0(1:n)))
    [x0Full, polLPs] = localBuildPolishedTradeOffStart( ...
        MILPproblem, Aeq1(:, 1:n), beq1, csenseeq1, Acpl(:, 1:n), bcpl, csensecpl, ...
        x0(1:n), groupOf, govIrr, govRev, isRev, nGroups, tol, warmCfg);
    nLPs = nLPs + polLPs;
    if ~isempty(x0Full)
        MILPproblem.x0 = x0Full;
    end
end

solution = solveCobraMILP(MILPproblem, 'timeLimit', solveTime, tolParams{:});
stat     = solution.stat;
if stat == 1 || stat == 3
    x          = solution.cont;
    z          = solution.int;
    selected   = z(1:nGroups) >= 0.5;
    selectedOn = ismember(groupOf, find(selected));
    binMask    = (direction ~= 0) | selectedOn | ungoverned;
    if userSetCutoff
        reacInd = binMask | (abs(x(1:n)) >= inclusionCutoff);
    else
        [reacInd, x, certLPs] = localResolveTradeOffSupportOptionC( ...
            x(1:n), binMask, inclusionCutoff, Aeq1(:, 1:n), beq1, csenseeq1, ...
            Acpl(:, 1:n), bcpl, csensecpl, lb(1:n), ub(1:n), warmCfg);
        nLPs = nLPs + certLPs;
    end
else
    x        = [];
    reacInd  = [];
    selected = [];
end
end

% =========================================================================
% tradeOff Warm-Start & Unsplit Feasibility Polisher
% =========================================================================
function warmCfg = localResolveTradeOffWarmConfig(solverTol)
FLOOR      = 1e-9;
useDefault = false;
if isempty(solverTol)
    effMilpTol    = 1e-8;
    effLpTol      = FLOOR;
    warmSolverArg = FLOOR;
elseif (ischar(solverTol) || (isstring(solverTol) && isscalar(solverTol))) && ...
        strcmpi(solverTol, 'default')
    milpFeas      = getCobraSolverParams('MILP', {'feasTol'}, struct());
    effMilpTol    = max(FLOOR, double(milpFeas));
    effLpTol      = FLOOR;
    warmSolverArg = FLOOR;
    useDefault    = true;
else
    effMilpTol    = max(FLOOR, double(solverTol));
    effLpTol      = effMilpTol;
    warmSolverArg = effLpTol;
end

warmCfg = struct();
warmCfg.effLpTol      = effLpTol;
warmCfg.effMilpTol    = effMilpTol;
warmCfg.warmSolverArg = warmSolverArg;
warmCfg.useDefault    = useDefault;
warmCfg.lpParams      = {'feasTol', effLpTol, 'optTol', effLpTol};
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
function [x0Full, polLPs] = localBuildPolishedTradeOffStart( ...
    MILPproblem, AeqFlux, beq, csenseeq, AcplFlux, bcpl, csensecpl, ...
    rawFlux, groupOf, govIrr, govRev, isRev, nGroups, tol, warmCfg)

x0Full = [];
polLPs = 0;
n      = size(AeqFlux, 2);
colLb  = MILPproblem.lb(1:n);
colUb  = MILPproblem.ub(1:n);
vIn    = reshape(double(rawFlux), n, 1);
nRev   = numel(govRev);

membersOf = cell(nGroups, 1);
for g = 1:nGroups
    membersOf{g} = find(groupOf == g);
end

Apol     = [AeqFlux; AcplFlux];
bpol     = [beq; bcpl];
sensepol = [csenseeq; csensecpl];
cuts     = [0.99 * tol, max(10.0 * warmCfg.effLpTol, 1e-11)];

for cIdx = 1:numel(cuts)
    cut     = cuts(cIdx);
    warm    = false(nGroups, 1);
    revSign = zeros(n, 1);

    for g = 1:nGroups
        mem = membersOf{g};
        if isempty(mem)
            continue;
        end
        forcedG = any((colLb(mem) > 0) | (colUb(mem) < 0));
        wantG   = forcedG || any(abs(vIn(mem)) >= cut);
        if ~wantG
            continue;
        end

        okGroup  = true;
        tmpSigns = zeros(numel(mem), 1);
        for k = 1:numel(mem)
            idx = mem(k);
            if ~isRev(idx)
                if colUb(idx) < tol
                    okGroup = false;
                    break;
                end
            else
                canF = (colUb(idx) >= tol);
                canR = (colLb(idx) <= -tol);
                if vIn(idx) > 0 && canF
                    tmpSigns(k) = 1;
                elseif vIn(idx) < 0 && canR
                    tmpSigns(k) = -1;
                elseif colUb(idx) < 0 && canR
                    tmpSigns(k) = -1;
                elseif canF && ~canR
                    tmpSigns(k) = 1;
                elseif canR && ~canF
                    tmpSigns(k) = -1;
                elseif canF
                    tmpSigns(k) = 1;
                else
                    okGroup = false;
                    break;
                end
            end
        end

        if okGroup
            warm(g) = true;
            for k = 1:numel(mem)
                if isRev(mem(k))
                    revSign(mem(k)) = tmpSigns(k);
                end
            end
        end
    end

    lbPol   = colLb;
    ubPol   = colUb;
    offMask = false(n, 1);

    if ~isempty(govIrr)
        onIrr = warm(groupOf(govIrr));
        irrOn = govIrr(onIrr);
        irrOff = govIrr(~onIrr);
        if ~isempty(irrOn)
            lbPol(irrOn) = max(tol, lbPol(irrOn));
        end
        if ~isempty(irrOff)
            lbPol(irrOff)  = 0.0;
            ubPol(irrOff)  = 0.0;
            offMask(irrOff) = true;
        end
    end

    aWarm = zeros(nRev, 1);
    bWarm = zeros(nRev, 1);
    if nRev > 0
        onRev = warm(groupOf(govRev));
        for k = 1:nRev
            idx = govRev(k);
            if ~onRev(k)
                lbPol(idx)   = 0.0;
                ubPol(idx)   = 0.0;
                offMask(idx) = true;
            elseif revSign(idx) >= 0
                aWarm(k)   = 1.0;
                lbPol(idx) = max(tol, lbPol(idx));
            else
                bWarm(k)   = 1.0;
                ubPol(idx) = min(-tol, ubPol(idx));
            end
        end
    end

    if any(lbPol > ubPol)
        continue;
    end

    % Fast 0-LP check when rawFlux already satisfies tradeOff constraints
    vTry          = vIn;
    vTry(offMask) = 0.0;
    vTry          = min(ubPol, max(lbPol, vTry));
    cand0         = [vTry; double(warm); aWarm; bWarm];
    if localMaxMilpStartViolation(MILPproblem, cand0) <= 0.5 * warmCfg.effMilpTol
        x0Full = cand0;
        return;
    end

    % 1-LP unsplit polish enforcing exact tradeOff bounds (v_off = 0, |v_on| >= tol)
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
        vPol          = min(ubPol, max(lbPol, vPol));
        candPol       = [vPol; double(warm); aWarm; bWarm];
        if localMaxMilpStartViolation(MILPproblem, candPol) <= 0.95 * warmCfg.effMilpTol
            x0Full = candPol;
        elseif strcmpi(getCobraSolver('MILP'), 'gurobi')
            x0Full = [nan(n, 1); double(warm); aWarm; bWarm];
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
function [reacInd, xOut, certLPs] = localResolveTradeOffSupportOptionC( ...
    xRaw, binMask, tradeCutoff, AeqFlux, beq, csenseeq, ...
    AcplFlux, bcpl, csensecpl, colLb, colUb, warmCfg)

n          = numel(binMask);
xOut       = reshape(double(xRaw), n, 1);
certLPs    = 0;
extraTrace = (~binMask) & (abs(xOut) >= tradeCutoff);

if ~any(extraTrace)
    reacInd = binMask;
    return;
end

Aall     = [AeqFlux; AcplFlux];
ball     = [beq; bcpl];
senseall = [csenseeq; csensecpl];
milpTolP = {'feasTol', warmCfg.effMilpTol, 'optTol', warmCfg.effMilpTol};

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
    if (solChk.stat ~= 1 || isempty(solChk.full)) && (warmCfg.effMilpTol > warmCfg.effLpTol)
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

allowMask         = binMask | extraTrace;
lbRes             = colLb;
ubRes             = colUb;
lbRes(~allowMask) = 0.0;
ubRes(~allowMask) = 0.0;
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
    reacInd          = binMask | (extraTrace & (abs(vRes) >= tradeCutoff));
    xOut             = vRes;
else
    reacInd = allowMask;
end
end
