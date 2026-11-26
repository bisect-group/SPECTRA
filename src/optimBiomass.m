function [reacInd, x, stat, nLPs] = optimBiomass(model, direction, weights, tol, steadyState, useCoupling, solverTol, inclusionCutoff)
% USAGE:
%   [reacInd, x, stat, nLPs] = optimBiomass(model, direction, weights, tol, steadyState, useCoupling, solverTol, inclusionCutoff)
%
% INPUTS:
%   model:           COBRA model structure (model.c ~= 0 identifies biomass/objective reactions).
%   direction:       Vector of length n with entries in {-1, 0, 1} specifying
%                    pinned core reaction directions (0 for free reactions).
%   weights:         Non-negative reaction weights (n x 1): rewards flux on
%                    objective reactions (model.c ~= 0) and penalizes absolute
%                    flux on all other reactions.
%   tol:             Minimum absolute flux required for pinned core reactions.
%   steadyState:     1 for stoichiometric steady state (S*v = 0),
%                    0 for topological accumulation (S*v >= 0).
%
% OPTIONAL INPUTS:
%   useCoupling:     Honor coupling constraints model.C, d, dsense (Default: false).
%   solverTol:       LP feasibility and optimality tolerance, or 'default' (Default: 1e-9).
%   inclusionCutoff: Minimum absolute flux to retain a non-core reaction (Default: derived).
%
% OUTPUTS:
%   reacInd:         Logical mask of reactions retained in the extracted model.
%   x:               Flux solution vector (n x 1).
%   stat:            LP solver status (1 = optimal).
%   nLPs:            Number of LPs solved.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('useCoupling', 'var') || isempty(useCoupling)
    useCoupling = false;
else
    useCoupling = logical(useCoupling);
end
if ~exist('solverTol', 'var')
    solverTol = [];
end
if ~exist('inclusionCutoff', 'var')
    inclusionCutoff = [];
end

[m, n] = size(model.S);
direction = reshape(double(direction), n, 1);

if ~isnumeric(weights) || numel(weights) ~= n || any(~isfinite(weights(:)))
    error('optimBiomass:badWeights', ...
        'weights must hold %d finite numeric entries, got %d.', n, numel(weights));
end
weights = double(weights(:));
if any(weights < 0)
    error('optimBiomass:negativeWeights', ...
        'weights must be non-negative for ''optimBiomass''; use ''tradeOff'' for signed weights.');
end

tolCfg = localResolveTolerances(tol, solverTol, inclusionCutoff);

lb = full(double(model.lb(:)));
ub = full(double(model.ub(:)));
colLb = lb;
colUb = ub;
colLb(direction == 1)  = max(tol, colLb(direction == 1));
colUb(direction == -1) = min(-tol, colUb(direction == -1));

isObj   = (model.c(:) ~= 0);
penMask = ~isObj;

penPos = penMask & (colLb >= 0) & (colUb > 0);
penNeg = penMask & (colUb <= 0) & (colLb < 0);
penRev = find(penMask & (colLb < 0) & (colUb > 0));
nSplit = numel(penRev);

colLb(penRev) = 0.0;
splitLb       = zeros(nSplit, 1);
splitUb       = -lb(penRev);

% Minimize -sum_{b in obj} w_b * v_b + sum_{i notin obj} w_i * |v_i|
colCost          = zeros(n, 1);
colCost(isObj)   = -weights(isObj);
colCost(penPos)  = weights(penPos);
colCost(penNeg)  = -weights(penNeg);
colCost(penRev)  = weights(penRev);
splitCost        = weights(penRev);

if nSplit > 0
    Asplit = [model.S, -model.S(:, penRev)];
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
        Acpl = [sparse(model.C), -sparse(model.C(:, penRev))];
    else
        Acpl = sparse(model.C);
    end
    Asplit   = [Asplit; Acpl];
    rhs      = [rhs; reshape(model.d, [], 1)];
    rowSense = [rowSense; reshape(model.dsense, [], 1)];
end

LP = struct( ...
    'A',      Asplit, ...
    'b',      rhs, ...
    'csense', rowSense, ...
    'lb',     [colLb; splitLb], ...
    'ub',     [colUb; splitUb], ...
    'c',      [colCost; splitCost], ...
    'osense', 1);

sol  = solveCobraLP(LP, tolCfg.lpParams{:});
nLPs = 1;
stat = sol.stat;
if stat ~= 1 || isempty(sol.full)
    x       = [];
    reacInd = [];
    return;
end

vOri = sol.full(1:n);
if nSplit > 0
    vOri(penRev) = vOri(penRev) - sol.full(n+1:end);
end
x = vOri;

[reacInd, certLPs] = localCertifySupport( ...
    model, vOri, direction ~= 0, steadyState, useCoupling, solverTol, tolCfg);
nLPs = nLPs + certLPs;
end

% =========================================================================
% Derived Numerical Tolerances & Support Certification
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
        error('optimBiomass:badInclusionCutoff', ...
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

% -------------------------------------------------------------------------
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
