function [reacInd, x, stat, selected] = tradeOff(model, direction, weights, tol, steadyState, solveTime, x0, prevSols, indicatorGroups, inclusionCutoff, solverTol, useCoupling)
% USAGE:
%   [reacInd, x, stat, selected] = tradeOff(model, direction, weights, tol, steadyState, solveTime, x0, prevSols, indicatorGroups, inclusionCutoff, solverTol, useCoupling)
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
%   x0:              Initial feasible flux vector for MIP warm start (Default: []).
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
else
    if ~isnumeric(inclusionCutoff) || ~isscalar(inclusionCutoff) || ...
            ~isfinite(inclusionCutoff) || inclusionCutoff <= 0
        error('tradeOff:badInclusionCutoff', ...
            'inclusionCutoff must be a positive finite scalar.');
    end
    inclusionCutoff = double(inclusionCutoff);
end
if ~exist('solverTol', 'var')
    solverTol = [];
end
if ~exist('useCoupling', 'var') || isempty(useCoupling)
    useCoupling = false;
else
    useCoupling = logical(useCoupling);
end

[m, n]    = size(model.S);
direction = reshape(double(direction), n, 1);
dir0      = (direction == 0);

if ~isnumeric(weights) || isempty(weights) || any(~isfinite(weights(:)))
    error('tradeOff:badWeights', 'weights must be finite numeric values (positive, zero, or negative).');
end
weights = double(weights(:));

lb    = full(double(model.lb(:)));
ub    = full(double(model.ub(:)));
isRev = (lb < 0);

[groupOf, nGroups, ungoverned] = resolveIndicatorGroups(indicatorGroups, dir0, n);
governed = (groupOf > 0);
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
    x0Flux = reshape(double(x0(1:n)), n, 1);
    warm   = false(nGroups, 1);
    if isempty(indicatorGroups)
        on = abs(x0Flux) >= 0.99 * tol;
        warm(groupOf(governed & on)) = true;
        aWarm = double(x0Flux(govRev) >= 0.99 * tol);
        bWarm = double(x0Flux(govRev) <= -0.99 * tol);
        MILPproblem.x0 = [x0Flux; double(warm); aWarm; bWarm];
    end
end

tolParams = milpTolerances(solverTol);
solution  = solveCobraMILP(MILPproblem, 'timeLimit', solveTime, tolParams{:});
stat      = solution.stat;
if stat == 1 || stat == 3
    x          = solution.cont;
    z          = solution.int;
    selected   = z(1:nGroups) >= 0.5;
    selectedOn = ismember(groupOf, find(selected));
    reacInd    = (direction ~= 0) | selectedOn | (abs(x(1:n)) >= inclusionCutoff) | ungoverned;
else
    x        = [];
    reacInd  = [];
    selected = [];
end
end
