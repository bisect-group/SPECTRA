function [Model, LPS] = spectraME2(model, core, tol, weights, nSol, remGene, solverTol)
% USAGE:
%   [Model, LPS] = spectraME2(model, core, tol, weights, nSol, remGene, solverTol)
%
% INPUTS:
%   model:     Consistent COBRA model structure (run spectraCC first).
%   core:      Core reactions required in the extracted model (indices,
%              logical mask, or reaction IDs).
%
% OPTIONAL INPUTS:
%   tol:       Minimum absolute flux required for core reactions (Default: 1e-4).
%   weights:   Non-negative weights for non-core reactions (Default: ones(n, 1)).
%   nSol:      Number of alternative solutions required (Default: 1).
%   remGene:   Remove unused genes from extracted model (Default: 0).
%   solverTol: LP feasibility and optimality tolerance, or 'default' (Default: 1e-9).
%
% OUTPUTS:
%   Model:     Extracted consistent thermodynamic model (or cell array if nSol > 1).
%   LPS:       Total number of LPs solved.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end
tol = double(tol);

[~, n] = size(model.S);
if ~exist('weights', 'var') || isempty(weights)
    weights = ones(n, 1);
else
    if ~isnumeric(weights) || numel(weights) ~= n || any(~isfinite(weights(:)))
        error('spectraME2:badWeights', 'weights must hold %d finite numeric values.', n);
    end
    weights = double(weights(:));
    if any(weights < 0)
        error('spectraME2:negativeWeights', 'weights must be non-negative for spectraME2.');
    end
end
if ~exist('nSol', 'var') || isempty(nSol)
    nSol = 1;
end
if ~exist('remGene', 'var') || isempty(remGene)
    remGene = 0;
end
if ~exist('solverTol', 'var')
    solverTol = [];
end

lpTolerances(solverTol);
if isempty(solverTol)
    effSolverTol = 1e-9;
elseif (ischar(solverTol) || (isstring(solverTol) && isscalar(solverTol))) && ...
        strcmpi(solverTol, 'default')
    effSolverTol = double(getCobraSolverParams('LP', 'feasTol'));
else
    effSolverTol = max(1e-9, double(solverTol));
end
lockTol = max(0.99 * tol, tol - 10.0 * effSolverTol);

origModel = model;
if isempty(core)
    coreMask = false(n, 1);
elseif iscell(core) || isstring(core) || ischar(core)
    if ischar(core), core = {core}; end
    coreMask = ismember(model.rxns, cellstr(core(:)));
elseif islogical(core) && numel(core) == n
    coreMask = core(:);
else
    coreMask = ismember((1:n)', core(:));
end

IrR    = model.ub <= 0;
tempUb = model.ub(IrR);
model.S(:, IrR) = -model.S(:, IrR);
if isfield(model, 'C') && ~isempty(model.C)
    model.C(:, IrR) = -model.C(:, IrR);
end
model.ub(IrR) = -model.lb(IrR);
model.lb(IrR) = -tempUb;
rev = (model.lb < 0);
model.rev = double(rev);

% Pure thermodynamic cone bounds
model.lb(~rev) = 0;
model.lb(rev)  = -Inf;
model.ub(:)    = Inf;

% Preserve original COBRA reaction orientation on the thermodynamic model
thermoOrigModel = origModel;
thermoOrigModel.lb(rev)          = -Inf;
thermoOrigModel.ub(rev)          = Inf;
thermoOrigModel.lb(~rev & ~IrR)  = 0;
thermoOrigModel.ub(~rev & ~IrR)  = Inf;
thermoOrigModel.lb(IrR)          = -Inf;
thermoOrigModel.ub(IrR)          = 0;
thermoOrigModel.rev              = double(rev);

[direction, pass0Ready, dirLPs] = localDiscoverThermoDirections( ...
    model, coreMask, weights, tol, solverTol, lockTol, false);
[reacInd, ~, ~, ~, extLPs] = minNet( ...
    model, direction, weights, tol, 1, 'minNetLP', [], pass0Ready, {}, {}, [], solverTol, false, true);
LPS = dirLPs + extLPs;

Model = removeRxns(thermoOrigModel, thermoOrigModel.rxns(~reacInd));
if remGene
    Model = removeUnusedGenes(Model);
end

if nSol > 1
    Model = {Model};
    if ~any(coreMask & rev)
        Model = repmat(Model, nSol, 1);
        return;
    end
    for j = 2:nSol
        [dirJ, ~, dirLPsJ] = localDiscoverThermoDirections( ...
            model, coreMask, weights, tol, solverTol, lockTol, true);
        [reacIndJ, ~, ~, ~, extLPsJ] = minNet( ...
            model, dirJ, weights, tol, 1, 'minNetLP', [], [], {}, {}, [], solverTol, false, true);
        LPS  = LPS + dirLPsJ + extLPsJ;
        ModJ = removeRxns(thermoOrigModel, thermoOrigModel.rxns(~reacIndJ));
        if remGene
            ModJ = removeUnusedGenes(ModJ);
        end
        Model = [Model; {ModJ}]; %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function [direction, pass0Ready, nLPs] = localDiscoverThermoDirections( ...
    model, coreMask, weights, tol, solverTol, lockTol, isAltSol)
n          = numel(model.rxns);
direction  = zeros(n, 1);
pass0Ready = [];
nLPs       = 0;
if ~any(coreMask)
    return;
end
isRev   = (model.lb < 0);
irrCore = coreMask & (~isRev);
revCore = coreMask & isRev;
direction(irrCore) = 1;
if ~any(revCore)
    return;
end

if ~isAltSol && any(irrCore)
    wAnc = weights;
    wAnc(coreMask) = 0;
    [~, vAnc, statAnc, ~, ancLPs] = minNet( ...
        model, direction, wAnc, tol, 1, 'minNetLP', [], [], {}, {}, [], solverTol, false, false);
    nLPs = nLPs + ancLPs;
    if statAnc == 1 && ~isempty(vAnc)
        settled = revCore & (abs(vAnc) >= 0.1 * tol);
        direction(settled) = sign(vAnc(settled));
        if all(direction(coreMask) ~= 0) && all(abs(vAnc(coreMask)) >= tol)
            pass0Ready = vAnc;
            return;
        elseif all(direction(coreMask) ~= 0) && all(abs(vAnc(coreMask)) >= lockTol)
            return;
        end
        direction(revCore & (abs(vAnc) < lockTol)) = 0;
    end
end

unsettled = revCore & (direction == 0);
stalls    = 0;
while any(unsettled)
    before = sum(unsettled);
    if isAltSol
        firstIsFwd = (randi(2) == 1);
    else
        firstIsFwd = true;
    end
    iterW = zeros(n, 1);
    iterW(unsettled) = 1.0 + 0.1 * rand(before, 1);
    passOrder = [firstIsFwd, ~firstIsFwd];
    for p = 1:2
        if ~any(unsettled), break; end
        mPinned = model;
        mPinned.lb(direction == 1)  = max(tol, mPinned.lb(direction == 1));
        mPinned.ub(direction == -1) = min(-tol, mPinned.ub(direction == -1));
        nLPs = nLPs + 1;
        if passOrder(p)
            [f, ~] = forward(mPinned, unsettled, tol, 1, false, solverTol, iterW);
        else
            [f, ~] = reverse(mPinned, unsettled, tol, 1, false, solverTol, iterW);
        end
        if any(isnan(f))
            error('spectraME2:infeasibleCore', ...
                'No flux distribution lets all %d core reactions carry >= tol=%g.', sum(coreMask), tol);
        end
        newly = unsettled & (abs(f) >= lockTol);
        direction(newly) = sign(f(newly));
        unsettled(newly) = false;
    end
    stalls = assertCoreProgress(unsettled, before, stalls, tol);
end
end
