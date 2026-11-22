function run_fastCoreWeighted_array(taskId, nTasks, dataRoot)
% run_fastCoreWeighted_array  Parallel SLURM/chunk runner for fastCoreWeighted
% across all 1,479 CCLE cell lines using the exact LocalGini core sets C and
% non-core weights W = 1 - RxnImp saved in Results_BuildDepMapModels.mat.
%
% Usage:
%   run_fastCoreWeighted_array(taskId, nTasks)
%   run_fastCoreWeighted_array(taskId, nTasks, '/path/to/DepMapModels_Recon')

if nargin < 3 || isempty(dataRoot)
    dataRoot = fileparts(mfilename('fullpath'));
end

gurobiHome = getenv('GUROBI_HOME');
if ~isempty(gurobiHome) && exist(fullfile(gurobiHome, 'matlab'), 'dir')
    addpath(fullfile(gurobiHome, 'matlab'));
end

initCobraToolbox(false);
changeCobraSolver('gurobi', 'all');
changeCobraSolverParams('LP', 'feasTol', 1e-8);
changeCobraSolverParams('LP', 'optTol', 1e-8);

chunkDir = fullfile(dataRoot, 'fastcore_weighted_chunks');
modelDir = fullfile(dataRoot, 'fastcore');
if ~exist(chunkDir, 'dir'); mkdir(chunkDir); end
if ~exist(modelDir, 'dir'); mkdir(modelDir); end

S = load(fullfile(dataRoot, 'Results_BuildDepMapModels.mat'), 'model', 'RxnImp', 'context');
model   = S.model;
RxnImp  = S.RxnImp;
context = S.context;

nTotal = size(RxnImp, 2);
nRxns  = numel(model.rxns);

chunkSize = ceil(nTotal / nTasks);
idxStart  = (taskId - 1) * chunkSize + 1;
idxEnd    = min(taskId * chunkSize, nTotal);
if idxStart > nTotal
    return;
end

indices   = idxStart:idxEnd;
nSlice    = numel(indices);
chunkFile = fullfile(chunkDir, sprintf('chunk_%03d.mat', taskId));

if exist(chunkFile, 'file')
    C = load(chunkFile, 'rxnBoolSlice', 'timesSlice', 'doneSlice', 'indices');
    rxnBoolSlice = C.rxnBoolSlice;
    timesSlice   = C.timesSlice;
    doneSlice    = C.doneSlice;
else
    rxnBoolSlice = false(nRxns, nSlice);
    timesSlice   = nan(nSlice, 1);
    doneSlice    = false(nSlice, 1);
end

for k = 1:nSlice
    idx = indices(k);
    ctxName = context{idx};
    if endsWith(ctxName, '.mat')
        outModelFile = fullfile(modelDir, ctxName);
    else
        outModelFile = fullfile(modelDir, [ctxName '.mat']);
    end

    if doneSlice(k) && exist(outModelFile, 'file')
        continue;
    end

    curr_core = RxnImp(:, idx);
    weights = zeros(nRxns, 1);
    weights(curr_core < 1) = 1 - curr_core(curr_core < 1);
    core = find(curr_core >= 1);

    t0 = tic;
    [~, A] = evalc('fastCoreWeighted(core, model, weights, 1e-4);');
    dt = toc(t0);

    rxnBool = false(nRxns, 1);
    rxnBool(A) = true;
    rxnBoolSlice(:, k) = rxnBool;
    timesSlice(k)      = dt;
    doneSlice(k)       = true;

    m = removeRxns(model, model.rxns(~rxnBool)); %#ok<NASGU>
    save(outModelFile, 'm');
    save(chunkFile, 'rxnBoolSlice', 'timesSlice', 'doneSlice', 'indices');
    fprintf('[%d/%d] idx=%d (%s): |C|=%d, |N_fcw|=%d, time=%.2f s\n', ...
        k, nSlice, idx, ctxName, numel(core), numel(A), dt);
end
end
