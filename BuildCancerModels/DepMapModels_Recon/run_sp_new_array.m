function run_sp_new_array(taskId, nTasks)
% run_sp_new_array  Rebuild slice of 1,479 CCLE CSMs using updated SPECTRA (commit 69d2281).
gurobiHome = getenv('GUROBI_HOME');
if ~isempty(gurobiHome) && exist(fullfile(gurobiHome, 'matlab'), 'dir')
    addpath(fullfile(gurobiHome, 'matlab'));
end

% Put NEW spectra_69d2281 FIRST on MATLAB path (never add old spectra-src)
spNewPath = '/scratch/pavan/q_4_4/spectra_69d2281/src';
addpath(spNewPath);
addpath(genpath('/scratch/pavan/cobratoolbox'));
initCobraToolbox(false, 'agent');
changeCobraSolver('gurobi', 'all');
changeCobraSolverParams('LP', 'feasTol', 1e-8);
changeCobraSolverParams('LP', 'optTol', 1e-8);

wSp = which('spectraME');
wMn = which('minNet');
fprintf('Using spectraME: %s\nUsing minNet:    %s\n', wSp, wMn);
assert(startsWith(wSp, spNewPath) && startsWith(wMn, spNewPath), ...
    'Wrong spectraME/minNet on MATLAB path!');

outRoot  = '/scratch/pavan/q_4_4/rerun_69d2281';
chunkDir = fullfile(outRoot, 'sp_chunks');
modelDir = fullfile(outRoot, 'models', 'spectra');
if ~exist(chunkDir, 'dir'); mkdir(chunkDir); end
if ~exist(modelDir, 'dir'); mkdir(modelDir); end

S = load('/scratch/pavan/q_4_4/Results_BuildDepMapModels.mat', 'model', 'RxnImp', 'context');
model   = S.model;
RxnImp  = S.RxnImp;
context = S.context;

nTotal = size(RxnImp, 2);
nRxns  = numel(model.rxns);

% Warm up Gurobi MEX once before timing loop so Task-first model isn't penalized by JIT/MEX load
c1 = RxnImp(:, 1);
w1 = zeros(nRxns, 1); w1(c1 < 1) = 1 - c1(c1 < 1);
rng(1, 'twister');
[~, ~] = spectraME(model, find(c1 >= 1), 1e-4, [], w1);

chunkSize = ceil(nTotal / nTasks);
idxStart  = (taskId - 1) * chunkSize + 1;
idxEnd    = min(taskId * chunkSize, nTotal);
if idxStart > nTotal
    fprintf('Task %d has no indices. Exiting.\n', taskId);
    return;
end

indices   = idxStart:idxEnd;
nSlice    = numel(indices);
chunkFile = fullfile(chunkDir, sprintf('chunk_%03d.mat', taskId));

if exist(chunkFile, 'file')
    C = load(chunkFile, 'rxnBoolSlice', 'timesSlice', 'lpsSlice', 'doneSlice', 'indices');
    rxnBoolSlice = C.rxnBoolSlice;
    timesSlice   = C.timesSlice;
    lpsSlice     = C.lpsSlice;
    doneSlice    = C.doneSlice;
else
    rxnBoolSlice = false(nRxns, nSlice);
    timesSlice   = nan(nSlice, 1);
    lpsSlice     = nan(nSlice, 1);
    doneSlice    = false(nSlice, 1);
end

fprintf('Task %d/%d: running indices %d..%d (%d models)\n', taskId, nTasks, idxStart, idxEnd, nSlice);

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

    rng(idx, 'twister');
    t0 = tic;
    [m, lps] = spectraME(model, core, 1e-4, [], weights);
    dt = toc(t0);

    rxnBool = ismember(model.rxns, m.rxns);
    assert(all(rxnBool(core)), 'Core reactions missing in %s!', ctxName);

    rxnBoolSlice(:, k) = rxnBool;
    timesSlice(k)      = dt;
    lpsSlice(k)        = double(lps);
    doneSlice(k)       = true;

    save(outModelFile, 'm');
    save(chunkFile, 'rxnBoolSlice', 'timesSlice', 'lpsSlice', 'doneSlice', 'indices');
    fprintf('[%d/%d] idx=%d (%s): |C|=%d, |N_sp|=%d, LPs=%d, time=%.2f s\n', ...
        k, nSlice, idx, ctxName, numel(core), numel(m.rxns), lps, dt);
end
fprintf('Task %d/%d FINISHED all %d models.\n', taskId, nTasks, nSlice);
end
