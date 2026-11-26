function run_timing_1479_array(taskId, nTasks)
% run_timing_1479_array  Apples-to-apples extraction runtime benchmark across
% all 1,479 CCLE CSMs on identical fir hardware & thread allocation (4 CPUs,
% default multi-threaded MATLAB + Gurobi, zero disk I/O inside tic/toc):
%   1. New SPECTRA_ME (commit 69d2281)
%   2. fastCoreWeighted (LocalGini weights)
%   3. Swiftcore (LocalGini weights + reduction)

gurobiHome = getenv('GUROBI_HOME');
if ~isempty(gurobiHome) && exist(fullfile(gurobiHome, 'matlab'), 'dir')
    addpath(fullfile(gurobiHome, 'matlab'));
end

% Add fastCoreWeighted from spectra-src FIRST, then put spectra_69d2281 on VERY TOP
addpath(genpath('/scratch/pavan/spectra-src/src'));
spNewPath = '/scratch/pavan/q_4_4/spectra_69d2281/src';
addpath(spNewPath, '-begin');
addpath(genpath('/scratch/pavan/cobratoolbox'));
initCobraToolbox(false, 'agent');
changeCobraSolver('gurobi', 'all');
changeCobraSolverParams('LP', 'feasTol', 1e-8);
changeCobraSolverParams('LP', 'optTol', 1e-8);

assert(startsWith(which('spectraME'), spNewPath), 'Wrong spectraME on path!');
assert(startsWith(which('minNet'), spNewPath), 'Wrong minNet on path!');
assert(~isempty(which('fastCoreWeighted')), 'fastCoreWeighted missing!');
assert(~isempty(which('swiftcore')), 'swiftcore missing!');

outDir = '/scratch/pavan/q_4_4/rerun_69d2281/timing_chunks';
if ~exist(outDir, 'dir'); mkdir(outDir); end
outFile = fullfile(outDir, sprintf('tchunk_%03d.mat', taskId));
if exist(outFile, 'file')
    fprintf('Timing chunk %s already exists. Exiting.\n', outFile);
    return;
end

S = load('/scratch/pavan/q_4_4/Results_BuildDepMapModels.mat', 'model', 'RxnImp', 'context');
model = S.model;
RxnImp = S.RxnImp;
nTotal = size(RxnImp, 2);
nRxns  = numel(model.rxns);

% Warm up JIT + Gurobi MEX + QR BLAS once on idx=1 before timing slice
c1 = RxnImp(:, 1);
w1 = zeros(nRxns, 1); w1(c1 < 1) = 1 - c1(c1 < 1);
core1 = find(c1 >= 1);
rng(1, 'twister'); [~, ~] = spectraME(model, core1, 1e-4, [], w1);
[~, ~] = evalc('fastCoreWeighted(core1, model, w1, 1e-4);');
[~] = swiftcore(model, core1, w1, 1e-10, 1);

chunkSize = ceil(nTotal / nTasks);
idxStart  = (taskId - 1) * chunkSize + 1;
idxEnd    = min(taskId * chunkSize, nTotal);
if idxStart > nTotal; return; end
indices = idxStart:idxEnd;
nSlice  = numel(indices);

spTimeSlice  = nan(nSlice, 1);
fcwTimeSlice = nan(nSlice, 1);
swTimeSlice  = nan(nSlice, 1);
spLpsSlice   = nan(nSlice, 1);

for k = 1:nSlice
    idx = indices(k);
    curr_core = RxnImp(:, idx);
    weights = zeros(nRxns, 1);
    weights(curr_core < 1) = 1 - curr_core(curr_core < 1);
    core = find(curr_core >= 1);

    rng(idx, 'twister');
    t0 = tic;
    [~, lps] = spectraME(model, core, 1e-4, [], weights);
    spTimeSlice(k) = toc(t0);
    spLpsSlice(k)  = double(lps);

    t0 = tic;
    [~, ~] = evalc('fastCoreWeighted(core, model, weights, 1e-4);');
    fcwTimeSlice(k) = toc(t0);

    t0 = tic;
    [~] = swiftcore(model, core, weights, 1e-10, 1);
    swTimeSlice(k) = toc(t0);

    fprintf('[%d/%d] idx=%d: sp=%.2fs (%d LPs) | fcw=%.2fs | sw=%.2fs\n', ...
        k, nSlice, idx, spTimeSlice(k), lps, fcwTimeSlice(k), swTimeSlice(k));
end

save(outFile, 'indices', 'spTimeSlice', 'fcwTimeSlice', 'swTimeSlice', 'spLpsSlice');
fprintf('Finished timing slice %d/%d.\n', taskId, nTasks);
end
