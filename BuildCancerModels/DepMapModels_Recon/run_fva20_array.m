function run_fva20_array(taskId)
% run_fva20_array  Run exact unconstrained FVA (optPercentage = 0) on one of the
% 20 named cancer/normal SPECTRA_ME models rebuilt with commit 69d2281.
gurobiHome = getenv('GUROBI_HOME');
if ~isempty(gurobiHome) && exist(fullfile(gurobiHome, 'matlab'), 'dir')
    addpath(fullfile(gurobiHome, 'matlab'));
end
addpath('/scratch/pavan/q_4_4/spectra_69d2281/src');
addpath(genpath('/scratch/pavan/cobratoolbox'));
initCobraToolbox(false, 'agent');
changeCobraSolver('gurobi', 'all');
changeCobraSolverParams('LP', 'feasTol', 1e-8);
changeCobraSolverParams('LP', 'optTol', 1e-8);

all20 = {
    'ACH-000376', 'ACH-000649', 'ACH-000739', 'ACH-000681', 'ACH-001295', ...
    'ACH-000108', 'ACH-000177', 'ACH-000019', 'ACH-000001', 'ACH-001961', ...
    'ACH-003150', 'ACH-003143', 'ACH-003155', 'ACH-000064', 'ACH-002466', ...
    'ACH-003161', 'ACH-001453', 'ACH-000642', 'ACH-000494', 'ACH-001854'
};

ach = all20{taskId};
outRoot = '/scratch/pavan/q_4_4/rerun_69d2281';
fvaDir  = fullfile(outRoot, 'fva');
if ~exist(fvaDir, 'dir'); mkdir(fvaDir); end

outMat = fullfile(fvaDir, [ach '.mat']);
if exist(outMat, 'file')
    fprintf('FVA for %s already exists. Exiting.\n', ach);
    return;
end

M = load(fullfile(outRoot, 'models', 'spectra', [ach '.mat']), 'm');
m = M.m;
nRxns = numel(m.rxns);
fprintf('=== FVA (optPercentage = 0) for %s (%d rxns, %d mets) ===\n', ach, nRxns, numel(m.mets));

t0 = tic;
[minFlux, maxFlux] = fluxVariability(m, 0, 'max', m.rxns);
dt = toc(t0);

rxns = m.rxns;
save(outMat, 'rxns', 'minFlux', 'maxFlux', 'dt');

% Also write CSV and XLSX matching archived format ['', 'minimum', 'maximum']
T = table(rxns(:), minFlux(:), maxFlux(:), 'VariableNames', {'rxn', 'minimum', 'maximum'});
writetable(T, fullfile(fvaDir, [ach '.csv']));
cellOut = [{'', 'minimum', 'maximum'}; [rxns(:), num2cell(minFlux(:)), num2cell(maxFlux(:))]];
writecell(cellOut, fullfile(fvaDir, [ach '.xlsx']), 'Sheet', 'Sheet1');
fprintf('Finished FVA for %s in %.1f s.\n', ach, dt);
end
