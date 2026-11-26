%% run_rhmc_new_chain.m
% Runs one RHMC PolytopeSampler chain (2,500 ESS) on a newly rebuilt SPECTRA_ME model (commit 69d2281).
psDir = '/scratch/pavan/SPECTRA_review_2/vendor/PolytopeSamplerMatlab';
addpath(genpath(fullfile(psDir, 'code')));
addpath(fullfile(psDir, 'bin'));

modelID  = getenv('MODEL_ID');
chain    = str2double(getenv('CHAIN'));
nSamples = str2double(getenv('NSAMPLES'));
maxTime  = str2double(getenv('MAXTIME'));
outDir   = getenv('OUTDIR');
modelDir = getenv('MODELDIR');

if isempty(modelID),  error('MODEL_ID not set'); end
if isnan(chain),      chain = 1; end
if isnan(nSamples),   nSamples = 2500; end
if isnan(maxTime),    maxTime = 16200; end
if isempty(outDir),   outDir = '/scratch/pavan/q_4_4/rerun_69d2281/sampling/chains'; end
if isempty(modelDir), modelDir = '/scratch/pavan/q_4_4/rerun_69d2281/models/spectra'; end
if ~exist(outDir, 'dir'); mkdir(outDir); end

outFile = fullfile(outDir, sprintf('%s_chain%d.mat', modelID, chain));
if exist(outFile, 'file')
    fprintf('%s exists, nothing to do.\n', outFile);
    return;
end

fprintf('=== %s chain %d | N=%d | maxTime=%ds ===\n', modelID, chain, nSamples, maxTime);

S = load(fullfile(modelDir, [modelID '.mat']));
m = S.m;
fprintf('model: %d reactions, %d metabolites\n', numel(m.rxns), numel(m.mets));

assert(all(m.csense == 'E'), 'Model has non-equality constraints!');
assert(all(m.b == 0), 'Model has non-zero b!');

P = struct();
P.Aeq = m.S;
P.beq = m.b(:);
P.lb  = m.lb(:);
P.ub  = m.ub(:);

opts = default_options();
opts.seed     = 1000 * chain + 7;
opts.maxTime  = maxTime;
opts.nWorkers = 1;
opts.module   = setdiff(opts.module, {'ProgressBar'}, 'stable');

t0 = tic;
o = sample(P, nSamples, opts);
elapsed = toc(t0);

samples = o.samples;
if iscell(samples); samples = samples{1}; end
fprintf('drew %d points in %.1f s (prepare %.1f s)\n', size(samples, 2), elapsed, o.prepareTime);

ess = effective_sample_size(samples);
fprintf('within-chain ESS: min %.1f | median %.1f | max %.1f\n', min(ess), median(ess), max(ess));

rxns = m.rxns;
assert(size(samples, 1) == numel(rxns), 'sample dimension mismatch!');

samples = single(samples);
sampleTime = elapsed;
prepareTime = o.prepareTime;
seed = opts.seed;

save(outFile, 'samples', 'rxns', 'ess', 'sampleTime', 'prepareTime', 'seed', ...
     'modelID', 'chain', 'nSamples', '-v7');
fprintf('wrote %s\n', outFile);
