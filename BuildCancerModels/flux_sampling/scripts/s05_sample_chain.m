%% s05_sample_chain.m
% Reviewer-2 response, Phase 2: MCMC flux sampling of one model, one chain.
%
% Reviewer 2 asked for "flux distributions of each reaction using an MCMC
% method", compared statistically between cancer and normal. This script draws
% one independent chain from the uniform distribution over the steady-state
% flux polytope
%
%       { v : S v = 0,  lb <= v <= ub }
%
% of one context-specific model, using constrained Riemannian Hamiltonian Monte
% Carlo (PolytopeSampler, shipped with the COBRA Toolbox as the 'RHMC' sampler).
% RHMC is used rather than CHRR because it terminates on an effective-sample-size
% criterion instead of a hand-tuned thinning interval, so the chain certifies its
% own mixing.
%
% NO extra constraints are applied: the polytope is exactly the one FVA saw when
% the published flux span ratios were computed, so the sampling result and the
% FSR result are directly comparable.
%
% Chains are run as separate processes with different seeds so that a
% between-chain Gelman-Rubin R-hat can be computed in the merge step; that is a
% stronger convergence statement than any single chain can make about itself.
%
% Driven by environment variables (set by the SLURM array wrapper):
%   MODEL_ID   e.g. ACH-000681
%   CHAIN      integer chain index, seeds the RNG
%   NSAMPLES   effective samples requested from this chain
%   MAXTIME    wall-clock cap in seconds
%   OUTDIR     where to write <MODEL_ID>_chain<CHAIN>.mat

root = '/home/pavan/scratch/SPECTRA_review_2';
psDir = '/home/pavan/scratch/SPECTRA_review_2/vendor/PolytopeSamplerMatlab';
addpath(genpath(fullfile(psDir,'code')));
addpath(fullfile(psDir,'bin'));

modelID  = getenv('MODEL_ID');
chain    = str2double(getenv('CHAIN'));
nSamples = str2double(getenv('NSAMPLES'));
maxTime  = str2double(getenv('MAXTIME'));
outDir   = getenv('OUTDIR');

if isempty(modelID),  error('MODEL_ID not set'); end
if isnan(chain),      chain = 1; end
if isnan(nSamples),   nSamples = 2500; end
if isnan(maxTime),    maxTime = 6*3600; end
if isempty(outDir),   outDir = fullfile(root,'sampling'); end
if ~exist(outDir,'dir'); mkdir(outDir); end

outFile = fullfile(outDir, sprintf('%s_chain%d.mat', modelID, chain));
if exist(outFile,'file')
    fprintf('%s exists, nothing to do.\n', outFile);
    return
end

fprintf('=== %s chain %d | N=%d | maxTime=%ds ===\n', modelID, chain, nSamples, maxTime);

S = load(fullfile(root,'models','sprintcore',[modelID '.mat']));
m = S.m;
fprintf('model: %d reactions, %d metabolites\n', numel(m.rxns), numel(m.mets));

% every constraint in these models is an equality with a zero right-hand side,
% checked in s04; assert rather than assume
assert(all(m.csense == 'E'), 'Model has non-equality constraints; P below would be wrong.');
assert(all(m.b == 0), 'Model has a non-zero b; P below would be wrong.');

P = struct();
P.Aeq = m.S;
P.beq = m.b(:);
P.lb  = m.lb(:);
P.ub  = m.ub(:);

opts = default_options();
opts.seed     = 1000*chain + 7;   % distinct, reproducible stream per chain
opts.maxTime  = maxTime;
opts.nWorkers = 1;                % parallelism comes from the SLURM array
opts.module   = setdiff(opts.module, {'ProgressBar'}, 'stable');  % keep logs readable

t0 = tic;
o = sample(P, nSamples, opts);
elapsed = toc(t0);

samples = o.samples;
if iscell(samples); samples = samples{1}; end
fprintf('drew %d points in %.1f s (prepare %.1f s)\n', ...
    size(samples,2), elapsed, o.prepareTime);

% within-chain ESS, so a bad chain can be spotted before the merge
ess = effective_sample_size(samples);
fprintf('within-chain ESS: min %.1f | median %.1f | max %.1f\n', ...
    min(ess), median(ess), max(ess));

rxns = m.rxns;
nRxns = numel(rxns);
assert(size(samples,1) == nRxns, 'sample dimension does not match reaction count');

% single precision halves the footprint; fluxes here are O(1e-6..1e3) and the
% downstream tests are rank-based, so ~7 significant digits is ample
samples = single(samples);
sampleTime = elapsed;
prepareTime = o.prepareTime;
seed = opts.seed;

save(outFile,'samples','rxns','ess','sampleTime','prepareTime','seed', ...
     'modelID','chain','nSamples','-v7');   % v7, not v7.3: readable by scipy.io.loadmat downstream
fprintf('wrote %s\n', outFile);
