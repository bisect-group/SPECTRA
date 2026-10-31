function [model, info] = bench_make_instance(k, nRemove, rep, cfg)
% bench_make_instance  Rebuild one benchmark instance deterministically.
%
%   [model, info] = bench_make_instance(7, 100, 3)
%
% Instances are regenerated on demand rather than stored. At the top end a
% single instance is ~184k reactions, so materialising all 1140 would run to
% terabytes; the expanded model for each cluster count is cached instead (19
% files) and the random deletion is replayed from its seed.
%
% The deletion is a pure function of (k, nRemove, rep, cfg.baseSeed), so every
% arm and every pass - timing, LP counting, any later re-analysis - sees byte
% for byte the same model.

if nargin < 4 || isempty(cfg), cfg = bench_config(); end

expPath = fullfile(cfg.modelDir, sprintf('expanded_k%02d.mat', k));
assert(exist(expPath, 'file') > 0, ...
    'Expanded model for k=%d not found at %s. Run bench_generate_models first.', ...
    k, expPath);
E = load(expPath);
model = E.model;

info.k        = k;
info.nRemove  = nRemove;
info.rep      = rep;
info.seed     = bench_seed(k, nRemove, rep, cfg);
info.nRxnsPre = numel(model.rxns);

if nRemove == 0
    % The fully consistent control arm. Replicates here are timing repeats on
    % an identical model, not different models.
    info.removed     = {};
    info.removedIdx  = [];
else
    rng(info.seed, 'twister');
    idx = sort(randsample(numel(model.rxns), nRemove));
    info.removedIdx = idx;
    info.removed    = model.rxns(idx);
    model = removeRxns(model, info.removed);
end

info.nRxns = numel(model.rxns);
info.nMets = numel(model.mets);
info.name  = sprintf('k%02d_rm%03d_rep%02d', k, nRemove, rep);
end
