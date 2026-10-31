function manifest = bench_generate_models(cfg)
% bench_generate_models  Build the expanded models and the instance manifest.
%
%   manifest = bench_generate_models()
%
% Reproduces the test set of Pires Pacheco et al. (2026): expanded Recon3D
% models of k clusters sharing a union compartment, with reactions deleted at
% random to plant known blocked reactions. Uses call_build_expanded_input_model
% from scFASTCORMICS unmodified.
%
% What this writes:
%   expanded_k<NN>.mat   one cached expanded model per cluster count (19 files)
%   manifest.mat/.csv    the 1140 instances, as a recipe
%
% What it deliberately does NOT write: the 1140 instances themselves. A single
% top-end instance is ~184k reactions, so storing them all would run to
% terabytes. Each is regenerated on demand by bench_make_instance from its
% seed, which is reproducible and far cheaper than the expansion itself.

if nargin < 1 || isempty(cfg), cfg = bench_config(); end
if ~exist(cfg.modelDir, 'dir'), mkdir(cfg.modelDir); end

%% Base model -------------------------------------------------------------
basePath = fullfile(cfg.baseDir, 'consRecon3D.mat');
assert(exist(basePath, 'file') > 0, ...
    'Base model not found at %s. Run bench_prepare_base first.', basePath);
S = load(basePath);
baseModel = S.consRecon3D;
nBase = numel(baseModel.rxns);
fprintf('Base model: %d reactions, %d metabolites, %d genes\n', ...
        nBase, numel(baseModel.mets), numel(baseModel.genes));

% Deduplicate genes before expanding, exactly as scFASTCORMICS.m:4 does
% before it calls call_build_expanded_input_model.
%
% This is required, not cosmetic. build_expanded_input_model strips transcript
% suffixes with regexprep(model.genes,'\.[0-9]+$','') , which collapses e.g.
% 8639.1 and 8639.2 onto the same id. COBRA's generateRules then does
%     genePos(j) = find(strcmp(model.genes, genes{i}{j}));
% which errors as soon as that find returns more than one index.
nGenes0   = numel(baseModel.genes);
baseModel = simplifyDuplicatedGenesFastbox(baseModel, 1);
fprintf('After simplifyDuplicatedGenesFastbox: %d genes (was %d)\n', ...
        numel(baseModel.genes), nGenes0);
assert(numel(baseModel.rxns) == nBase, ...
    'Gene simplification changed the reaction count (%d -> %d); the size law assumes it does not.', ...
    nBase, numel(baseModel.rxns));
fprintf('\n');

%% Expand once per cluster count ------------------------------------------
ks    = cfg.clusterCounts;
sizes = zeros(numel(ks), 1);

for i = 1:numel(ks)
    k = ks(i);
    cachePath = fullfile(cfg.modelDir, sprintf('expanded_k%02d.mat', k));
    if exist(cachePath, 'file')
        E = load(cachePath, 'model');
        sizes(i) = numel(E.model.rxns);
        fprintf('k=%2d  cached            %7d rxns\n', k, sizes(i));
    else
        fprintf('k=%2d  expanding ...     ', k); t0 = tic;
        model = call_build_expanded_input_model(baseModel, k, 0);
        sizes(i) = numel(model.rxns);
        fprintf('%7d rxns  (%.1f s)\n', sizes(i), toc(t0));
        save(cachePath, 'model', '-v7.3');
    end
end

%% Verify the affine size law ---------------------------------------------
% compute_expanded_input_model.m:177 gives
%   n(k) = k*(|internal| + |intercluster|) + |trans| + |Ex|
% with |trans| = |Ex| = numel(exmets), both independent of k.
if cfg.verifySizeLaw && numel(ks) >= 2
    a = (sizes(2) - sizes(1)) / (ks(2) - ks(1));
    b = sizes(1) - a * ks(1);
    maxErr = max(abs(a * ks(:) + b - sizes));
    fprintf('\nSize law : n(k) = %.0f*k + %.0f   (max deviation %.0f rxns)\n', a, b, maxErr);
    fprintf('           a = |internal|+|intercluster| = nBase - |exRxns|\n');
    fprintf('           b = |trans| + |Ex|            = 2*|exmets|\n');
    fprintf('Achieved : %d - %d reactions\n', min(sizes), max(sizes));
    fprintf('Paper    : 21200 - 183920 reactions\n');
    if maxErr > 0.5
        warning(['Reaction count is not affine in k (max deviation %.0f). ' ...
                 'The size model needs revisiting.'], maxErr);
    end
end

%% Manifest ----------------------------------------------------------------
rows = {};
idx  = 0;
for i = 1:numel(ks)
    k = ks(i);
    for r = 1:numel(cfg.removalCounts)
        nRemove = cfg.removalCounts(r);
        for rep = 1:cfg.removalReps
            idx = idx + 1;
            % Predicted post-deletion size; the deletion itself happens in
            % bench_make_instance.
            rows(end+1, :) = { idx, ...
                sprintf('k%02d_rm%03d_rep%02d', k, nRemove, rep), ...
                k, nRemove, rep, bench_seed(k, nRemove, rep, cfg), ...
                sizes(i) - nRemove }; %#ok<AGROW>
        end
    end
end

manifest = cell2table(rows, 'VariableNames', ...
    {'id','name','clusters','nRemoved','rep','seed','nRxns'});
% Per stage: the expanded-model cache is shared, but the instance list is not.
save(fullfile(cfg.modelDir, sprintf('manifest_%s.mat', cfg.stage)), 'manifest');
writetable(manifest, fullfile(cfg.modelDir, sprintf('manifest_%s.csv', cfg.stage)));

fprintf('\n%d instances in manifest (paper reports 1140).\n', height(manifest));
fprintf('  %d cluster counts x %d removal levels x %d replicates\n', ...
        numel(ks), numel(cfg.removalCounts), cfg.removalReps);
fprintf('Note: at nRemoved=0 the replicates are timing repeats on an\n');
fprintf('      identical model, not independent models.\n');
fprintf('Manifest: %s\n', ...
        fullfile(cfg.modelDir, sprintf('manifest_%s.csv', cfg.stage)));
end
