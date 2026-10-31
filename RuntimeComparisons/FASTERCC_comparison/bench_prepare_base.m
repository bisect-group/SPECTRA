function consRecon3D = bench_prepare_base(reconPath, cfg)
% bench_prepare_base  Build the flux-consistent Recon3D core used as the
% per-cluster unit of the benchmark.
%
%   consRecon3D = bench_prepare_base(reconPath)
%
% FASTERCC's benchmark replicates "10,600 reactions from Recon3D" per cluster,
% i.e. the flux-consistent core of Recon3D. This computes that core once and
% caches it, so the expansion step does not recompute it per cluster count.
%
% INPUT
%   reconPath  path to a .mat holding a Recon3D COBRA model. Defaults to
%              $RECON3D_MAT, else base/Recon3D.mat.
%
% Note: call_build_expanded_input_model runs fastcc_fastcore on whatever it is
% given anyway, so passing an already-consistent model is idempotent. Doing it
% here makes the base size explicit and checkable against the paper's 10,600.

if nargin < 2 || isempty(cfg), cfg = bench_config(); end
if ~exist(cfg.baseDir, 'dir'), mkdir(cfg.baseDir); end

outPath = fullfile(cfg.baseDir, 'consRecon3D.mat');
if exist(outPath, 'file')
    S = load(outPath); consRecon3D = S.consRecon3D;
    fprintf('Loaded cached base model: %d reactions\n', numel(consRecon3D.rxns));
    return
end

if nargin < 1 || isempty(reconPath)
    reconPath = getenv('RECON3D_MAT');
    if isempty(reconPath), reconPath = fullfile(cfg.baseDir, 'Recon3D.mat'); end
end
assert(exist(reconPath, 'file') > 0, ...
    ['Recon3D not found at %s. Set RECON3D_MAT to a .mat holding Recon3D ' ...
     '(e.g. the VMH distribution).'], reconPath);

%% Load whatever COBRA model is inside -------------------------------------
tmp = load(reconPath);
fn  = fieldnames(tmp);
model = [];
for i = 1:numel(fn)
    c = tmp.(fn{i});
    if isstruct(c) && all(isfield(c, {'S','rxns','mets','lb','ub'}))
        model = c; break
    end
end
assert(~isempty(model), 'No COBRA model structure found in %s', reconPath);
fprintf('Loaded %s: %d reactions, %d metabolites\n', ...
        reconPath, numel(model.rxns), numel(model.mets));

%% Flux-consistent core ----------------------------------------------------
changeCobraSolver(cfg.solver, 'LP', 0, -1);
fprintf('Running fastcc to extract the consistent core ... ');
t0 = tic;
A = fastcc(model, cfg.epsilon);
fprintf('%.1f s\n', toc(t0));

consRecon3D = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), A)));
n = numel(consRecon3D.rxns);
fprintf('Consistent core: %d reactions (paper states 10,600)\n', n);
if abs(n - 10600) > 200
    warning(['Consistent core is %d reactions, not ~10,600. The expanded ' ...
             'model sizes will not match the paper''s 21,200-183,920 span. ' ...
             'Check which Recon3D release is being used.'], n);
end

save(outPath, 'consRecon3D', '-v7.3');
fprintf('Saved %s\n', outPath);
end
