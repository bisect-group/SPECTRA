function T = bench_paper_runtime(cfg)
% bench_paper_runtime  SPECTRA paper's consistency-check runtime benchmark,
% extended to FASTERCC.
%
%   T = bench_paper_runtime()
%
% This reproduces the design of RuntimeComparisons/ForSPECTRA_CCResults.m on
% the paper_codes branch - random subnetworks of the flux-consistent Recon3D
% core at increasing size - and adds FASTERCC as a fourth method.
%
% It is deliberately a SEPARATE benchmark from the FASTERCC-paper one
% (bench_run_one / expanded multi-cluster models). Two different questions:
%   this one  - how does runtime scale with model size, 20 .. 10,600 rxns
%   that one  - how does runtime scale with blocked-reaction count, at the
%               very large sizes FASTERCC was designed for (up to ~184k)
%
% Differences from the original script, all deliberate:
%   1. FASTERCC added (revopt = 0; see below).
%   2. Seeded. The original called randsample unseeded, so its subnetworks
%      were not reproducible. Seeds derive from cfg.baseSeed, so the exact
%      same 530 subnetworks can be regenerated.
%   3. Checkpointed after every size, so a wall-clock timeout costs one size
%      rather than the whole run, and the job can resume.
%   4. Gurobi throughout. The original already used Gurobi; FASTERCC now does
%      too, via fastercc_lp_backend/ (its shipped cplexlp layer is unusable
%      here - see bench_config).
%
% FASTERCC is run with revopt = 0. revopt = 1 cannot change the consistent set
% (fasterCC.m:322 finalises A, the revopt block starts at :330) - it only
% tightens bounds, via a one-reaction-at-a-time LP loop that cost ~40x the
% runtime in smoke testing. Timing it against methods that return only a
% consistent set would charge FASTERCC for work the others never do. Set
% cfg.paperIncludeRevopt1 = true to measure it anyway.

if nargin < 1 || isempty(cfg), cfg = bench_config('smoke'); end

outDir = fullfile(cfg.root, 'results', 'paper_runtime');
if ~exist(outDir, 'dir'), mkdir(outDir); end
ckpt = fullfile(outDir, 'checkpoint.mat');

%% Base model -------------------------------------------------------------
S = load(fullfile(cfg.baseDir, 'consRecon3D.mat'));
model = S.consRecon3D;
nBase = numel(model.rxns);
fprintf('Base model: %d reactions, %d metabolites\n', nBase, numel(model.mets));

changeCobraSolver(cfg.solver, 'all', 0);
fprintf('Solver: %s\n', cfg.solver);

tol = cfg.epsilon;
n   = 20:20:nBase;              % 530 sizes for a 10,600-reaction base
fprintf('Sizes: %d .. %d step 20  (%d models)\n\n', n(1), n(end), numel(n));

arms = {'FastCC','SwiftCC','SPECTRA_CC','FASTERCC'};
if isfield(cfg,'paperIncludeRevopt1') && cfg.paperIncludeRevopt1
    arms{end+1} = 'FASTERCC_revopt1';
end

%% Resume from checkpoint if present --------------------------------------
if exist(ckpt, 'file')
    C = load(ckpt);
    R = C.R; startIdx = C.nextIdx;
    fprintf('Resuming from checkpoint at size index %d/%d\n\n', startIdx, numel(n));
else
    R = struct('size', {}, 'arm', {}, 'runtime', {}, 'lps', {}, ...
               'nCons', {}, 'seed', {});
    startIdx = 1;
end

%% Sweep ------------------------------------------------------------------
for i = startIdx:numel(n)
    seed = cfg.baseSeed + i;
    rng(seed, 'twister');
    rxn_ids = sort(randsample(nBase, n(i)));
    m = removeRxns(model, model.rxns(setdiff(1:nBase, rxn_ids)));
    % swiftcc takes rev explicitly; set it after subsetting so it cannot be
    % stale relative to the trimmed bounds.
    m.rev = double(m.lb < 0);

    fprintf('[%3d/%3d] n=%5d (%d rxns after trim)\n', i, numel(n), n(i), numel(m.rxns));

    for a = 1:numel(arms)
        arm = arms{a};
        t = NaN; lps = NaN; nCons = NaN;
        try
            switch arm
                case 'FastCC'
                    tic; [A, ~, ~, lps] = fastcc(m, tol); t = toc;
                    nCons = numel(A);
                case 'SwiftCC'
                    tic; A = swiftcc(m.S, m.rev, cfg.solver); t = toc;
                    nCons = nnz(A);     % returns a 0-1 indicator vector
                case 'SPECTRA_CC'
                    tic; [A, lps] = spectraCC(m, tol, 'stoichiometry', seed); t = toc;
                    nCons = numel(A);
                case 'FASTERCC'
                    tic; [A, ~, nLP] = fasterCC(m, tol, 0, 0); t = toc;
                    nCons = numel(A); lps = nLP.total;
                case 'FASTERCC_revopt1'
                    tic; [A, ~, nLP] = fasterCC(m, tol, 0, 1); t = toc;
                    nCons = numel(A); lps = nLP.total;
            end
            fprintf('    %-18s %8.3f s  |cons|=%6d  LPs=%s\n', ...
                    arm, t, nCons, num2str(lps));
        catch ME
            fprintf(2, '    %-18s FAILED: %s\n', arm, ME.message);
        end

        R(end+1) = struct('size', n(i), 'arm', arm, 'runtime', t, ...
                          'lps', lps, 'nCons', nCons, 'seed', seed); %#ok<AGROW>
    end

    % Checkpoint after every size, so a timeout costs one size, not the run.
    nextIdx = i + 1; %#ok<NASGU>
    save(ckpt, 'R', 'nextIdx', '-v7.3');
end

%% Results ----------------------------------------------------------------
T = struct2table(R);
T.Properties.VariableNames = {'SizeOfModel','Algo','Runtime','LPS','nConsistent','Seed'};
T = T(:, {'Algo','SizeOfModel','Runtime','LPS','nConsistent','Seed'});

writetable(T, fullfile(outDir, 'paper_runtime.csv'));
save(fullfile(outDir, 'paper_runtime.mat'), 'T', 'n', 'arms', '-v7.3');
fprintf('\nWrote %s  (%d rows)\n', fullfile(outDir,'paper_runtime.csv'), height(T));

%% Summary ----------------------------------------------------------------
fprintf('\n=== median runtime (s) over all sizes ===\n');
for a = 1:numel(arms)
    v = T.Runtime(strcmp(T.Algo, arms{a}) & ~isnan(T.Runtime));
    fprintf('  %-18s %8.3f   (total %8.1f s)\n', arms{a}, median(v), sum(v));
end

fprintf('\n=== median LPs solved ===\n');
for a = 1:numel(arms)
    v = T.LPS(strcmp(T.Algo, arms{a}) & ~isnan(T.LPS));
    if isempty(v), fprintf('  %-18s  n/a\n', arms{a});
    else,          fprintf('  %-18s %8.1f\n', arms{a}, median(v)); end
end

% Consistent-set agreement. SwiftCC is reported separately because the
% FASTERCC paper excludes it for over-reporting consistent reactions, and
% this benchmark is where that claim can actually be checked.
fprintf('\n=== consistent-set size vs FastCC (negative = fewer) ===\n');
ref = T(strcmp(T.Algo,'FastCC'), {'SizeOfModel','nConsistent'});
for a = 1:numel(arms)
    if strcmp(arms{a},'FastCC'), continue; end
    sub = T(strcmp(T.Algo, arms{a}), {'SizeOfModel','nConsistent'});
    J = innerjoin(ref, sub, 'Keys','SizeOfModel');
    d = J.nConsistent_sub - J.nConsistent_ref;
    fprintf('  %-18s identical on %d/%d sizes, max diff %+d\n', ...
            arms{a}, sum(d == 0), numel(d), maxabs(d));
end
end

function v = maxabs(d)
d = d(~isnan(d));
if isempty(d), v = 0; else, [~,k] = max(abs(d)); v = d(k); end
end
