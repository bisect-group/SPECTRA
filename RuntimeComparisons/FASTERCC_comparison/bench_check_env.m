function ok = bench_check_env(cfg)
% bench_check_env  Verify the cluster has everything the benchmark needs.
%
%   ok = bench_check_env()
%
% Run this FIRST on FIR. It fails loudly rather than letting an array job die
% one task at a time.

if nargin < 1 || isempty(cfg), cfg = bench_config(); end

ok = true;
fprintf('=== benchmark environment check (stage: %s) ===\n\n', cfg.stage);

%% MATLAB ------------------------------------------------------------------
v = version('-release');
fprintf('MATLAB release        : %s\n', v);
if str2double(v(1:4)) > 2020
    % The paper pins R2019b because the IBM CPLEX MATLAB API ships only up to
    % CPLEX 12.10. In practice the R2019-era cplexlp.p + cplexlink12100.mexa64
    % DO load and solve under newer MATLAB: verified on fir 2026-10-23 with
    % MATLAB 2023b.2 (cplexlp smoke test PASS, solveCobraLP stat=1). So this
    % is a note, not a warning - the cplexlp smoke test below is what decides.
    fprintf(['  note: newer than the paper''s R2019b. The CPLEX 12.10 MATLAB\n' ...
             '  API still loads here; the cplexlp smoke test below confirms it.\n']);
end

%% Solver ------------------------------------------------------------------
hasCplexlp = exist('cplexlp', 'file') > 0;
fprintf('cplexlp on path       : %s\n', tf(hasCplexlp));
if cfg.requireCplexlp && ~hasCplexlp
    fprintf(2, ['  FATAL: cplexlp not found - FASTERCC cannot run.\n' ...
                '  Add <CPLEX_Studio>/cplex/matlab/<arch> to the path.\n']);
    ok = false;
end

try
    changeCobraSolver(cfg.solver, 'LP', 0, -1);
    fprintf('COBRA LP solver       : %s (set ok)\n', cfg.solver);
catch ME
    fprintf(2, '  FATAL: could not set COBRA solver %s: %s\n', cfg.solver, ME.message);
    ok = false;
end

%% Required functions ------------------------------------------------------
groups = struct( ...
    'COBRA',    {{'initCobraToolbox','solveCobraLP','removeRxns','fastcc'}}, ...
    'FASTERCC', {{'fasterCC','findDeadEndsFastbox','structureAnalyseFastbox', ...
                  'fixIrr_rFASTCORMICS','restoreBounds','LP3_rFASTCORMICS', ...
                  'LP7_rFASTCORMICS'}}, ...
    'SPECTRA',  {{'spectraCC','forwardcc','reverse'}}, ...
    'scFASTCORMICS', {{'call_build_expanded_input_model','build_expanded_input_model'}});

fns = fieldnames(groups);
for g = 1:numel(fns)
    fprintf('\n%s:\n', fns{g});
    for i = 1:numel(groups.(fns{g}))
        name  = groups.(fns{g}){i};
        found = exist(name, 'file') > 0;
        fprintf('  %-34s %s\n', name, tf(found));
        if ~found, ok = false; end
    end
end

%% fasterCC LP-count patch -------------------------------------------------
% fastercc_lpcount.patch adds a 3rd output, nLP. Without it, LP counts for
% the FASTERCC arms are unavailable and bench_run_one will error.
fprintf('\nfasterCC LP-count patch:\n');
if exist('fasterCC', 'file') > 0
    try
        n = nargout('fasterCC');
        fprintf('  fasterCC outputs                  %d\n', n);
        if n >= 3
            fprintf('  patch applied                     YES\n');
        else
            fprintf(2, ['  patch applied                     NO\n' ...
                        '  FATAL: apply fastercc_lpcount.patch in %s\n' ...
                        '    cd $FASTERCC_SRC && patch -p0 < <harness>/fastercc_lpcount.patch\n'], ...
                        cfg.fasterccSrc);
            if cfg.requireFastercc3Out, ok = false; end
        end
    catch ME
        fprintf(2, '  could not inspect fasterCC: %s\n', ME.message);
    end
end

%% Smoke test: cplexlp actually solves ------------------------------------
if hasCplexlp
    fprintf('\ncplexlp smoke test    : ');
    try
        x = cplexlp([-1; -1], [1 1], 1, [], [], [0; 0], [1; 1]);
        if abs(sum(x) - 1) < 1e-6, fprintf('PASS\n');
        else, fprintf(2, 'FAIL (unexpected solution)\n'); ok = false; end
    catch ME
        fprintf(2, 'FAIL (%s)\n', ME.message); ok = false;
    end
end

%% Directories -------------------------------------------------------------
fprintf('\nDirectories:\n');
for d = {cfg.baseDir, cfg.modelDir, cfg.resultDir, cfg.logDir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
    fprintf('  %s\n', d{1});
end

fprintf('\nCores visible         : %d\n', feature('numcores'));
fprintf('\n=== %s ===\n', ternary(ok, 'ENVIRONMENT OK', 'ENVIRONMENT INCOMPLETE'));
end

% ------------------------------------------------------------------ helpers
function s = tf(b),            if b, s = 'FOUND'; else, s = 'MISSING'; end, end
function s = ternary(c, a, b), if c, s = a;       else, s = b;         end, end
