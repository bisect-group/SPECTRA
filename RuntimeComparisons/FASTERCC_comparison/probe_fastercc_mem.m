function probe_fastercc_mem(cfg)
% probe_fastercc_mem  Where inside fasterCC does the k=20 OOM happen?
%
% All twelve k=20 instances died OUT_OF_MEMORY at 64 GB after spectracc had
% already finished. The SLURM record says only that the process was killed;
% it does not say which part of fasterCC was holding the memory.
%
% fasterCC is NOT modified. Instead its stages are reproduced here in order,
% calling its own stock functions, with RSS sampled between each. The stage
% that is running when RSS crosses the cap is the answer.
%
% RSS comes from /proc/self/statm, which is the resident set the kernel will
% actually kill on, rather than MATLAB's own accounting of what it thinks it
% has allocated.
%
% Three sizes are probed so the growth can be attributed to a stage rather
% than read off a single failure: k=10 completes within 64 GB, k=20 does not,
% and k=5 anchors the low end.

if nargin < 1 || isempty(cfg), cfg = bench_config('smoke'); end
bench_addpaths(cfg);
changeCobraSolver(cfg.solver, 'all', 0);
changeCobraSolverParams('LP', 'feasTol', cfg.feasTol);
changeCobraSolverParams('LP', 'optTol',  cfg.optTol);
eps_ = cfg.epsilon;

for k = [5 10 20]
    fprintf('\n================ k = %d ================\n', k);
    [model, ~] = bench_make_instance(k, 500, 1, cfg);
    model_org = model;                                  %#ok<NASGU>
    fprintf('%-34s %8d rxns   RSS %7.1f GB\n', 'loaded', numel(model.rxns), rss());

    % --- fasterCC.m:39-55, save and relax bounds --------------------------
    orig_lb = model.lb; orig_ub = model.ub;             %#ok<NASGU>
    changed_lb = model.lb > 0; changed_ub = model.ub < 0; %#ok<NASGU>
    model.lb(model.lb > 0) = 0;
    model.ub(model.ub < 0) = 0;
    fprintf('%-34s %8s        RSS %7.1f GB\n', 'bounds relaxed', '', rss());

    % --- fasterCC.m:62-67, the structural pass ----------------------------
    t0 = tic; model = findDeadEndsFastbox(model);
    fprintf('%-34s %8d rxns   RSS %7.1f GB  %6.1f s\n', ...
            'findDeadEndsFastbox', numel(model.rxns), rss(), toc(t0));

    t0 = tic; model = structureAnalyseFastbox(model);
    fprintf('%-34s %8d rxns   RSS %7.1f GB  %6.1f s\n', ...
            'structureAnalyseFastbox', numel(model.rxns), rss(), toc(t0));

    t0 = tic; model = fixIrr_rFASTCORMICS(model);
    fprintf('%-34s %8d rxns   RSS %7.1f GB  %6.1f s\n', ...
            'fixIrr_rFASTCORMICS', numel(model.rxns), rss(), toc(t0));

    % --- fasterCC.m:96, the first batch LP7 -------------------------------
    % This is the prime suspect: LP7 on the FULL irreversible set builds a
    % problem with one extra variable and two extra constraint rows per
    % reaction in J, so its matrix is several times the size of S.
    J = find(model.lb >= 0);
    fprintf('%-34s %8d in J    RSS %7.1f GB\n', 'about to call LP7 (batch)', numel(J), rss());
    t0 = tic;
    try
        V = LP7_rFASTCORMICS(J, model, eps_); %#ok<NASGU>
        fprintf('%-34s %8s        RSS %7.1f GB  %6.1f s\n', ...
                'LP7 batch returned', '', rss(), toc(t0));
        clear V
    catch ME
        fprintf('%-34s FAILED after %.1f s: %s\n', 'LP7 batch', toc(t0), ME.message);
    end
    clear model model_org orig_lb orig_ub changed_lb changed_ub J
    fprintf('%-34s %8s        RSS %7.1f GB\n', 'cleared', '', rss());
end
end

% -------------------------------------------------------------------------
function g = rss()
% Resident set in GB, read from /proc - the number the OOM killer uses.
fid = fopen('/proc/self/statm', 'r');
if fid < 0, g = NaN; return; end
v = sscanf(fgetl(fid), '%f'); fclose(fid);
g = v(2) * 4096 / 1024^3;     % field 2 is resident pages
end
