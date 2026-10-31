function diag_fastercc_crash(cfg)
% diag_fastercc_crash  Is the fastercc_r0 failure ours or theirs?
%
% The exp-1 pilot lost fastercc_r0 at rm=10000 with
%     solveCobraLP: gurobi reporting OPTIMAL but no solution
% and that level is also where the model goes fully blocked (|cons| = 0).
%
% Two candidate causes, and they need different responses:
%   1. feasTol = 1e-9. We introduced that in this round. If the failure
%      disappears at COBRA's default, it is a tolerance we chose and the
%      right fix is to loosen it - a self-inflicted artefact must not be
%      reported as a FASTERCC defect.
%   2. The degenerate model. If it fails at every tolerance, it is a genuine
%      edge case in fasterCC on fully-blocked input, and reportable.
%
% Also reports where along the sweep the model actually collapses, since
% levels past that point compare methods on an empty answer and may not be
% worth running at all.

if nargin < 1 || isempty(cfg), cfg = bench_config('smoke'); end
bench_addpaths(cfg);
changeCobraSolver(cfg.solver, 'all', 0);

S    = load(fullfile(cfg.baseDir, 'consRecon3D.mat'));
base = S.consRecon3D;

%% 1. the crash, at three tolerances -------------------------------------
seed = cfg.baseSeed + 1000*10000 + 1;    % matches recon_rm10000_rep01
rng(seed, 'twister');
drop  = sort(randsample(numel(base.rxns), 10000));
model = removeRxns(base, base.rxns(drop));
fprintf('\nrm=10000 instance: %d rxns, %d mets\n\n', ...
        numel(model.rxns), numel(model.mets));

for tol = [1e-9 1e-7 1e-6]
    changeCobraSolverParams('LP', 'feasTol', tol);
    changeCobraSolverParams('LP', 'optTol',  tol);
    fprintf('feasTol=%.0e  ', tol);
    try
        t0 = tic; [A,~,nLP] = fasterCC(model, cfg.epsilon, 0, 0); t = toc(t0);
        fprintf('fastercc_r0 OK  %.2f s  |cons|=%d  LPs=%d\n', t, numel(A), nLP.total);
    catch ME
        fprintf('fastercc_r0 FAILED: %s\n', ME.message);
    end
    try
        t0 = tic; [c,l] = spectraCC(model, cfg.epsilon); t = toc(t0);
        fprintf('              spectracc   OK  %.2f s  |cons|=%d  LPs=%d\n', ...
                t, numel(c), l);
    catch ME
        fprintf('              spectracc   FAILED: %s\n', ME.message);
    end
end

%% 2. where does the model collapse? -------------------------------------
% Cheap: spectracc only, one draw per level. Levels whose consistent set is
% empty compare the methods on nothing.
changeCobraSolverParams('LP', 'feasTol', cfg.feasTol);
changeCobraSolverParams('LP', 'optTol',  cfg.optTol);
fprintf('\n%-10s %8s %10s\n', 'removed', 'nRxns', '|cons|');
for r = 0:500:10000
    rng(cfg.baseSeed + 1000*r + 1, 'twister');
    m = base;
    if r > 0
        d = sort(randsample(numel(m.rxns), r));
        m = removeRxns(m, m.rxns(d));
    end
    try
        c = spectraCC(m, cfg.epsilon);
        fprintf('%-10d %8d %10d\n', r, numel(m.rxns), numel(c));
    catch ME
        fprintf('%-10d %8d   FAILED %s\n', r, numel(m.rxns), ME.message);
    end
end
end
