% Adjudicate the experiment-1 arm disagreements with FVA.
%
% 12 reactions across 4 instances, so the FVA itself is 24 LPs; the cost is
% rebuilding each instance.
%
% The solver tolerance is deliberately NOT set. These instances ran at the
% solver default, and FVA has to answer under the same numerics it is
% adjudicating - scoring a 1e-6 run with a 1e-9 FVA would measure the
% tolerance difference rather than the methods.
cfg = bench_config('smoke');
bench_addpaths(cfg);
d = dir(fullfile(cfg.root,'results','*dtol*_a4'));
acc = {};
for i = 1:numel(d)
    T = check_disagreements_fva(cfg, fullfile(d(i).folder, d(i).name));
    if ~isempty(T), acc{end+1} = T; end %#ok<SAGROW>
end
if isempty(acc)
    fprintf('\nNo disagreements found.\n');
    return
end
A = vertcat(acc{:});
writetable(A, fullfile(cfg.root,'exp1_dtol_disagreement_fva.csv'));
fprintf('\n===== POOLED OVER ALL EXPERIMENT-1 INSTANCES =====\n');
fprintf('%-16s %12s %12s\n','arm','overclaim','underclaim');
arms = unique(A.arm);
for a = 1:numel(arms)
    s = A(strcmp(A.arm, arms{a}), :);
    fprintf('%-16s %12d %12d\n', arms{a}, sum(s.overclaim), sum(s.underclaim));
end
