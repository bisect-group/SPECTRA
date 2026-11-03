% Adjudicate the experiment-1 disagreements from the feasTol = 1e-9 pass.
%
% The default-tolerance pass left SPECTRA_CC overclaiming 5 reactions and
% SwiftCC 5, with FASTERCC underclaiming 2. The open question is whether
% that is the tolerance or the methods: the community sweeps showed
% SPECTRA's answer is least stable at loose tolerances, so the overclaims
% may simply be solver slop read as flux.
%
% This scores the OTHER pass - same construction, same arms, feasTol pinned
% to 1e-9 - and the FVA is run at 1e-9 too, because an adjudicator has to
% answer under the numerics it is adjudicating.
cfg = bench_config('smoke');
bench_addpaths(cfg);
changeCobraSolver(cfg.solver, 'all', 0);
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('LP', 'optTol',  1e-9);

d = dir(fullfile(cfg.root,'results','variants_papersizes_*_a4'));
d = d(~contains({d.name}, 'dtol'));          % the pinned-1e-9 family only
fprintf('scanning %d directories (feasTol 1e-9 family)\n', numel(d));

acc = {};
for i = 1:numel(d)
    T = check_disagreements_fva(cfg, fullfile(d(i).folder, d(i).name));
    if ~isempty(T), acc{end+1} = T; end %#ok<SAGROW>
end
if isempty(acc), fprintf('\nNo disagreements found.\n'); return; end
A = vertcat(acc{:});
writetable(A, fullfile(cfg.root,'exp1_tight_disagreement_fva.csv'));
fprintf('\n===== POOLED, feasTol 1e-9 =====\n');
fprintf('%-16s %12s %12s\n','arm','overclaim','underclaim');
arms = unique(A.arm);
for a = 1:numel(arms)
    s = A(strcmp(A.arm, arms{a}), :);
    fprintf('%-16s %12d %12d\n', arms{a}, sum(s.overclaim), sum(s.underclaim));
end
