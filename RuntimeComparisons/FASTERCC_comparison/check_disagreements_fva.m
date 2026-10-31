function T = check_disagreements_fva(cfg, resultDir)
% check_disagreements_fva  Adjudicate stored arm disagreements with FVA.
%
%   T = check_disagreements_fva(cfg, 'results/variants_recon_..._a3')
%
% The benchmark runs deliberately carry no FVA: it is ~90% of wall time and
% would make replicates unaffordable. Instead each instance stores the
% reactions its arms disagreed about (res.disagreeRxns / res.disagreeFlags).
% This is the follow-up pass that decides who was right, and it is cheap
% because it scores only those reactions - 2 LPs each, not 2 per reaction in
% the model.
%
% A disagreement means at least one arm is WRONG, and the two ways of being
% wrong are not equally serious:
%
%   overclaim   arm says consistent, FVA says blocked. UNSOUND. A downstream
%               reconstruction keeps a reaction that can carry no flux.
%   underclaim  arm says blocked, FVA says consistent. Sound but incomplete:
%               a real reaction is silently deleted from the model.
%
% The instance is rebuilt from res.removedRxns, the stored deletion list,
% rather than replayed from the seed, so the result does not depend on this
% MATLAB's RNG matching the one that produced the run.

if nargin < 1 || isempty(cfg), cfg = bench_config('smoke'); end
if nargin < 2 || isempty(resultDir)
    error('check_disagreements_fva:dir', 'Pass the results directory to scan.');
end
if ~isfolder(resultDir)
    resultDir = fullfile(cfg.root, resultDir);
end

changeCobraSolver(cfg.solver, 'all', 0);
tol = cfg.epsilon;

S    = load(fullfile(cfg.baseDir, 'consRecon3D.mat'));
base = S.consRecon3D;

files = dir(fullfile(resultDir, '*.mat'));
files = files(~strcmp({files.name}, 'checkpoint.mat'));

rows = {};
for i = 1:numel(files)
    L = load(fullfile(resultDir, files(i).name));
    if ~isfield(L, 'res'), continue; end
    r = L.res;
    if ~isfield(r, 'nDisagree') || r.nDisagree == 0, continue; end
    if ~isfield(r, 'disagreeRxns')
        fprintf('%s: %d disagreements but no stored reaction list - skipped.\n', ...
                r.name, r.nDisagree);
        continue;
    end

    fprintf('\n=== %s : %d disagreeing of %d reactions ===\n', ...
            r.name, r.nDisagree, r.nRxns);

    % Rebuild the exact instance.
    model = base;
    if ~isempty(r.removedRxns)
        model = removeRxns(model, r.removedRxns);
    end
    % If this fails the stored result does not describe this model and no
    % verdict from it would mean anything - better to stop than to report a
    % confident answer about the wrong network.
    assert(isequal(model.rxns, r.rxns), 'check_disagreements_fva:rebuild', ...
        ['Rebuilt %s has a different reaction set than the stored run ' ...
         '(%d vs %d reactions).'], r.name, numel(model.rxns), numel(r.rxns));

    t0 = tic;
    [minF, maxF] = fluxVariability(model, 0, 'max', r.disagreeRxns);
    tF = toc(t0);

    % FVA is the ground truth: a reaction is consistent iff it can carry a
    % flux of at least tol in some direction.
    truth = (abs(minF) >= tol) | (abs(maxF) >= tol);
    fprintf('FVA on %d reactions in %.1f s\n\n', numel(r.disagreeRxns), tF);

    % Per-reaction detail.
    hdr = sprintf('%-18s %11s %11s %11s', 'reaction', 'FVA min', 'FVA max', 'FVA says');
    for a = 1:numel(r.arms), hdr = [hdr sprintf('%16s', r.arms{a})]; end %#ok<AGROW>
    disp(hdr);
    for j = 1:numel(r.disagreeRxns)
        if truth(j), verdict = 'consistent'; else, verdict = 'BLOCKED'; end
        line = sprintf('%-18s %11.4g %11.4g %11s', ...
                       r.disagreeRxns{j}, minF(j), maxF(j), verdict);
        for a = 1:numel(r.arms)
            if r.disagreeFlags(j,a), s = 'consistent'; else, s = 'BLOCKED'; end
            if r.disagreeFlags(j,a) ~= truth(j), s = [s ' <-WRONG']; end %#ok<AGROW>
            line = [line sprintf('%16s', s)]; %#ok<AGROW>
        end
        disp(line);
    end

    % Per-arm tally.
    fprintf('\n%-18s %12s %12s\n', 'arm', 'overclaim', 'underclaim');
    for a = 1:numel(r.arms)
        claimed = r.disagreeFlags(:,a);
        over    = sum(claimed  & ~truth);
        under   = sum(~claimed &  truth);
        fprintf('%-18s %12d %12d\n', r.arms{a}, over, under);
        rows(end+1,:) = {r.name, r.arms{a}, r.nRxns, r.nDisagree, over, under}; %#ok<AGROW>
    end
end

if isempty(rows)
    fprintf('\nNo stored disagreements found in %s\n', resultDir);
    T = table();
    return
end

T = cell2table(rows, 'VariableNames', ...
    {'instance','arm','nRxns','nDisagree','overclaim','underclaim'});
out = fullfile(resultDir, 'disagreement_fva.csv');
writetable(T, out);
fprintf('\nWrote %s\n', out);
end
