function res = bench_run_one(k, nRemove, rep, cfg)
% bench_run_one  Run every arm once on one benchmark instance.
%
%   res = bench_run_one(10, 500, 1)
%
% ONE run per arm yields BOTH metrics:
%   runtime   wall-clock around the call, nothing else in the timed region
%   LP count  spectraCC's 2nd output; fasterCC's 3rd output (nLP), added by
%             fastercc_lpcount.patch - counter increments only
%
% Writes cfg.resultDir/<name>.mat holding runtimes, LP counts and the full
% consistent set from every arm, so later comparisons need no re-run.
%
% Per-arm timeouts are enforced by the SCHEDULER's wall limit, not in MATLAB:
% a running LP cannot be safely interrupted from the same process.

if nargin < 4 || isempty(cfg), cfg = bench_config(); end
if ~exist(cfg.resultDir, 'dir'), mkdir(cfg.resultDir); end

[model, info] = bench_make_instance(k, nRemove, rep, cfg);

outPath = fullfile(cfg.resultDir, [info.name '.mat']);
if exist(outPath, 'file')
    fprintf('%s already done, skipping\n', info.name);
    S = load(outPath); res = S.res; return
end

res          = info;
res.epsilon  = cfg.epsilon;
res.solver   = cfg.solver;
res.host     = getenv('HOSTNAME');
res.matlab   = version('-release');
res.started  = datestr(now, 31);

changeCobraSolver(cfg.solver, 'LP', 0, -1);

fprintf('=== %s : %d rxns, %d mets, %d removed ===\n', ...
        info.name, info.nRxns, info.nMets, info.nRemove);

arms = cfg.arms(cell2mat(cfg.arms(:,2)), 1);
for a = 1:numel(arms)
    arm = arms{a};
    fprintf('  %-20s ', arm);

    consIdx = []; lps = NaN; lp7 = NaN; lp3 = NaN; t = NaN; errMsg = '';

    try
        switch arm
            case 'spectracc'
                % Stock spectraCC. The seed argument is part of its published
                % signature; it only makes the random objective reproducible.
                tic;
                [consIdx, lps] = spectraCC(model, cfg.epsilon, 'stoichiometry', info.seed);
                t = toc;

            case 'spectracc_relaxed'
                % Same stock function, handed a model whose bounds have been
                % relaxed so the zero flux vector is feasible - exactly what
                % fasterCC does to its own input (fasterCC.m:40-46). Relaxing
                % happens OUTSIDE the timed region so it is not charged to
                % either method.
                relaxed = model;
                relaxed.lb(relaxed.lb > 0) = 0;
                relaxed.ub(relaxed.ub < 0) = 0;
                tic;
                [consIdx, lps] = spectraCC(relaxed, cfg.epsilon, 'stoichiometry', info.seed);
                t = toc;

            case 'fastercc_r0'
                tic;
                [consIdx, ~, nLP] = fasterCC(model, cfg.epsilon, 0, 0);
                t = toc;
                [lps, lp7, lp3] = unpackLP(nLP);

            case 'fastercc_r1'
                tic;
                [consIdx, ~, nLP] = fasterCC(model, cfg.epsilon, 0, 1);
                t = toc;
                [lps, lp7, lp3] = unpackLP(nLP);

            otherwise
                error('bench_run_one:unknownArm', 'Unknown arm "%s"', arm);
        end

        fprintf('%9.2f s   |cons| = %6d', t, numel(consIdx));
        if ~isnan(lps), fprintf('   LPs = %d', lps); end
        fprintf('\n');

    catch ME
        errMsg = sprintf('%s: %s', ME.identifier, ME.message);
        fprintf('FAILED (%s)\n', ME.message);
    end

    res.(arm).time   = t;
    res.(arm).lps    = lps;
    res.(arm).lp7    = lp7;
    res.(arm).lp3    = lp3;
    res.(arm).nCons  = numel(consIdx);
    res.(arm).error  = errMsg;
    if cfg.storeIndices
        res.(arm).consIdx = consIdx(:);
    end
    if cfg.storeConsistentSets && ~isempty(consIdx)
        % Names, not indices: survives any later re-indexing of the model.
        res.(arm).consRxns = model.rxns(consIdx);
    else
        res.(arm).consRxns = {};
    end
end

res.finished = datestr(now, 31);
save(outPath, 'res', '-v7.3');
fprintf('  -> %s\n', outPath);
end

% ------------------------------------------------------------------ helpers
function [total, lp7, lp3] = unpackLP(nLP)
% fasterCC returns nLP only once fastercc_lpcount.patch is applied.
if isstruct(nLP)
    total = nLP.total; lp7 = nLP.lp7; lp3 = nLP.lp3;
else
    total = nLP; lp7 = NaN; lp3 = NaN;
end
end
