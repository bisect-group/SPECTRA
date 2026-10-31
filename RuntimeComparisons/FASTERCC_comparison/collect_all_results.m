function T = collect_all_results(cfg)
% collect_all_results  One long-format table over every benchmark run so far.
%
%   T = collect_all_results(cfg)
%
% Walks results/variants_* and flattens every stored instance into one row per
% (run, instance, arm). Runs differ in arm set, replicate count, whether FVA
% was scored and which fields existed when they were written, so every field
% is read defensively - a run that predates a field contributes NaN for it
% rather than aborting the collection.
%
% Writes all_results.csv next to the results tree, and prints the tables the
% write-up needs.

if nargin < 1 || isempty(cfg), cfg = bench_config('smoke'); end
resRoot = fullfile(cfg.root, 'results');

% Arm order for runs written before res.arms was stored. This is the order
% bench_variants used at the time; it is recoverable from the directory name
% because the order is in the tag.
legacySpectraFirst  = {'spectracc','spectracc_prep','spectracc_dir','fastercc_r0','fastercc_r1'};
legacyFasterccFirst = {'fastercc_r0','fastercc_r1','spectracc','spectracc_prep','spectracc_dir'};

d = dir(fullfile(resRoot, 'variants_*'));
d = d([d.isdir]);

rows = {};
for k = 1:numel(d)
    runName = d(k).name;
    files   = dir(fullfile(resRoot, runName, '*.mat'));
    files   = files(~strcmp({files.name}, 'checkpoint.mat'));
    if isempty(files), continue; end

    hasFVA = isempty(strfind(runName, 'nofva')); %#ok<STREMP>
    if ~isempty(strfind(runName, 'fastercc_first')) %#ok<STREMP>
        order = "fastercc_first"; legacy = legacyFasterccFirst;
    else
        order = "spectra_first";  legacy = legacySpectraFirst;
    end

    for i = 1:numel(files)
        L = load(fullfile(resRoot, runName, files(i).name));
        if ~isfield(L, 'res'), continue; end
        r = L.res;
        if isfield(r, 'arms'), arms = r.arms; else, arms = legacy; end

        for a = 1:numel(arms)
            if ~isfield(r, arms{a}), continue; end
            e = r.(arms{a});
            if ~isempty(e.error), continue; end
            rows(end+1,:) = { string(runName), order, hasFVA, ...
                string(r.name), string(r.kind), r.nRxns, r.nRemove, ...
                string(arms{a}), e.time, pick(e,'tmin'), pick(e,'tmax'), ...
                e.lps, e.nCons, ...
                pick(e,'overclaim'), pick(e,'underclaim'), pick(e,'dirMismatch'), ...
                pick(r,'nDisagree') }; %#ok<AGROW>
        end
    end
end

T = cell2table(rows, 'VariableNames', {'run','order','fva','instance','kind', ...
    'nRxns','nRemoved','arm','timeSec','tmin','tmax','nLPs','nCons', ...
    'overclaim','underclaim','dirMismatch','nDisagree'});

% Blocked fraction is the variable that turned out to explain the runtime
% crossover, so derive it once here rather than in every downstream table.
T.blockedPct = 100 * (1 - T.nCons ./ T.nRxns);

out = fullfile(cfg.root, 'all_results.csv');
writetable(T, out);
fprintf('Wrote %s (%d rows)\n\n', out, height(T));

print_all(T);
end

% -------------------------------------------------------------------------
function v = pick(s, f)
% Field value, or NaN if this run predates the field / never scored it.
% full(): several of these counters are computed from sparse logical
% indexing and stay sparse, which sprintf refuses outright.
if isfield(s, f) && ~isempty(s.(f)), v = full(double(s.(f))); else, v = NaN; end
end

% -------------------------------------------------------------------------
function print_all(T)
armOrder = {'spectracc','spectracc_prep','spectracc_dir','fastercc_r0','fastercc_r1'};

runs = unique(T.run, 'stable');
for k = 1:numel(runs)
    S = T(T.run == runs(k), :);
    fprintf('\n================================================================\n');
    fprintf('%s\n', runs(k));
    fprintf('  order=%s   FVA=%d   instances=%d\n', S.order(1), S.fva(1), ...
            numel(unique(S.instance)));
    fprintf('================================================================\n');

    inst = unique(S.instance, 'stable');
    fprintf('%-22s %8s %8s %7s | %10s %9s %7s\n', ...
            'instance / arm','nRxns','blocked%','removed','time(s)','[min,max]','LPs');
    for i = 1:numel(inst)
        I = S(S.instance == inst(i), :);
        for a = 1:numel(armOrder)
            J = I(I.arm == armOrder{a}, :);
            if isempty(J), continue; end
            if a == 1 || isempty(I(I.arm == armOrder{max(a-1,1)}, :))
                lbl = sprintf('%s', inst(i));
            else
                lbl = '';
            end
            fprintf('%-22s %8d %8.1f %7d | %10.2f %4.0f-%-4.0f %7d\n', ...
                sprintf('%s %s', lbl, armOrder{a}), J.nRxns(1), J.blockedPct(1), ...
                J.nRemoved(1), J.timeSec(1), J.tmin(1), J.tmax(1), J.nLPs(1));
        end
    end

    if S.fva(1)
        fprintf('\n  accuracy vs FVA (summed over instances)\n');
        fprintf('  %-18s %11s %11s %12s\n','arm','overclaim','underclaim','dirMismatch');
        for a = 1:numel(armOrder)
            J = S(S.arm == armOrder{a}, :);
            if isempty(J), continue; end
            fprintf('  %-18s %11s %11s %12s\n', armOrder{a}, ...
                nstr(sum(J.overclaim)), nstr(sum(J.underclaim)), nstr(sum(J.dirMismatch)));
        end
    end
end

%% The cross-run table: LP count against blocked fraction -----------------
fprintf('\n\n================================================================\n');
fprintf('LP COUNT vs BLOCKED FRACTION  (all runs pooled, sorted)\n');
fprintf('This is the cross-run result: fastercc_r0 scales with blocked\n');
fprintf('fraction, spectracc does not.\n');
fprintf('================================================================\n');
U = T(ismember(T.arm, {'spectracc','spectracc_prep','fastercc_r0'}), :);
[~, idx] = sort(U.blockedPct);
U = U(idx, :);
inst = unique(U.instance, 'stable');
fprintf('%-20s %8s %9s | %10s %10s %10s\n', ...
        'instance','nRxns','blocked%','sp LPs','prep LPs','fcc LPs');
seen = strings(0);
for i = 1:numel(inst)
    I = U(U.instance == inst(i), :);
    key = sprintf('%s_%d', inst(i), I.nRxns(1));
    if any(seen == key), continue; end
    seen(end+1) = key; %#ok<AGROW>
    fprintf('%-20s %8d %9.1f | %10s %10s %10s\n', inst(i), I.nRxns(1), ...
        I.blockedPct(1), lp(I,'spectracc'), lp(I,'spectracc_prep'), lp(I,'fastercc_r0'));
end
end

% -------------------------------------------------------------------------
function s = lp(I, arm)
J = I(I.arm == arm, :);
if isempty(J), s = '-'; else, s = sprintf('%d', round(median(J.nLPs))); end
end

function s = nstr(x)
x = full(x);
if isnan(x), s = '-'; else, s = sprintf('%d', round(x)); end
end
