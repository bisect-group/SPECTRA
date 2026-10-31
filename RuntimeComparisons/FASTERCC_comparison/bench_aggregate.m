function T = bench_aggregate(cfg)
% bench_aggregate  Collect per-instance results into one table and summarise.
%
%   T = bench_aggregate()
%   T = bench_aggregate(bench_config('pilot'))
%
% Writes results_long_<stage>.csv (one row per instance x arm) and prints the
% summary needed to compare against the paper's claims:
%   - median runtime and LP count by number of blocked reactions
%   - median per-instance speedup relative to the reference arm
%     (paper vs FASTCC: 2.3-5.2x on consistent networks, 17.5-31.5x at 500)
%   - consistent-set agreement between arms
%
% Agreement is reported, not adjudicated. With FASTCC and SWIFTCC excluded
% there is no independent third method here, so a disagreement between
% SPECTRA_CC and FASTERCC cannot be resolved from this benchmark alone.
% Disagreeing instances are listed for separate FVA follow-up.

if nargin < 1 || isempty(cfg), cfg = bench_config(); end

files = dir(fullfile(cfg.resultDir, '*.mat'));
assert(~isempty(files), 'No results in %s', cfg.resultDir);

arms = cfg.arms(cell2mat(cfg.arms(:,2)), 1);
rows = {};

for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name));
    r = S.res;
    for a = 1:numel(arms)
        arm = arms{a};
        if ~isfield(r, arm), continue; end
        e = r.(arm);
        rows(end+1, :) = { r.name, arm, r.k, r.nRemove, r.rep, ...
                           r.nRxns, r.nMets, e.time, e.lps, e.lp7, e.lp3, ...
                           e.nCons, string(e.error) }; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'instance','arm','clusters','nRemoved','rep','nRxns','nMets', ...
     'runtimeSec','nLPs','nLP7','nLP3','nConsistent','error'});
outCsv = fullfile(cfg.root, sprintf('results_long_%s.csv', cfg.stage));
writetable(T, outCsv);
fprintf('Wrote %s  (%d rows)\n\n', outCsv, height(T));

ok = T.error == "";
levels = unique(T.nRemoved)';

%% Runtime and LP count by condition --------------------------------------
printMatrix('median runtime (s)', T, ok, arms, levels, 'runtimeSec');
printMatrix('median LPs solved',  T, ok, arms, levels, 'nLPs');

%% Speedup vs the reference arm -------------------------------------------
ref = cfg.verifyAgainst;
fprintf('\n=== median per-instance speedup vs %s (>1 = faster than %s) ===\n', ref, ref);
printHeader(levels);
for a = 1:numel(arms)
    if strcmp(arms{a}, ref), continue; end
    line = sprintf('%-22s', arms{a});
    for L = levels
        sub = T(T.nRemoved == L & ok, :);
        A = sub(strcmp(sub.arm, arms{a}), {'instance','runtimeSec'});
        B = sub(strcmp(sub.arm, ref),     {'instance','runtimeSec'});
        % Pair within instance: a ratio of medians is not a median of ratios.
        J = innerjoin(A, B, 'Keys', 'instance');
        line = [line sprintf('%12s', fmt(median(J.runtimeSec_B ./ J.runtimeSec_A)))]; %#ok<AGROW>
    end
    disp(line);
end

%% Set agreement -----------------------------------------------------------
fprintf('\n=== consistent-set agreement vs %s ===\n', ref);
disagree = {};
for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name)); r = S.res;
    if ~isfield(r, ref) || isempty(r.(ref).consRxns), continue; end
    refSet = r.(ref).consRxns;
    for a = 1:numel(arms)
        arm = arms{a};
        if strcmp(arm, ref) || ~isfield(r, arm) || isempty(r.(arm).consRxns), continue; end
        onlyRef = numel(setdiff(refSet, r.(arm).consRxns));
        onlyArm = numel(setdiff(r.(arm).consRxns, refSet));
        if onlyRef || onlyArm
            disagree(end+1, :) = {r.name, arm, onlyRef, onlyArm}; %#ok<AGROW>
        end
    end
end

if isempty(disagree)
    fprintf('All arms returned identical consistent sets on every instance.\n');
else
    D = cell2table(disagree, 'VariableNames', {'instance','arm','onlyInRef','onlyInArm'});
    writetable(D, fullfile(cfg.root, sprintf('disagreements_%s.csv', cfg.stage)));
    fprintf(2, '%d instance/arm pairs disagree with %s.\n', height(D), ref);
    fprintf(['Neither arm is an oracle and FASTCC/SWIFTCC are not in this\n' ...
             'run, so these need FVA adjudication separately.\n']);
    disp(head(D, 10));
end
end

% ------------------------------------------------------------------ helpers
function printMatrix(title, T, ok, arms, levels, col)
fprintf('=== %s, by blocked-reaction count ===\n', title);
printHeader(levels);
for a = 1:numel(arms)
    line = sprintf('%-22s', arms{a});
    for L = levels
        v = T.(col)(strcmp(T.arm, arms{a}) & T.nRemoved == L & ok);
        line = [line sprintf('%12s', fmt(median(v(~isnan(v)))))]; %#ok<AGROW>
    end
    disp(line);
end
fprintf('\n');
end

function printHeader(levels)
hdr = sprintf('%-22s', 'arm');
for L = levels, hdr = [hdr sprintf('%12s', sprintf('rm=%d', L))]; end %#ok<AGROW>
disp(hdr);
end

function s = fmt(x)
if isempty(x) || all(isnan(x)), s = '-'; else, s = sprintf('%.2f', x); end
end
