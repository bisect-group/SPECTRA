% Emit every INDIVIDUAL timing for the experiment-1 paper sizes.
%
% The per-run CSVs carry runtimeSec = median(times) for each draw, so a plot
% built from them can only ever be a mean of medians. bench_variants also
% stores the raw vector in res.<arm>.times, which is what this pulls out, so
% the figure can average all 3 draws x 3 timings = 9 measurements per cell.
cfg = bench_config('smoke');
d = dir(fullfile(cfg.root,'results','variants_papersizes_*_a4'));
fid = fopen(fullfile(cfg.root,'exp1_papersizes_timings.csv'),'w');
fprintf(fid,'tag,instance,arm,nRxns,rep,timing,runtimeSec,nLPs,nConsistent\n');
nrow = 0;
for k = 1:numel(d)
    f = dir(fullfile(d(k).folder,d(k).name,'*.mat'));
    for j = 1:numel(f)
        if strcmp(f(j).name,'checkpoint.mat'), continue; end
        S = load(fullfile(f(j).folder,f(j).name));
        res = S.res;
        nm  = erase(f(j).name,'.mat');
        rep = sscanf(nm,'n%*5d_rep%d');
        for a = 1:numel(res.arms)
            arm = res.arms{a};
            if ~isfield(res,arm), continue; end
            r = res.(arm);
            if ~isfield(r,'times'), continue; end
            t = r.times(:)';
            for q = 1:numel(t)
                if isnan(t(q)), continue; end
                fprintf(fid,'%s,%s,%s,%d,%d,%d,%.6f,%d,%d\n', ...
                        d(k).name, nm, arm, numel(res.rxns), rep, q, t(q), ...
                        full(r.lps), r.nCons);
                nrow = nrow + 1;
            end
        end
    end
end
fclose(fid);
fprintf('wrote %d timing rows -> exp1_papersizes_timings.csv\n', nrow);
