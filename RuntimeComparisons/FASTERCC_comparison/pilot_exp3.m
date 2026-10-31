% Pilot for experiment 3. Exercises the whole community path end to end for
% the first time: read the frozen pool, load 50 AGORA2 models, join them with
% createCommModel, find the [u] exchanges, impose a random medium, then run
% both arms.
%
% One draw, one timing - this is a correctness check, not a measurement. What
% it needs to establish before 200-member communities are worth queueing:
%   - createCommModel runs here at all (it was only just deployed, and its
%     getCobraComParams dependency comes from COBRA's SteadyCom tree)
%   - the community exchanges are found - the code looks for '(u)' and falls
%     back to '[u]', and if AGORA2 + this joiner produce neither, the medium
%     would silently apply to nothing and every reaction would stay consistent
%   - how long the join itself takes, which is pure overhead charged to no arm
%   - peak memory, since 200 members is ~4x this
cfg = bench_config('smoke');
cfg.epsilon = 1e-8;          % community setting: below the 1e-6 media floor
bench_addpaths(cfg);
bench_variants(cfg, 'community', 'spectra_first', ...
    struct('fva', false, 'reps', 1, 'instReps', 1, 'ks', 50));
