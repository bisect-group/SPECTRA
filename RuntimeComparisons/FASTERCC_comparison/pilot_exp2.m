% Pilot for experiment 2: measures the real scaling exponent at k = 5, 10, 20
% so the k = 30 and 50 wall-times can be set from data rather than from an
% extrapolation, and reports MaxRSS so the 64G request can be confirmed
% before the large instances are committed.
cfg = bench_config('smoke');
bench_addpaths(cfg);
bench_variants(cfg, 'expanded', 'spectra_first', ...
    struct('fva', false, 'reps', 1, 'instReps', 1, ...
           'ks', [5 10 20], 'levels', [500 5000]));
