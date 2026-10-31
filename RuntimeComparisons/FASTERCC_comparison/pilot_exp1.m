% Pilot for experiment 1: sizes fastcc and swiftcc on Recon3D, which have
% never been timed here, and confirms the two new arms run at all.
% Three levels spanning the range, one draw, one timing - enough to size the
% full 21-level job without spending its budget.
cfg = bench_config('smoke');
bench_addpaths(cfg);
bench_variants(cfg, 'recon', 'spectra_first', ...
    struct('fva', false, 'reps', 1, 'instReps', 1, ...
           'arms', {{'spectracc','fastercc_r0','fastcc','swiftcc'}}, ...
           'levels', [0 5000 10000]));
