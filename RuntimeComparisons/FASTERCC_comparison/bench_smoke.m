function T = bench_smoke(reconPath)
% bench_smoke  End-to-end sample run on the smallest instances.
%
%   T = bench_smoke()
%   T = bench_smoke('/path/to/Recon3D.mat')
%
% RUN THIS FIRST on FIR, before the pilot and before any array job. It
% exercises every moving part on tiny models, so failures surface in minutes
% rather than after a large array has been queued:
%
%   1. environment   cplexlp loads, COBRA solver sets, fasterCC patch applied
%   2. base model    Recon3D -> flux-consistent core (expect ~10,600 rxns)
%   3. expansion     k = 1, 2 via call_build_expanded_input_model
%   4. size law      checks n(k) = a*k + b, and that n(2) == 2*nBase
%   5. all four arms spectracc, spectracc_relaxed, fastercc_r0, fastercc_r1
%   6. LP counts     from the same single run of each arm
%   7. agreement     do the arms return the same consistent set
%
% Nothing here is a timing result worth reporting: k = 1 and 2 are far below
% the sizes where the methods separate. This checks that the pipeline runs.

cfg = bench_config('smoke');

fprintf('\n########## 1. environment ##########\n');
bench_addpaths(cfg);
if ~bench_check_env(cfg)
    error('bench_smoke:env', 'Environment incomplete - fix the above first.');
end

fprintf('\n########## 2. base model ##########\n');
if nargin < 1, reconPath = ''; end
base  = bench_prepare_base(reconPath, cfg);
nBase = numel(base.rxns);

fprintf('\n########## 3-4. expansion and size law ##########\n');
manifest = bench_generate_models(cfg);

% n(2) should be exactly 2*nBase: the k-independent Trans_/Ex_ reactions
% (numel(exmets) each) replace the per-cluster exchange reactions. This is
% what makes the paper's stated 21,200 minimum come out exactly.
n2 = unique(manifest.nRxns(manifest.clusters == 2 & manifest.nRemoved == 0));
if ~isempty(n2)
    fprintf('n(2) = %d, 2*nBase = %d  -> %s\n', n2, 2*nBase, ...
            ternary(n2 == 2*nBase, 'MATCH', 'MISMATCH, revisit the size model'));
end

fprintf('\n########## 5-6. arms, runtime and LP counts ##########\n');
for i = 1:height(manifest)
    bench_run_one(manifest.clusters(i), manifest.nRemoved(i), manifest.rep(i), cfg);
end

fprintf('\n########## 7. aggregate ##########\n');
T = bench_aggregate(cfg);

fprintf('\nSmoke test complete. Artifacts under %s\n', cfg.root);
fprintf(['If all four arms ran, reported LP counts, and agree on the\n' ...
         'consistent set, move on to bench_config(''pilot'').\n']);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
