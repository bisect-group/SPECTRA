function cfg = bench_config(stage)
% bench_config  Single source of truth for the SPECTRA_CC vs FASTERCC benchmark.
%
%   cfg = bench_config()          % defaults to 'pilot'
%   cfg = bench_config('smoke')   % 4 instances, pipeline check
%   cfg = bench_config('pilot')   % 12 instances, first real numbers
%   cfg = bench_config('full')    % size set after the pilot - see below
%
% Every other script reads its settings from here, so a run is fully described
% by one file. Edit this, not the scripts.
%
% Reference being reproduced and extended:
%   Pires Pacheco, Gonzalez & Sauter (2026), "FASTERCC: Accelerating Flux
%   Consistency Testing and Context-Specific Reconstruction for Large-Scale
%   Metabolic Network Models", bioRxiv 2026.03.19.712885.

if nargin < 1 || isempty(stage), stage = 'pilot'; end
cfg.stage = lower(stage);

%% ---------------------------------------------------------------- numerics
cfg.epsilon = 1e-4;        % FASTERCC uses 1e-4; SPECTRA_CC default too
cfg.solver  = 'gurobi';    % COBRA solver for BOTH methods

% Solver tolerances, applied identically to every arm so none of them runs in
% a different numerical regime. 1e-9 is Gurobi's hard floor for FeasibilityTol.
cfg.feasTol = 1e-9;
cfg.optTol  = 1e-9;

% Independent random instances per condition - the SCIENTIFIC replicate,
% distinct from cfg.timingReps which re-times one instance. Three is the
% standing choice: enough to see draw-to-draw spread, cheap enough to afford
% at the sizes the cluster sweep reaches.
cfg.instReps  = 3;
cfg.timingReps = 3;

% AGORA2 community settings (experiment 3).
% NOTE ON THE UNIVERSE: the model store holds the 5,500 reconstructions
% reachable from the shared Drive folder, whose listing endpoint caps there
% and stops mid-Salmonella. Initials span A-S only, so the tail - including
% Staphylococcus, Streptococcus, Treponema, Veillonella - is absent. Samples
% drawn here are from those 5,500, not from all of AGORA2, and the methods
% section must say so.
cfg.agoraDir        = getenv('AGORA2_DIR');
cfg.mediaOpenProb   = 0.30;   % fraction of community uptakes left open
cfg.mediaMinInflow  = 1e-6;   % log-uniform floor on an open uptake
cfg.mediaMaxInflow  = 10;     % and its ceiling

% Why Gurobi and not CPLEX, despite the paper using CPLEX:
% fasterCC hard-wires cplexlp, the IBM CPLEX MATLAB API, which exists only up
% to CPLEX 12.10. That 2019-era MEX segfaults intermittently inside
% cplexlink12100.mexa64 under MATLAB 2023b - the oldest MATLAB on fir, and
% there is no R2019b here. Three crashes in ~6 runs, while a 75-solve stress
% test passed clean: non-deterministic, which is worse across a long array run
% than an outright failure.
%
% So FASTERCC's LP layer is routed through solveCobraLP by
% fastercc_lp_backend/, and both methods run on Gurobi.
%
% REPORTING CONSEQUENCE: these timings are "FASTERCC's algorithm on a common
% Gurobi backend", not "FASTERCC as published". The change cuts both ways - it
% also removes the per-LP API overhead that SPECTRA_CC paid and FASTERCC did
% not, so the comparison is cleaner, but it is a different claim and must be
% described as such.
cfg.requireCplexlp = false;

%% -------------------------------------------------------- test-set recipe
% Derived from the expansion code, not the paper's prose, because the prose is
% self-inconsistent ("20 expanded models composed of 1-10 clusters" cannot
% reach the stated 183,920-reaction maximum).
%
% compute_expanded_input_model.m:177 builds
%     rxns = [internal_rxns x k ; intercluster_transport x k ; trans ; Ex]
% where trans = Trans_<exmet> and Ex = Ex_<exmet> are both numel(exmets) and
% do NOT scale with k. So reaction count is affine in cluster count k:
%
%     n(k) = k*(nBase - |exRxns|) + 2*|exmets|
%
% With |exRxns| ~ |exmets|, k=2 collapses to exactly 2*nBase = 21,200 (the
% paper's stated minimum) and k=19 gives ~183.9k (its stated maximum). And
% 19 k-values x 6 removal levels x 10 replicates = 1140, the paper's total.
%
% Removal levels are the paper's five plus 0. The 0 arm is the fully
% consistent control, and it matters: the paper's speedup is 17.5-31.5x at
% 500 blocked reactions but only 2.3-5.2x at 0, so the two ends bracket the
% whole claim.
% Canonical removal list. Stages select a SUBSET of this into
% cfg.removalCounts, but bench_seed always indexes the canonical list, so an
% instance name denotes the same model in every stage.
% APPENDED, never reordered. bench_seed indexes into this list, so inserting
% a value in sorted position would renumber every level after it and silently
% redefine every cached instance. The experiment-2 levels go on the end,
% where they cannot disturb indices 1-6.
cfg.allRemovalCounts = [0 10 50 100 200 500 2000 5000 10000];
cfg.removalCounts    = cfg.allRemovalCounts;
cfg.baseSeed         = 20260101;

switch cfg.stage
    case 'smoke'
        % Smallest thing that exercises every moving part.
        cfg.clusterCounts = [1 2];
        cfg.removalCounts = [0 50];
        cfg.removalReps   = 1;          % 4 instances

    case 'pilot'
        % Both ends of the size range and both ends of the blocked range, so
        % the first real speedup numbers land where the paper's claim is
        % strongest and weakest. Also calibrates SLURM wall-time and memory.
        % k = 2, 11, 20 spans 21,200 -> 183,920 reactions.
        cfg.clusterCounts = [2 11 20];
        cfg.removalCounts = [0 500];
        cfg.removalReps   = 2;          % 12 instances

    case 'full'
        % SIZE DELIBERATELY NOT SET: chosen once the pilot gives real
        % per-instance wall-clock on FIR. Candidate was
        %   cfg.clusterCounts = [1 2 4 7 10 13 16 19];
        %   cfg.removalReps   = 3;      % 8 x 6 x 3 = 144 instances
        error('bench_config:fullNotSized', ...
              ['The full run is not sized yet. Run the pilot first, then ' ...
               'set clusterCounts/removalReps here from its timings.']);

    otherwise
        error('bench_config:stage', ...
              'stage must be smoke, pilot or full (got "%s")', stage);
end

% Have bench_generate_models fit the affine size law and report the achieved
% span against the paper's 21,200-183,920.
cfg.verifySizeLaw = true;

%% ------------------------------------------------------------------- arms
%   spectracc          stock spectraCC on the instance as generated.
%
%   spectracc_relaxed  stock spectraCC on the same instance after relaxing
%                      bounds so the zero flux vector is feasible:
%                          model.lb(model.lb > 0) = 0;
%                          model.ub(model.ub < 0) = 0;
%                      This is exactly what fasterCC does to its own input
%                      before anything else (fasterCC.m:40-46), restoring
%                      bounds at the end via restoreBounds. Running SPECTRA
%                      both ways shows whether a measured gap is the algorithm
%                      or merely that FASTERCC solves a relaxed problem.
%                      SPECTRA code is untouched - only its input differs.
%
%   fastercc_r0 /      the published function's own documented 4th argument:
%   fastercc_r1            fasterCC(model, epsilon, printLevel, revopt)
%                      The paper never says which it timed. revopt=1 runs only
%                      AFTER the consistent set A is final (fasterCC.m:322 vs
%                      :330), so it cannot change A; it only tightens bounds
%                      via a one-reaction-at-a-time LP loop. revopt=0 is the
%                      like-for-like arm against a plain consistency check;
%                      revopt=1 is reported separately as the full method.
%
% FASTCC and SWIFTCC are deliberately excluded: at these sizes they are the
% multi-hour arms. The cost is that no independent third method is in the run,
% so arm disagreements cannot be adjudicated from this benchmark alone -
% bench_aggregate flags them for separate FVA follow-up rather than picking a
% winner.
cfg.arms = { ...
    'spectracc',         true;  ...
    'spectracc_relaxed', true;  ...
    'fastercc_r0',       true;  ...
    'fastercc_r1',       true;  ...
};

%% -------------------------------------------------------------- LP counts
% Both metrics come from ONE run of each arm; nothing is executed twice.
%   spectraCC  returns its LP count as its 2nd output (stock behaviour).
%   fasterCC   returns it as a 3rd output, nLP, added by fastercc_lpcount.patch
%              - 12 added lines, all counter increments, no behaviour change.
%              Apply the patch to the rFASTCORMICS checkout before running;
%              bench_check_env verifies it is applied.
cfg.requireFastercc3Out = true;

%% --------------------------------------------------------- stored outputs
% Persist the actual consistent sets, not just timings, so any later
% comparison needs no re-run.
cfg.storeConsistentSets = true;          % reaction NAMES (index-independent)
cfg.storeIndices        = true;          % and indices into the instance model
cfg.verifyAgainst       = 'fastercc_r0'; % reference arm for set comparison

%% ------------------------------------------------------------------ limits
% Enforced by the scheduler's wall limit, not inside MATLAB: a running LP
% cannot be safely interrupted from the same process.
cfg.armTimeoutSec = 6 * 3600;

%% ------------------------------------------------------------------- paths
root = getenv('BENCH_ROOT');
if isempty(root)
    root = fileparts(mfilename('fullpath'));
end
cfg.root      = root;
cfg.baseDir   = fullfile(root, 'base');                 % consistent Recon3D core
% Expanded models depend only on k, so every stage shares one cache.
cfg.modelDir  = fullfile(root, 'models');
% Results are per stage, so a pilot never overwrites a full run.
cfg.resultDir = fullfile(root, 'results', cfg.stage);
cfg.logDir    = fullfile(root, 'logs');

% External code that must be on the MATLAB path (set by env.sh on the cluster).
cfg.spectraSrc  = getenv('SPECTRA_SRC');        % spectra/src
cfg.fasterccSrc = getenv('FASTERCC_SRC');       % rFASTCORMICS/FASTERCC
cfg.scfastSrc   = getenv('SCFASTCORMICS_SRC');  % scFASTCORMICS
cfg.cobraSrc    = getenv('COBRA_SRC');          % cobratoolbox

end
