function T = bench_variants(cfg, regime, order, opts)
% bench_variants  Runtime and accuracy of three SPECTRA_CC variants against
% FASTERCC, on small and large models.
%
%   T = bench_variants(cfg, 'small')   % Recon3D with random reaction removal
%   T = bench_variants(cfg, 'large')   % expanded Recon3D clusters
%   T = bench_variants(cfg, 'small', 'fastercc_first')   % reversed arm order
%
% ARM ORDER IS A CONTROL, NOT COSMETIC. Whichever arm runs first on a given
% model absorbs the one-off costs: MATLAB JIT warm-up, Gurobi environment
% startup, and first-touch paging of the model. With the default order SPECTRA
% pays all of that and FASTERCC still wins, so running FASTERCC first tests
% whether its lead survives when it carries the warm-up instead. If the gap
% closes under 'fastercc_first', the gap was partly warm-up, not algorithm.
%
% ARMS
%   spectracc        stock spectraCC.
%   spectracc_prep   FASTERCC's structural preprocessing (prep_fastercc, which
%                    calls FASTERCC's own stock functions) then stock
%                    spectraCC. Isolates how much of FASTERCC's advantage is
%                    the structural pass rather than the LP strategy - a
%                    separation the FASTERCC paper never makes, because it
%                    never applies its preprocessing to any baseline.
%   spectracc_dir    spectraCCdir: forward and reverse tracked separately, so
%                    it classifies forward / reverse / reversible / blocked.
%                    The SPECTRA analogue of fasterCC revopt=1.
%   fastercc_r0      FASTERCC, consistency only. The like-for-like baseline.
%   fastercc_r1      FASTERCC with its reversibility pass. The only fair
%                    comparator for spectracc_dir, since both return
%                    directionality and r0 does not.
%
% spectracc_prep is charged the preprocessing time as well as the solve.
% Charging only the solve would make the preprocessing look free.
%
% GROUND TRUTH - two strategies, by regime
%   small   full FVA. At ~10,600 reactions that is ~21,200 LPs, affordable,
%           and it scores every reaction including direction labels.
%   large   TARGETED FVA. A full FVA at ~184k reactions would be ~368k LPs,
%           which is not affordable and would dwarf the benchmark it is meant
%           to check. Instead: run the arms, find the reactions they actually
%           DISAGREE about, and FVA only those. Agreement between independent
%           methods is weak evidence but disagreement is exactly where truth
%           is needed, so this spends the LP budget where it decides something.
%
% Two error modes are reported SEPARATELY and must never be summed:
%   overclaim  - called feasible, FVA says blocked. UNSOUND, a real bug.
%   underclaim - called blocked, FVA says feasible. A miss: the batch loop
%                stops at the first unproductive LP with no one-by-one
%                fallback. Fixable via opts.confirmBlocked, at a cost in LPs.

if nargin < 1 || isempty(cfg),    cfg = bench_config('smoke'); end
if nargin < 2 || isempty(regime), regime = 'small'; end
if nargin < 3 || isempty(order),  order  = 'spectra_first'; end
if nargin < 4 || isempty(opts),   opts   = struct(); end
if ~isfield(opts,'fva'),  opts.fva  = true; end   % false = runtime only
if ~isfield(opts,'reps'),     opts.reps     = cfg.timingReps; end
if ~isfield(opts,'instReps'), opts.instReps = cfg.instReps;   end

% opts.repOnly pins a SINGLE draw, so a SLURM array can put one instance in
% one task. Without it every task would regenerate the whole draw set and
% they would duplicate each other's work. The draw number still feeds
% bench_seed, so instance k20_rm05000_rep02 is the same model whether it was
% built as task 2 of an array or as the second iteration of a loop.
if ~isfield(opts,'repOnly'), opts.repOnly = []; end

% TWO KINDS OF REPLICATE, and they measure different things.
%   opts.reps      re-times the SAME model. Captures machine noise only -
%                  contention, paging, JIT. Measured at up to 2.5x.
%   opts.instReps  draws a DIFFERENT random set of reactions to delete at the
%                  same removal level. Captures how much the answer depends on
%                  WHICH reactions are blocked, which is the scientific
%                  question and the one the FASTERCC paper averages over with
%                  its 10 replicates per condition.
% Reporting only the first would attribute instance-to-instance spread to the
% machine, or vice versa. They are nested, so both are kept.

% FVA is ~90%% of wall time and is DETERMINISTIC, so it never needs repeating.
% Dropping it makes replicates cheap, which matters: measured run-to-run
% variance reached 2.5x on bit-identical seeded LPs (k02_rm500: 32.85s vs
% 12.91s), which is larger than most inter-arm gaps. Without replicates no
% runtime difference under ~2x is reportable.

tag = [regime '_' order];
if ~opts.fva,      tag = [tag '_nofva']; end
if opts.reps > 1,  tag = [tag sprintf('_rep%d', opts.reps)]; end
if opts.instReps > 1, tag = [tag sprintf('_x%d', opts.instReps)]; end
% opts.tagSuffix separates concurrent SLURM array tasks. Without it every
% task of an array writes to one results directory and they overwrite each
% other's checkpoint.
if isfield(opts,'tagSuffix') && ~isempty(opts.tagSuffix)
    tag = [tag '_' opts.tagSuffix];
end
changeCobraSolver(cfg.solver, 'all', 0);

% Tighten the solver for EVERY arm. 1e-9 is Gurobi's hard floor for
% FeasibilityTol, so this is as strict as the solver allows.
%
% Why it matters here and not only in the community runs: the one unsound
% result found so far - spectraCC calling four provably blocked reactions
% consistent at rm=1020 - is consistent with slop being read as flux. A
% looser feasibility tolerance lets an LP return a vertex that violates
% S*v = 0 by more than the detection cutoff, and the reaction then looks
% like it carries flux. Setting this identically for every arm means no arm
% is advantaged by a different numerical regime.
%
% opts.feasTol / opts.epsilon override the configured numerics for the
% tolerance sweeps, and nothing else. A swept run therefore differs from the
% standing configuration in exactly the one number being swept - the instance
% is seeded from (size, rep) alone, so the model handed to the arm is byte
% for byte the model the 1e-9 run saw.
if isfield(opts,'feasTol') && ~isempty(opts.feasTol)
    cfg.feasTol = opts.feasTol;
    cfg.optTol  = opts.feasTol;
end
changeCobraSolverParams('LP', 'feasTol', cfg.feasTol);
changeCobraSolverParams('LP', 'optTol',  cfg.optTol);

eps_ = cfg.epsilon;
if isfield(opts,'epsilon') && ~isempty(opts.epsilon), eps_ = opts.epsilon; end

% Guard the threshold against the tolerance. spectraCC counts a reaction as
% carrying flux at 0.99*epsilon; if that sits too close to feasTol, solver
% noise is indistinguishable from real flux. 100x is the margin the default
% epsilon=1e-4 enjoys against feasTol=1e-9 by a wide margin; epsilon=1e-8
% leaves only ~10x, which is thin but deliberate for the community runs.
if 0.99*eps_ < 10*cfg.feasTol
    warning('bench_variants:thinMargin', ...
        ['Detection cutoff %.3g is within 10x of feasTol %.3g. Fluxes at ' ...
         'that scale cannot be distinguished from solver noise.'], ...
        0.99*eps_, cfg.feasTol);
end

%% Instance list - no replicates --------------------------------------
switch lower(regime)
    case 'small'
        S = load(fullfile(cfg.baseDir,'consRecon3D.mat'));
        base     = S.consRecon3D;
        removals = [0 10 50 100 200 500];
        fvaMode  = 'full';
        inst = {};
        for r = removals
            inst{end+1} = struct('name',sprintf('recon_rm%03d',r), ...
                                 'kind','small','k',0,'nRemove',r,'rep',1); %#ok<AGROW>
        end

    case 'recon'
        % Recon3D alone, replicated. The paper-style design: fix the model,
        % vary WHICH reactions are removed, and average over draws.
        S = load(fullfile(cfg.baseDir,'consRecon3D.mat'));
        base     = S.consRecon3D;

        % A SWEEP, not a few points. 21 levels from 0 to 10,000 deleted out of
        % 10,600, so the surviving model runs from fully intact down to ~600
        % reactions. This is the axis the FASTERCC paper's speedup claim
        % actually moves along (2.3-5.2x at 0 blocked, 17.5-31.5x at 500), and
        % it only sampled up to 500. Carrying it to 10,000 shows where, or
        % whether, the advantage turns over - and the measured answer so far is
        % that it reverses once past a few hundred blocked reactions.
        removals = 0:500:10000;

        % opts.levels narrows the sweep, so a pilot can size the full job
        % using the SAME code path rather than a separate script that might
        % drift from it.
        if isfield(opts,'levels') && ~isempty(opts.levels)
            removals = opts.levels;
        end

        % This regime is intended to run with opts.fva = false, which forces
        % fvaMode to 'none' below. Disagreements are still DETECTED and
        % written out in full (res.disagreeRxns / res.disagreeFlags), they
        % are simply not adjudicated here - that is a separate, later pass on
        % the stored reaction lists, so no FVA is ever paid for during a
        % timing run. 'targeted' is the fallback if someone does pass
        % opts.fva = true; it would fire only on disagreeing reactions, never
        % the full model, which at ~65 min per instance would be ~55 h here.
        fvaMode  = 'targeted';

        inst = {};
        for r = removals
            % rm000 deletes nothing, so every "draw" would be the SAME model.
            % Replicating it would report machine noise as if it were
            % draw-to-draw spread and make the control look falsely tight.
            % One instance; opts.reps still re-times it.
            nRep = opts.instReps; if r == 0, nRep = 1; end
            for rp = 1:nRep
                inst{end+1} = struct( ...
                    'name',sprintf('recon_rm%04d_rep%02d',r,rp), ...
                    'kind','small','k',0,'nRemove',r,'rep',rp); %#ok<AGROW>
            end
        end

    case 'large'
        base     = [];
        ks       = [2 5 10];
        removals = [0 500];
        fvaMode  = 'targeted';

        % Expanding a cluster model is expensive and cached per k. Only k=1,2
        % exist from the smoke run, so make sure k=5 and k=10 are built before
        % the sweep rather than failing partway through.
        gcfg = cfg;
        gcfg.stage         = 'variants_large';
        gcfg.clusterCounts = ks;
        gcfg.removalCounts = removals;
        gcfg.removalReps   = 1;
        bench_generate_models(gcfg);

        inst = {};
        for k = ks
            for r = removals
                inst{end+1} = struct('name',sprintf('k%02d_rm%03d',k,r), ...
                                     'kind','large','k',k,'nRemove',r,'rep',1); %#ok<AGROW>
            end
        end

    case 'papersizes'
        % Experiment 1, as the SPECTRA paper plots it.
        %
        % ForSPECTRA_CCResults.m sweeps n = 20:20:10600 and KEEPS n randomly
        % chosen reactions:
        %     rxn_ids = sort(randsample(numel(model.rxns), n(i)));
        %     m = removeRxns(model, model.rxns(setdiff(1:end, rxn_ids)));
        % but Runtime_SPECTRA_CC.ipynb plots only 10 of those 530 points:
        %     ids = [int(i) for i in np.linspace(0, 528, 10)]
        % which lands on the sizes below. Matching the published figure
        % exactly means these timings drop straight into it, with fastercc_r0
        % added as the fourth bar.
        %
        % Note the axis is KEPT reactions, not removed ones, so it runs the
        % opposite way to the 'recon' sweep: 20 is the tiny, almost entirely
        % blocked end and 10580 is the near-intact model.
        S = load(fullfile(cfg.baseDir,'consRecon3D.mat'));
        base     = S.consRecon3D;
        sizes    = [20 1180 2360 3540 4700 5880 7060 8220 9400 10580];
        if isfield(opts,'levels') && ~isempty(opts.levels), sizes = opts.levels; end
        fvaMode  = 'targeted';

        inst = {};
        for ns = sizes
            for rp = 1:opts.instReps
                inst{end+1} = struct( ...
                    'name',sprintf('n%05d_rep%02d',ns,rp), ...
                    'kind','keepn','k',0,'nRemove',numel(base.rxns)-ns, ...
                    'nKeep',ns,'rep',rp); %#ok<AGROW>
            end
        end

    case 'expanded'
        % Experiment 2: the cluster-size axis, k = 5 .. 50.
        %
        % DRAWS TAPER WITH k, and the reason is arithmetic, not taste. The
        % measured scaling exponent from k=2 -> k=5 is ~2.54 (2.28x the
        % reactions, 8.1x the time), so one k=50 instance is ~8 h per arm.
        % Three draws x four removal levels x two arms at that size is ~190 h,
        % which no wall allows. Spending the draws where they are affordable
        % and taking single points at the top is the only design that reaches
        % k=50 at all; the top two rows carry no error bar and must be
        % reported as single measurements.
        base     = [];
        ks       = [5 10 20 30 50];
        % Draws per k. Built AFTER any opts.ks override so a pilot asking for
        % an unlisted k does not index a missing key; anything not named here
        % falls back to a single draw, which is the safe direction - it
        % under-replicates rather than overrunning the wall.
        bigK     = [30 50];
        drawsFor = @(k) (opts.instReps * ~ismember(k, bigK)) + ismember(k, bigK);
        removals = [500 2000 5000 10000];
        if isfield(opts,'ks')     && ~isempty(opts.ks),     ks       = opts.ks;     end
        if isfield(opts,'levels') && ~isempty(opts.levels), removals = opts.levels; end
        fvaMode  = 'targeted';

        % bench_generate_models builds one expanded model per k; the removal
        % counts reach it only through the manifest it writes, and it seeds
        % each manifest row with bench_seed, which rejects any level outside
        % cfg.allRemovalCounts. Hand it the canonical list so generation
        % cannot fail on a level it never actually uses.
        gcfg = cfg;
        gcfg.stage         = 'variants_expanded';
        gcfg.clusterCounts = ks;
        gcfg.removalCounts = cfg.allRemovalCounts;
        gcfg.removalReps   = 1;
        bench_generate_models(gcfg);

        inst = {};
        for k = ks
            for r = removals
                reps = 1:drawsFor(k);
                if ~isempty(opts.repOnly), reps = opts.repOnly; end
                for rp = reps
                    inst{end+1} = struct( ...
                        'name',sprintf('k%02d_rm%05d_rep%02d',k,r,rp), ...
                        'kind','large','k',k,'nRemove',r,'rep',rp); %#ok<AGROW>
                end
            end
        end

    case 'community'
        % Experiment 3: AGORA2 communities under random media. Instance
        % construction lives in bench_make_community, which also writes the
        % media it drew, so an instance is reproducible from its seed alone.
        base     = [];
        sizes    = [50 100 200];
        if isfield(opts,'ks') && ~isempty(opts.ks), sizes = opts.ks; end
        fvaMode  = 'targeted';
        inst = {};
        for ns = sizes
            reps = 1:opts.instReps;
            if ~isempty(opts.repOnly), reps = opts.repOnly; end
            for rp = reps
                inst{end+1} = struct( ...
                    'name',sprintf('com%03d_rep%02d',ns,rp), ...
                    'kind','community','k',ns,'nRemove',0,'rep',rp); %#ok<AGROW>
            end
        end

    otherwise
        error('bench_variants:regime', ...
              ['regime must be one of: small, recon, papersizes, large, ' ...
               'expanded, community (got "%s")'], regime);
end

% ARM SET, per regime.
%
%   recon      spectracc, fastercc_r0, fastcc, swiftcc. The full field, run
%              only here because this is the only regime small enough that
%              fastcc and swiftcc finish - at 10,600 reactions. swiftcc's
%              QR factorisation is also the memory risk, and it is bounded
%              here in a way it would not be at k=50.
%   expanded,  spectracc vs fastercc_r0 only. The directional arms were
%   community  dropped once full FVA showed them adding nothing as
%              consistency checkers, and spectracc_prep was dropped once
%              three independent runs showed FASTERCC's structural
%              preprocessing gives SPECTRA no measurable speedup
%              (56.7 s vs 54.7 s over a 20-level sweep).
%
% Pass opts.arms to override; the directional arms and spectracc_prep still
% exist and can be re-run, they are simply not on the default path.
if isfield(opts,'arms') && ~isempty(opts.arms)
    arms = opts.arms;
else
    switch lower(regime)
        case {'recon','papersizes'}
            arms = {'spectracc','fastercc_r0','fastcc','swiftcc'};
        case {'expanded','community'}
            arms = {'spectracc','fastercc_r0'};
        otherwise
            arms = {'spectracc','spectracc_prep','fastercc_r0'};
    end
end
if strcmpi(order,'fastercc_first')
    isF = startsWith(arms, {'fastercc','fastcc','swiftcc'});
    arms = [arms(isF), arms(~isF)];
end

% The arm count goes in the tag. Without it a 3-arm run would land in the
% same directory as the completed 5-arm runs and resume from their
% checkpoint, silently producing a table whose rows have different arm sets.
tag = [tag sprintf('_a%d', numel(arms))];
outDir = fullfile(cfg.root, 'results', ['variants_' tag]);
if ~exist(outDir,'dir'), mkdir(outDir); end
ckpt = fullfile(outDir, 'checkpoint.mat');

fprintf('Arm order (%s): %s\n', order, strjoin(arms, ' -> '));
if ~opts.fva
    fvaMode = 'none';
    fprintf('FVA DISABLED - runtime only, %d replicate(s) per arm\n', opts.reps);
end

%% Resume -----------------------------------------------------------------
if exist(ckpt,'file')
    C = load(ckpt); R = C.R; startIdx = C.nextIdx;
    fprintf('Resuming at instance %d/%d\n', startIdx, numel(inst));
else
    R = {}; startIdx = 1;
end

%% Sweep ------------------------------------------------------------------
for ii = startIdx:numel(inst)
    I = inst{ii};

    if strcmp(I.kind,'small')
        % Seed from (nRemove, rep), NOT from the position ii. The old
        % positional seed made an instance's identity depend on the length and
        % order of the instance list, so adding a removal level or a replicate
        % silently redrew every other instance, and recon_rm050 in one run was
        % a different model from recon_rm050 in the next.
        %
        % Deliberately NOT bench_seed: that function indexes nRemove into
        % cfg.allRemovalCounts and errors on anything outside it, which rules
        % out the 20:500:10000 sweep. Using nRemove directly needs no list, so
        % any removal level is addressable. Multiplying by 1000 keeps draws
        % separated, since rep only ever reaches instReps, so (nRemove, rep)
        % maps to a distinct seed for every instance in both regimes that use
        % this branch.
        seed  = cfg.baseSeed + 1000*I.nRemove + I.rep;
        model = base;
        removedNames = {};
        if I.nRemove > 0
            rng(seed,'twister');
            drop = sort(randsample(numel(model.rxns), I.nRemove));
            removedNames = model.rxns(drop);
            model = removeRxns(model, removedNames);
        end
    elseif strcmp(I.kind,'keepn')
        % Keep nKeep reactions, mirroring ForSPECTRA_CCResults.m: it samples
        % the reactions to KEEP and deletes the complement. Sampling the
        % complement instead would draw a different subset for the same seed,
        % so the sampling call is kept the way the paper ran it.
        seed = cfg.baseSeed + 1000*I.nKeep + I.rep;
        rng(seed,'twister');
        keep  = sort(randsample(numel(base.rxns), I.nKeep));
        drop  = setdiff((1:numel(base.rxns))', keep);
        removedNames = base.rxns(drop);
        model = removeRxns(base, removedNames);

    elseif strcmp(I.kind,'community')
        [model, mi]  = bench_make_community(I.k, I.rep, cfg);
        seed         = mi.seed;
        removedNames = {};        % nothing is deleted; the medium does the work
        res_extra    = mi;        % species list, medium stats - stored below
    else
        [model, mi]  = bench_make_instance(I.k, I.nRemove, I.rep, cfg);
        seed         = mi.seed;
        removedNames = mi.removed;
    end
    n = numel(model.rxns);
    fprintf('\n=== [%d/%d] %s : %d rxns, %d mets ===\n', ...
            ii, numel(inst), I.name, n, numel(model.mets));

    res = struct('name',I.name,'kind',I.kind,'k',I.k,'nRemove',I.nRemove, ...
                 'seed',seed,'nRxns',n,'nMets',numel(model.mets));

    %% arms ---------------------------------------------------------------
    % PARTIAL RESULTS ARE WRITTEN AFTER EVERY ARM, and a marker file is laid
    % down before each one starts.
    %
    % An OOM kill is SIGKILL: MATLAB cannot catch it, run a cleanup, or record
    % anything. Saving only at the end of an instance therefore lost
    % everything when an arm was killed - the k=20 run had spectracc finish in
    % ~26 min and then died in fastercc_r0, and nothing at all was kept.
    %
    % With this, a killed arm leaves its marker behind and every earlier arm
    % is already on disk, so "fastercc_r0 could not complete in 64 GB at
    % k=20" becomes a recorded result rather than a lost run.
    instFile = fullfile(outDir, [I.name '.mat']);
    feas = false(n, numel(arms));     % membership matrix, for disagreement
    clsAll = cell(numel(arms),1);
    for a = 1:numel(arms)
        arm = arms{a};
        marker = fullfile(outDir, sprintf('%s.%s.running', I.name, arm));
        mfid = fopen(marker, 'w');
        if mfid > 0
            fprintf(mfid, '%s started %s\n', arm, datestr(now, 31)); %#ok<TNOW1,DATST>
            fclose(mfid);
        end
        cons = []; cls = []; lps = NaN; t = NaN; err = '';
        times = nan(1, opts.reps);
        for rp = 1:opts.reps
        try
            switch arm
                % spectraCC is called WITHOUT a seed, so it uses its own
                % default RNG behaviour exactly as a user would get it.
                % Pinning the seed would have made SPECTRA the only arm whose
                % randomness was controlled while fasterCC ran free - and
                % fasterCC has now been shown to be nondeterministic (same
                % model, LP counts of 14, 9, 13). Seeding one side and not the
                % other is not a like-for-like comparison. The instance itself
                % is still deterministic: the reaction DELETION is seeded
                % above, so every arm sees the identical model.
                case 'spectracc'
                    tic; [cons,lps] = spectraCC(model, eps_); t = toc;
                case 'spectracc_prep'
                    tic;
                    [prep, map] = prep_fastercc(model);
                    [cP, lps]   = spectraCC(prep, eps_);
                    cons = map(cP);
                    t = toc;
                case 'spectracc_dir'
                    tic; [cls,iD] = spectraCCdir(model, eps_, 'stoichiometry'); t = toc;
                    cons = iD.consistent; lps = iD.LPS;
                case 'fastercc_r0'
                    tic; [A,~,nLP] = fasterCC(model, eps_, 0, 0); t = toc;
                    cons = A(:); lps = nLP.total;

                % FASTCC, the original. COBRA's copy is patched to return its
                % LP count as a 4th output - counter increments only, no
                % behaviour change.
                case 'fastcc'
                    tic; [A,~,~,nl] = fastcc(model, eps_, 0); t = toc;
                    cons = A(:); lps = nl;

                % SWIFTCC. Takes (S, rev), not a model.
                %
                % rev is model.lb < 0, matching ForSPECTRA_CCResults.m:5 -
                % the convention the SPECTRA paper's own runtime comparison
                % used. It is NOT the same as (lb < 0 & ub > 0): a
                % reverse-only reaction counts as reversible here. Changing
                % it would make these timings incomparable with the published
                % ones, so it stays as the paper ran it.
                %
                % LP count is 1, also as ForSPECTRA_CCResults.m:39 records it
                % (ones(530,1)). It is literally true - swiftcc solves one LP,
                % blocked(), for the irreversible part - but it is not
                % comparable with the other arms, because the reversible part
                % is settled by a QR factorisation of S' rather than by LPs.
                % Compare this arm on RUNTIME; its LP column means something
                % different from everyone else's.
                case 'swiftcc'
                    rev = model.lb < 0;
                    tic; sc = swiftcc(model.S, rev, cfg.solver); t = toc;
                    cons = sc(:); lps = 1;
                case 'fastercc_r1'
                    tic; [A,mOut,nLP] = fasterCC(model, eps_, 0, 1); t = toc;
                    cons = A(:); lps = nLP.total;
                    cls = repmat(categorical({'blocked'}), n, 1);
                    for q = 1:numel(A)
                        lo = mOut.lb(q); hi = mOut.ub(q);
                        if     lo < -eps_ && hi > eps_, cls(A(q)) = categorical({'reversible'});
                        elseif hi >  eps_,              cls(A(q)) = categorical({'forward'});
                        elseif lo < -eps_,              cls(A(q)) = categorical({'reverse'});
                        end
                    end
            end
            times(rp) = t;
        catch ME
            err = sprintf('%s: %s', ME.identifier, ME.message);
            fprintf(2,'  %-16s FAILED: %s\n', arm, ME.message);
            break
        end
        end
        if isempty(err)
            if opts.reps > 1
                fprintf('  %-16s %9.2f s  [min %6.2f max %6.2f of %d]  |cons|=%7d  LPs=%s\n', ...
                        arm, median(times,'omitnan'), min(times), max(times), ...
                        opts.reps, numel(cons), num2str(lps));
            else
                fprintf('  %-16s %9.2f s  |cons|=%7d  LPs=%s\n', ...
                        arm, t, numel(cons), num2str(lps));
            end
        end
        feas(cons, a) = true;
        clsAll{a} = cls;
        % Report the MEDIAN across replicates, and keep the raw vector: with
        % 2.5x observed variance, a single draw is not a measurement.
        res.(arm) = struct('time',median(times,'omitnan'),'times',times, ...
                           'tmin',min(times),'tmax',max(times), ...
                           'lps',lps,'nCons',numel(cons),'error',err, ...
                           'consRxns',{model.rxns(cons)});
        if ~isempty(cls), res.(arm).cls = cls; end

        % Persist what we have, then clear the marker. Order matters: if the
        % save is interrupted the marker stays and the arm reads as killed,
        % which is the safe direction.
        res.armsAttempted = arms(1:a);
        save(instFile, 'res', '-v7.3');
        if exist(marker, 'file'), delete(marker); end
    end

    %% ground truth -------------------------------------------------------
    % Which reactions do the arms actually disagree about?
    ran      = ~cellfun(@(a) ~isempty(res.(a).error), arms);
    disagree = find(any(feas(:,ran) ~= feas(:,find(ran,1)), 2));
    % Plus direction disagreements, but ONLY if both directional arms are in
    % this run. They are off by default now, so find() returns empty and
    % indexing clsAll with it would error - hence the isscalar guards.
    iDir = find(strcmp(arms,'spectracc_dir'));
    iR1  = find(strcmp(arms,'fastercc_r1'));
    if isscalar(iDir) && isscalar(iR1) && ...
            ~isempty(clsAll{iDir}) && ~isempty(clsAll{iR1})
        disagree = union(disagree, find(clsAll{iDir} ~= clsAll{iR1}));
    end
    fprintf('  arms disagree on %d / %d reactions\n', numel(disagree), n);
    res.nDisagree = numel(disagree);

    % PERSIST THE DISAGREEMENT ITSELF, not just its size. With FVA switched
    % off there is no adjudication in this run, so the disagreement is the
    % whole scientific output and a bare count cannot be followed up - you
    % could see "17 reactions" and have no way to learn which, or which arm
    % was the odd one out, without re-running the sweep.
    %
    % Names, not indices: removal shifts indices, so an index means nothing
    % outside the instance that produced it.
    res.arms          = arms;

    % The numerics this instance actually ran under. The tag carries them too,
    % but a tag is a directory name and gets renamed; this travels with the
    % result, so a pooled table cannot silently mix tolerances.
    res.epsilon = eps_;
    res.feasTol = cfg.feasTol;

    res.disagreeRxns  = model.rxns(disagree);
    res.disagreeFlags = feas(disagree, :);   % rows: reactions, cols: arms

    % Full per-arm membership, as a logical mask rather than a list of names.
    % n x nArms logical is ~32 KB per instance here against ~2 MB for the
    % name lists, and loses nothing: res.rxns below gives the key. This is
    % what makes any later re-analysis possible without recomputation.
    res.rxns = model.rxns;
    res.feas = feas;

    % Which reactions were deleted to build this instance. Redundant with the
    % seed, which regenerates them, but a stored list needs no rerun and no
    % matching RNG version to interpret.
    res.removedRxns = removedNames;

    % Community instances carry their species list and the medium actually
    % drawn. Without these the instance is not reconstructible from the file
    % alone, and min/maxOpenInflow is what makes the tolerance margin
    % auditable rather than assumed.
    if strcmp(I.kind,'community') && exist('res_extra','var')
        res.community = res_extra;
    end

    switch fvaMode
        case 'none'
            minF = []; maxF = []; idxF = []; tF = 0;
            fprintf('  FVA skipped (runtime-only run)\n');
        case 'full'
            fprintf('  %-16s ', 'FVA (full)');
            t0 = tic; [minF,maxF] = fluxVariability(model, 0, 'max'); tF = toc(t0);
            idxF = (1:n)';
        case 'targeted'
            if isempty(disagree)
                fprintf('  FVA skipped: arms agree everywhere\n');
                minF = []; maxF = []; idxF = []; tF = 0;
            else
                fprintf('  %-16s ', sprintf('FVA (%d rxns)', numel(disagree)));
                t0 = tic;
                [minF,maxF] = fluxVariability(model, 0, 'max', model.rxns(disagree));
                tF = toc(t0);
                idxF = disagree(:);
            end
    end

    if ~isempty(idxF)
        fprintf('%8.2f s\n', tF);
        fwdT = maxF >  eps_;
        revT = minF < -eps_;
        clsT = repmat(categorical({'blocked'}), numel(idxF), 1);
        clsT(fwdT & ~revT) = categorical({'forward'});
        clsT(~fwdT & revT) = categorical({'reverse'});
        clsT(fwdT &  revT) = categorical({'reversible'});
        feasT = fwdT | revT;

        res.fva = struct('time',tF,'idx',idxF,'min',minF,'max',maxF, ...
                         'cls',clsT,'feasible',feasT, ...
                         'mode',fvaMode,'nScored',numel(idxF));

        for a = 1:numel(arms)
            if ~isempty(res.(arms{a}).error), continue; end
            fp = feas(idxF, a);
            res.(arms{a}).overclaim  = sum(fp & ~feasT);
            res.(arms{a}).underclaim = sum(~fp &  feasT);
            if ~isempty(clsAll{a})
                both = feasT & fp;
                res.(arms{a}).dirMismatch = sum(clsAll{a}(idxF(both)) ~= clsT(both));
            end
        end
    end

    R{end+1} = res; %#ok<AGROW>
    nextIdx = ii + 1; %#ok<NASGU>
    save(ckpt, 'R', 'nextIdx', '-v7.3');
    save(instFile, 'res', '-v7.3');
end

%% Table ------------------------------------------------------------------
rows = {};
for ii = 1:numel(R)
    r = R{ii};
    for a = 1:numel(arms)
        if ~isfield(r, arms{a}), continue; end
        e = r.(arms{a});
        rows(end+1,:) = { r.name, arms{a}, r.kind, r.k, r.nRemove, r.nRxns, ...
                          e.time, getf(e,'tmin'), getf(e,'tmax'), numel(e.times), ...
                          e.lps, e.nCons, r.nDisagree, ...
                          getf(e,'overclaim'), getf(e,'underclaim'), ...
                          getf(e,'dirMismatch'), string(e.error) }; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', ...
    {'instance','arm','kind','clusters','nRemoved','nRxns', ...
     'runtimeSec','runtimeMin','runtimeMax','nReps','nLPs', ...
     'nConsistent','nDisagree','overclaim','underclaim','dirMismatch','error'});
csv = fullfile(cfg.root, sprintf('variants_%s.csv', tag));
writetable(T, csv);
fprintf('\nWrote %s (%d rows)\n', csv, height(T));

summarise(T, arms);
end

% ------------------------------------------------------------------ helpers
function v = getf(s, f)
if isfield(s,f) && ~isempty(s.(f)), v = double(s.(f)); else, v = NaN; end
end

function summarise(T, arms)
ok     = T.error == "";
levels = unique(T.nRemoved)';
hdr = sprintf('%-18s', 'arm');
for L = levels, hdr = [hdr sprintf('%11s', sprintf('rm=%d',L))]; end %#ok<AGROW>

fprintf('\n=== runtime (s), by blocked-reaction count ===\n'); disp(hdr);
for a = 1:numel(arms)
    line = sprintf('%-18s', arms{a});
    for L = levels
        v = T.runtimeSec(strcmp(T.arm,arms{a}) & T.nRemoved==L & ok);
        line = [line sprintf('%11s', fmt(v))]; %#ok<AGROW>
    end
    disp(line);
end

% Spread, printed only when there is actually more than one measurement per
% cell. A bare median hid the finding that mattered on the k=2 runs: the arms
% had near-identical MINIMA and very different medians, i.e. they differed in
% tail behaviour under contention, not in throughput. n is printed so a reader
% can see how much the median is worth.
nPer = arrayfun(@(L) max(arrayfun(@(a) ...
          sum(strcmp(T.arm,arms{a}) & T.nRemoved==L & ok), 1:numel(arms))), levels);
if any(nPer > 1)
    fprintf('\n=== runtime spread (s) ===\n');
    fprintf('%-18s %7s %4s %9s %9s %9s\n','arm','rm','n','median','min','max');
    for a = 1:numel(arms)
        for L = levels
            v = T.runtimeSec(strcmp(T.arm,arms{a}) & T.nRemoved==L & ok);
            v = v(~isnan(v));
            if isempty(v), continue; end
            fprintf('%-18s %7d %4d %9.2f %9.2f %9.2f\n', ...
                    arms{a}, L, numel(v), median(v), min(v), max(v));
        end
    end
end

fprintf('\n=== LPs solved ===\n'); disp(hdr);
for a = 1:numel(arms)
    line = sprintf('%-18s', arms{a});
    for L = levels
        v = T.nLPs(strcmp(T.arm,arms{a}) & T.nRemoved==L & ok);
        line = [line sprintf('%11s', fmt(v(~isnan(v))))]; %#ok<AGROW>
    end
    disp(line);
end

fprintf('\n=== accuracy vs FVA (totals; FVA scores all reactions in the\n');
fprintf('    small regime, only the disputed ones in the large regime) ===\n');
fprintf('%-18s %12s %12s %12s\n','arm','overclaim','underclaim','dirMismatch');
fprintf('%-18s %12s %12s %12s\n','','(UNSOUND)','(missed)','');
for a = 1:numel(arms)
    m = strcmp(T.arm,arms{a}) & ok;
    fprintf('%-18s %12s %12s %12s\n', arms{a}, ...
        fmt0(sum(T.overclaim(m),'omitnan')), ...
        fmt0(sum(T.underclaim(m),'omitnan')), ...
        fmt0(sum(T.dirMismatch(m),'omitnan')));
end
end

function s = fmt(x)
x = x(~isnan(x));
if isempty(x), s = '-'; else, s = sprintf('%.2f', median(x)); end
end
function s = fmt0(x)
if isnan(x), s = '-'; else, s = sprintf('%d', round(x)); end
end
