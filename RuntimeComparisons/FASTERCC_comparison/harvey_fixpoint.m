function harvey_fixpoint(arm, cfg)
% harvey_fixpoint  Is consistency checking idempotent on Harvey?
%
%   harvey_fixpoint('spectracc')
%   harvey_fixpoint('spectracc_relaxed')
%   harvey_fixpoint('fastercc_r0')
%
% Runs one consistency checker on the Harvey whole-body model, builds a new
% model from the reactions it called consistent, and runs it again - until
% the answer stops changing or the iteration cap is hit.
%
% WHY THIS IS A TEST AND NOT A FORMALITY
% A flux-consistent subnetwork is by definition entirely flux consistent, so
% the second pass MUST return every reaction it is given. Any reaction the
% second pass drops was reported consistent by the first pass and is not,
% which is an overclaim - the error class that matters, because a model
% built on it carries dead reactions the checker certified. Equally, a
% second pass that keeps shrinking for several rounds says the method
% converges rather than decides, and the round count is the cost of
% trusting it. This is the strongest soundness check available without FVA,
% and it needs no ground truth.
%
% WHY HARVEY
% Every model benchmarked so far admits the zero flux vector. Harvey does
% not: 167 of its 81,094 reactions carry bounds that exclude zero - 88
% Diet_EX_ uptakes forced by setDietConstraints, 49 Kidney_EX_ blood
% exchanges and 12 BBB_ transports from physiologicalConstraintsHMDBbased,
% plus water/CO2/urea excretion, Heart_DM_atp_c_ and the whole-body
% objective pinned at 1. That makes it the only model here where "relax the
% bounds to admit zero" is not a no-op, which is what the second arm tests.
%
% THE ARMS
%   spectracc          spectraCC on the model exactly as distributed.
%   spectracc_relaxed  spectraCC after lb(lb>0)=0 and ub(ub<0)=0, applied
%                      ONCE at the start. This is precisely what fasterCC
%                      does to its own input (fasterCC.m:40-46) before it
%                      does anything else, restoring the bounds at the end
%                      through restoreBounds. So the relaxed arm is SPECTRA
%                      answering the question FASTERCC actually answers, and
%                      the gap between the two SPECTRA arms is the price of
%                      the forced fluxes rather than an algorithmic
%                      difference. SPECTRA's own code is untouched - only
%                      its input differs.
%   fastercc_r0        fasterCC(model, epsilon, 0, 0), the like-for-like arm.
%
% NUMERICS: epsilon = 1e-4. SOLVER TOLERANCES ARE NOT SET - the solver runs
% at whatever default it ships with, which is the point of this run and the
% one way it differs from the rest of the benchmark, where feasTol is pinned
% to 1e-9. The effective value is read back with getCobraSolverParams and
% recorded in every saved file, so the run is reproducible without this
% script having chosen it. Read it from the output rather than assuming it:
% the community sweeps showed SPECTRA's answer is only reproducible at
% 1e-9, so a set that keeps shrinking here may be the tolerance rather than
% the method, and the two are not separable from this run alone.
%
% Results are written after EVERY iteration, so a wall-clock kill keeps
% everything already computed. That matters: an iteration here is tens of
% minutes and the round count is not known in advance.

if nargin < 2 || isempty(cfg), cfg = bench_config('smoke'); end

EPSILON = 1e-4;
MAXITER = 12;      % a cap, not an expectation; 2 is the sound answer

outDir = fullfile(cfg.root, 'results', 'harvey_fixpoint', arm);
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% solver ----------------------------------------------------------------
% Solver selected, tolerances DELIBERATELY NOT TOUCHED. No
% changeCobraSolverParams call anywhere in this file.
changeCobraSolver(cfg.solver, 'all', 0);
[feasTolEff, optTolEff] = getCobraSolverParams('LP', {'feasTol', 'optTol'});
fprintf('arm=%s  epsilon=%.0e  solver=%s\n', arm, EPSILON, cfg.solver);
fprintf('solver tolerances left at their defaults: feasTol=%g optTol=%g\n\n', ...
        feasTolEff, optTolEff);

%% model -----------------------------------------------------------------
% loadPSCMfile resolves the newest Harvey on the path. The file it finds is
% byte-identical (md5 5d91cda703b0f9bfafed14cbf61371fd) to the copy under
% test/models/mat in the COBRA Toolbox checkout, so this is the distributed
% model, not a local variant.
model = loadPSCMfile('Harvey');
checkCouplingWidth(model);

n0 = numel(model.rxns);
fprintf('Harvey loaded: %d rxns, %d mets\n', n0, numel(model.mets));
fprintf('  lb > 0 : %d\n  ub < 0 : %d\n', ...
        full(sum(model.lb > 0)), full(sum(model.ub < 0)));

if strcmp(arm, 'spectracc_relaxed')
    % The relaxation, applied once. Later iterations inherit it through
    % removeRxns, which preserves bounds, so it is not reapplied.
    nRelax = full(sum(model.lb > 0)) + full(sum(model.ub < 0));
    model.lb(model.lb > 0) = 0;
    model.ub(model.ub < 0) = 0;
    fprintf('  RELAXED %d bounds so the zero vector is feasible\n', nRelax);
end
fprintf('\n');

%% iterate ---------------------------------------------------------------
R = struct('iter', {}, 'nRxnsIn', {}, 'nMetsIn', {}, 'nCons', {}, ...
           'nDropped', {}, 'lps', {}, 'time', {}, 'stable', {}, 'err', {});

stable = false;
failed = false;
it     = 0;
for it = 1:MAXITER
    nIn = numel(model.rxns);
    fprintf('--- iteration %d: %d rxns, %d mets ---\n', it, nIn, numel(model.mets));

    cons = []; lps = NaN; t = NaN; errMsg = '';
    try
        switch arm
            % spectraCC is called UNSEEDED, as everywhere else in this
            % benchmark: seeding one side and not fasterCC would not be a
            % like-for-like comparison.
            case {'spectracc', 'spectracc_relaxed'}
                tic; [cons, lps] = spectraCC(model, EPSILON); t = toc;
            case 'fastercc_r0'
                tic; [A, ~, nLP] = fasterCC(model, EPSILON, 0, 0); t = toc;
                cons = A(:); lps = nLP.total;
            otherwise
                error('harvey_fixpoint:arm', 'unknown arm "%s"', arm);
        end
    catch ME
        errMsg = ME.message;
        failed = true;
        fprintf('  FAILED after %.1f s: %s\n', t, errMsg);
    end

    cons = cons(:);
    % err as a string, not a char row: struct2table turns char fields of
    % differing length into a ragged column and errors, which would lose the
    % summary on exactly the iteration that failed.
    rec = struct('iter', it, 'nRxnsIn', nIn, 'nMetsIn', numel(model.mets), ...
                 'nCons', numel(cons), 'nDropped', nIn - numel(cons), ...
                 'lps', lps, 'time', t, 'stable', false, 'err', string(errMsg));

    if ~failed
        fprintf('  %8.1f s   |cons| = %d   dropped = %d   LPs = %s\n', ...
                t, numel(cons), nIn - numel(cons), num2str(lps));
        % Idempotent iff the pass returns everything it was given.
        stable     = isequal(sort(cons), (1:nIn)');
        rec.stable = stable;
    end

    % Persist BEFORE deciding what to do next, so this iteration's record
    % survives whatever the next one does to the job.
    consNames = {};
    if ~isempty(cons) && max(cons) <= nIn, consNames = model.rxns(cons); end
    S = struct('arm', arm, 'iter', it, 'epsilon', EPSILON, ...
               'feasTol', feasTolEff, 'optTol', optTolEff, ...
               'rec', rec, 'consRxns', {consNames}, 'rxnsIn', {model.rxns}, ...
               'lb', model.lb, 'ub', model.ub);
    save(fullfile(outDir, sprintf('iter%02d.mat', it)), 'S', '-v7.3');
    R(end+1) = rec; %#ok<AGROW>
    writetable(struct2table(R), fullfile(outDir, 'summary.csv'));

    if failed, break; end
    if stable
        fprintf('\nSTABLE at iteration %d: the pass returned all %d reactions.\n', it, nIn);
        break
    end
    if isempty(cons)
        fprintf('\nCOLLAPSED at iteration %d: consistent set empty, nothing to iterate on.\n', it);
        break
    end

    % Build the next model from the consistent reactions only.
    drop  = setdiff((1:nIn)', cons);
    model = removeRxns(model, model.rxns(drop));
end

if ~stable && ~failed && it == MAXITER
    fprintf('\nNOT STABLE after %d iterations - still shrinking. Cap reached.\n', MAXITER);
end

%% report ----------------------------------------------------------------
fprintf('\n=== %s: %d iteration(s) ===\n', arm, numel(R));
fprintf('%5s %10s %10s %10s %8s %10s %7s\n', ...
        'iter', 'rxns in', '|cons|', 'dropped', 'LPs', 'time s', 'stable');
for i = 1:numel(R)
    fprintf('%5d %10d %10d %10d %8s %10.1f %7d\n', ...
            R(i).iter, R(i).nRxnsIn, R(i).nCons, R(i).nDropped, ...
            num2str(R(i).lps), R(i).time, R(i).stable);
end
fprintf('\ntotal time %.1f s, total LPs %s\n', ...
        sum([R.time], 'omitnan'), num2str(sum([R.lps], 'omitnan')));
fprintf('results -> %s\n', outDir);
end

% -------------------------------------------------------------------------
function checkCouplingWidth(m)
% Harvey carries coupling constraints (C, d, ctrs, dsense). removeRxns
% checks size(model.C,2) == numel(model.rxns) and errors if they disagree -
% which is how a malformed C surfaces, as "length of rowBool must equal
% size(A,2)" several minutes into a run rather than here. Fail up front.
%
% Note the consistency checkers themselves read only S, lb and ub, so the
% coupling constraints do not participate in the answer; they only have to
% survive removeRxns between iterations.
if isfield(m, 'C') && ~isempty(m.C) && size(m.C, 2) ~= numel(m.rxns)
    error('harvey_fixpoint:badC', ...
        ['model.C is %dx%d but the model has %d reactions. removeRxns ' ...
         'would fail mid-iteration.'], size(m.C,1), size(m.C,2), numel(m.rxns));
end
end
