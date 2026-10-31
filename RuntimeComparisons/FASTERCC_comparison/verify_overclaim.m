function verify_overclaim(cfg)
% verify_overclaim  Independently re-test the reactions SPECTRA overclaimed.
%
% The disagreement FVA on recon_rm1020_rep01 reported that spectraCC calls
% HMR_0015, HMR_0016, EX_M02560[e] and EX_M02561[e] consistent while FVA puts
% both their min and max flux at exactly 0. That is the unsound direction, so
% before it is believed it has to survive a check that shares as little
% machinery as possible with the run that produced it.
%
% Three independent routes to the same question:
%   1. call spectraCC directly and test NAME membership in its output, which
%      bypasses the benchmark's index bookkeeping entirely - if the harness
%      were mismapping indices, this is where it would show
%   2. maximise and minimise each reaction with optimizeCbModel, a different
%      COBRA entry point than fluxVariability
%   3. report the stoichiometric degree of each reaction's metabolites, since
%      a reaction whose metabolites appear nowhere else cannot carry steady
%      state flux and should be blocked on inspection alone

if nargin < 1 || isempty(cfg), cfg = bench_config('smoke'); end
changeCobraSolver(cfg.solver, 'all', 0);
tol = cfg.epsilon;

resDir = fullfile(cfg.root, 'results', 'variants_recon_spectra_first_nofva_rep5_a3');
L = load(fullfile(resDir, 'recon_rm1020_rep01.mat'));
r = L.res;

S     = load(fullfile(cfg.baseDir, 'consRecon3D.mat'));
model = removeRxns(S.consRecon3D, r.removedRxns);
assert(isequal(model.rxns, r.rxns), 'verify_overclaim:rebuild', 'Model mismatch.');

names = r.disagreeRxns;

%% 1. spectraCC, direct ---------------------------------------------------
fprintf('\n--- spectraCC called directly, membership tested by NAME ---\n');
[cons, lps] = spectraCC(model, tol, 'stoichiometry');
consNames = model.rxns(cons);
fprintf('spectraCC: %d consistent of %d, %d LPs\n', numel(cons), numel(model.rxns), lps);
for j = 1:numel(names)
    fprintf('  %-18s in spectraCC output: %d\n', names{j}, ...
            ismember(names{j}, consNames));
end

%% 2. optimizeCbModel, per reaction --------------------------------------
fprintf('\n--- optimizeCbModel, each reaction as the objective ---\n');
fprintf('%-18s %12s %12s %12s\n', 'reaction', 'max', 'min', 'verdict');
for j = 1:numel(names)
    i = find(strcmp(model.rxns, names{j}), 1);
    m = changeObjective(model, names{j});
    sMax = optimizeCbModel(m, 'max');
    sMin = optimizeCbModel(m, 'min');
    vMax = NaN; vMin = NaN;
    if sMax.stat == 1, vMax = sMax.x(i); end
    if sMin.stat == 1, vMin = sMin.x(i); end
    if max(abs([vMax vMin])) >= tol, v = 'consistent'; else, v = 'BLOCKED'; end
    fprintf('%-18s %12.6g %12.6g %12s\n', names{j}, vMax, vMin, v);
end

%% 3. stoichiometric neighbourhood ---------------------------------------
% A reaction can only carry steady state flux if every metabolite it touches
% can be balanced by some other reaction. deg==1 means this reaction is the
% only one touching that metabolite, which forces its flux to zero.
fprintf('\n--- stoichiometric degree of each touched metabolite ---\n');
% full(): S is sparse, so the row sums come back sparse too, and fprintf
% refuses sparse arguments.
deg = full(sum(model.S ~= 0, 2));
for j = 1:numel(names)
    i    = find(strcmp(model.rxns, names{j}), 1);
    mets = find(model.S(:,i) ~= 0);
    fprintf('%-18s bounds [%g %g]  mets:', names{j}, model.lb(i), model.ub(i));
    for q = mets'
        fprintf('  %s(deg %d)', model.mets{q}, deg(q));
    end
    if any(deg(mets) == 1)
        fprintf('   <- dead end, must be blocked');
    end
    fprintf('\n');
end
end
