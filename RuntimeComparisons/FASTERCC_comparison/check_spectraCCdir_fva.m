function R = check_spectraCCdir_fva(cfg)
% check_spectraCCdir_fva  Validate spectraCCdir's directional labels against
% FVA, and against fasterCC revopt=1.
%
%   R = check_spectraCCdir_fva()
%
% FVA is ground truth: it computes each reaction's true min and max flux, so
% it settles directionality exactly. This is the same check the FASTERCC paper
% applied to its own revopt=1 reversibility correction (its Table 3, FVA on
% Recon3D / iHuman1 / iMM1415 / iMM1865), so running it here puts the two
% methods on the same evidential footing.
%
% TWO MODELS, because the two error modes need separating:
%
%   consistent   full consRecon3D, 10,600 reactions, entirely flux consistent
%                (fastcc confirms this). FVA therefore finds ZERO blocked
%                reactions, so every judgement is forward vs reverse vs
%                reversible. This isolates direction discrimination.
%
%   blocked500   the same model with 500 reactions deleted at random, which
%                strands pathways and plants genuinely blocked reactions.
%                This is the paper's own mechanism for creating them, and it
%                isolates blocked detection.
%
% A uniformly random SUBNETWORK was rejected for this job: deleting most of a
% metabolic network leaves mostly rubble, so the test set would be dominated
% by blocked reactions and would barely exercise direction labelling at all.
% Random deletion of a small fraction (blocked500) plants blocked reactions
% without destroying the network.
%
% WHAT COUNTS AS A BUG
%   spectraCCdir can only ever UNDER-report feasibility. Every positive label
%   comes from an actual feasible flux vector, so forward/reverse/reversible
%   are sound by construction. Only 'blocked' is a claim about absence, and
%   spectraCC's loop has no one-by-one fallback, so a reaction can be feasible
%   yet unproven. Hence:
%     claimed feasible, FVA says blocked  -> UNSOUND, a real bug
%     claimed blocked, FVA says feasible  -> a miss, fixable via confirmBlocked
%   The report keeps these apart; they mean very different things.

if nargin < 1 || isempty(cfg), cfg = bench_config('smoke'); end

eps_   = cfg.epsilon;
outDir = fullfile(cfg.root, 'results', 'dir_fva');
if ~exist(outDir,'dir'), mkdir(outDir); end

S = load(fullfile(cfg.baseDir,'consRecon3D.mat'));
base = S.consRecon3D;
changeCobraSolver(cfg.solver,'all',0);

% Case 2: plant blocked reactions by deleting 500 reactions at random.
nRemove = 500;
rng(cfg.baseSeed, 'twister');
drop = sort(randsample(numel(base.rxns), nRemove));
blocked500 = removeRxns(base, base.rxns(drop));

cases = { 'consistent', base; ...
          'blocked500', blocked500 };

R = struct();
for c = 1:size(cases,1)
    name  = cases{c,1};
    model = cases{c,2};
    n     = numel(model.rxns);

    fprintf('\n################ %s : %d reactions, %d metabolites ################\n', ...
            name, n, numel(model.mets));

    %% spectraCCdir ------------------------------------------------------
    fprintf('spectraCCdir ... ');
    t0 = tic;
    [clsSpec, infoSpec] = spectraCCdir(model, eps_, 'stoichiometry', cfg.baseSeed);
    tSpec = toc(t0);
    fprintf('%.1f s, %d LPs\n', tSpec, infoSpec.LPS);
    report('spectraCCdir', clsSpec);

    %% spectraCCdir with one-by-one confirmation -------------------------
    fprintf('spectraCCdir (confirmBlocked) ... ');
    t0 = tic;
    [clsSpecC, infoSpecC] = spectraCCdir(model, eps_, 'stoichiometry', ...
                                cfg.baseSeed, [], struct('confirmBlocked', true));
    tSpecC = toc(t0);
    fprintf('%.1f s, %d LPs\n', tSpecC, infoSpecC.LPS);
    report('spectraCCdir+confirm', clsSpecC);

    %% fasterCC revopt=1 -------------------------------------------------
    fprintf('fasterCC revopt=1 ... ');
    t0 = tic;
    [A, mOut] = fasterCC(model, eps_, 0, 1);
    tFast = toc(t0);
    fprintf('%.1f s\n', tFast);

    % Direction comes from the bounds fasterCC tightened on the surviving set.
    clsFast = repmat(categorical({'blocked'}), n, 1);
    for i = 1:numel(A)
        lo = mOut.lb(i); hi = mOut.ub(i);
        if     lo < -eps_ && hi > eps_, clsFast(A(i)) = categorical({'reversible'});
        elseif hi >  eps_,              clsFast(A(i)) = categorical({'forward'});
        elseif lo < -eps_,              clsFast(A(i)) = categorical({'reverse'});
        end
    end
    report('fasterCC revopt=1', clsFast);

    %% FVA ground truth --------------------------------------------------
    fprintf('FVA (%d LPs) ... ', 2*n);
    t0 = tic;
    [minF, maxF] = fluxVariability(model, 0, 'max');
    tFVA = toc(t0);
    fprintf('%.1f s\n', tFVA);

    clsFVA = repmat(categorical({'blocked'}), n, 1);
    fwd = maxF >  eps_;
    rev = minF < -eps_;
    clsFVA(fwd & ~rev) = categorical({'forward'});
    clsFVA(~fwd & rev) = categorical({'reverse'});
    clsFVA(fwd &  rev) = categorical({'reversible'});
    report('FVA (truth)', clsFVA);

    %% Agreement ---------------------------------------------------------
    entry = struct('n', n, ...
        'times', struct('spectraCCdir',tSpec,'spectraCCdirConfirm',tSpecC, ...
                        'fasterCC_r1',tFast,'FVA',tFVA), ...
        'LPS',   struct('spectraCCdir',infoSpec.LPS,'spectraCCdirConfirm',infoSpecC.LPS));

    methods = {'spectraCCdir', clsSpec; ...
               'spectraCCdir+confirm', clsSpecC; ...
               'fasterCC revopt=1', clsFast};

    for mI = 1:size(methods,1)
        mName = methods{mI,1}; cls = methods{mI,2};
        fprintf('\n--- %s vs FVA ---\n', mName);
        agree = sum(cls == clsFVA);
        fprintf('exact label match                : %d/%d (%.2f%%)\n', agree, n, 100*agree/n);

        feasT = clsFVA ~= 'blocked';
        feasP = cls    ~= 'blocked';
        overclaim  = sum(feasP & ~feasT);
        underclaim = sum(~feasP &  feasT);
        both       = feasT & feasP;
        dirMismatch = sum(cls(both) ~= clsFVA(both));

        fprintf('claimed feasible, FVA blocked    : %d   <- UNSOUND if > 0\n', overclaim);
        fprintf('claimed blocked, FVA feasible    : %d   <- missed, not unsound\n', underclaim);
        fprintf('direction disagreement (shared)  : %d/%d\n', dirMismatch, sum(both));

        key = matlab.lang.makeValidName(mName);
        entry.(key) = struct('exact',agree,'overclaim',overclaim, ...
                             'underclaim',underclaim,'dirMismatch',dirMismatch, ...
                             'confusion',confusion(cls, clsFVA));
    end

    R.(name) = entry;
    save(fullfile(outDir, sprintf('dir_fva_%s.mat', name)), ...
         'clsSpec','clsSpecC','clsFast','clsFVA','minF','maxF','entry','-v7.3');
end

save(fullfile(outDir,'dir_fva_summary.mat'), 'R', '-v7.3');
fprintf('\nSaved results under %s\n', outDir);
end

% ------------------------------------------------------------------ helpers
function report(name, cls)
cats = {'forward','reverse','reversible','blocked'};
line = sprintf('  %-22s', name);
for c = 1:numel(cats)
    line = [line sprintf('%s=%-6d ', cats{c}, sum(cls == cats{c}))]; %#ok<AGROW>
end
disp(line);
end

function M = confusion(a, b)
% rows = method label, cols = FVA label, order forward/reverse/reversible/blocked
cats = {'forward','reverse','reversible','blocked'};
M = zeros(numel(cats));
for i = 1:numel(cats)
    for j = 1:numel(cats)
        M(i,j) = sum(a == cats{i} & b == cats{j});
    end
end
end
