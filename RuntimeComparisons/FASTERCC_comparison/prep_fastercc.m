function [prep, map, info] = prep_fastercc(model)
% prep_fastercc  FASTERCC's structural preprocessing, applied on its own.
%
%   [prep, map, info] = prep_fastercc(model)
%
% Runs exactly the preprocessing block that fasterCC performs before it solves
% anything (fasterCC.m:40-65), by CALLING ITS STOCK FUNCTIONS in the same
% order. Nothing in rFASTCORMICS is copied or modified - this is composition,
% not a fork:
%
%   1. relax bounds so the zero flux vector is feasible (lb>0 -> 0, ub<0 -> 0)
%   2. findDeadEndsFastbox      remove type I and II dead ends
%   3. structureAnalyseFastbox  turn structurally one-way reversibles irreversible
%   4. fixIrr_rFASTCORMICS      orient irreversibles forward (lb>=0, ub>0)
%
% WHY THIS EXISTS
% FASTERCC's speedup has two separable sources: (a) this structural pass,
% which shrinks the problem before any LP is solved, and (b) the LP strategy
% that follows. The paper never applies (a) to any baseline, so it cannot show
% how much of its reported speedup is actually (b). Feeding a preprocessed
% model to stock spectraCC separates them.
%
% OUTPUTS
%   prep   preprocessed model - FEWER reactions than the input
%   map    numel(prep.rxns) x 1, index of each surviving reaction in the
%          ORIGINAL model. Use it to translate any result back:
%              consistentInOriginal = map(consistentInPrep);
%          Reactions dropped by step 2 are blocked by construction and are
%          simply absent from map.
%   info   .nRemoved  reactions dropped by the structural pass
%          .nDerev    reversible -> irreversible conversions
%          .time      seconds spent preprocessing (charge this to the caller)

t0 = tic;

rxns0 = model.rxns;
nRev0 = nnz(model.lb < 0 & model.ub > 0);

% 1. Relax bounds exactly as fasterCC does before its structural pass.
model.lb(model.lb > 0) = 0;
model.ub(model.ub < 0) = 0;

% 2-4. Stock FASTERCC functions, in fasterCC's own order.
model = findDeadEndsFastbox(model);
model = structureAnalyseFastbox(model);
model = fixIrr_rFASTCORMICS(model);

prep = model;

% Map survivors back to original indices by NAME, so this stays correct
% regardless of how the stock functions reindex internally.
[tf, loc] = ismember(prep.rxns, rxns0);
assert(all(tf), 'prep_fastercc:unmapped', ...
    'Preprocessing produced %d reaction(s) absent from the input model.', sum(~tf));
map = loc(:);

info.nRemoved = numel(rxns0) - numel(prep.rxns);
info.nDerev   = nRev0 - nnz(prep.lb < 0 & prep.ub > 0);
info.time     = toc(t0);
end
