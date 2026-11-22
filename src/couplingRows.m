function [A, b, csense] = couplingRows(model, nAux, useCoupling)
% Returns a model's coupling constraints as rows of a solver problem
%
% USAGE:
%
%    [A, b, csense] = couplingRows(model, nAux, useCoupling)
%
% INPUTS:
%    model:   COBRA model structure. Coupling constraints live in the
%             fields C (one row per constraint, one column per reaction),
%             d (the right-hand side) and dsense (its direction, as 'L',
%             'G' or 'E' per row).
%    nAux:    Number of auxiliary variables the problem carries after the
%             fluxes, so that the rows can be padded to its width.
%             (Default: 0)
%    useCoupling: Whether to return the constraints at all. False gives no
%             rows, so the problem is built from S, lb and ub alone.
%             (Default: false)
%
% OUTPUTS:
%    A:       The coupling rows, [C, 0], as a sparse matrix of nCtrs rows
%             and numel(model.rxns)+nAux columns. Empty with the right
%             width for a model that has no coupling constraints.
%    b:       Their right-hand side.
%    csense:  Their direction, one character per row.
%
% NOTE:
%    Every formulation in SPECTRA builds its own problem out of model.S,
%    so a constraint that is not in S is not in the problem unless it is
%    added here.
%
%    Off by default, because flux consistency and the extractions built on
%    it are conventionally defined over the stoichiometry alone - which is
%    what fastcc, FASTCORE, FASTERCC and swiftcore all answer, and what a
%    comparison against them has to answer too. A coupling constraint is a
%    modelling device rather than stoichiometry, so a reaction blocked only
%    by one is not blocked in the network; it is blocked in that
%    parameterisation.
%
%    Switch it on when the parameterisation is the point. A community model
%    ties each unit's reactions to its own anchor through C, and an
%    extraction that ignores C is free to run a unit's metabolism while the
%    unit is not growing, which is the whole thing the coupling forbids:
%    gapfillCommunity and minimalMicrobiome therefore ask for it. Switch it
%    on, too, whenever the answer has to agree with optimizeCbModel on the
%    same model, since that honours C always. On a toy community whose
%    consumer cannot grow, FBA puts a transport's range at [-0.01, 0] while
%    a check ignoring C reports it able to carry 0.1.
%
%    A model without the fields, which is every model that has never been
%    through coupleRxnList2Rxn, gives no rows either way.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: coupleRxnList2Rxn, minNet, tradeOff, growthOptim

if nargin < 2 || isempty(nAux)
    nAux = 0;
end
if nargin < 3 || isempty(useCoupling)
    useCoupling = false;
end
n = numel(model.rxns);

if ~useCoupling || ~isfield(model, 'C') || isempty(model.C)
    A = sparse(0, n + nAux);
    b = zeros(0, 1);
    csense = repmat('L', 0, 1);
    return
end

C = model.C;
if size(C, 2) ~= n
    error('couplingRows:badC', ...
        ['model.C has %d columns but the model has %d reactions, so the ' ...
         'coupling constraints cannot be read.'], size(C, 2), n);
end
nCtrs = size(C, 1);

if ~isfield(model, 'd') || numel(model.d) ~= nCtrs
    error('couplingRows:badD', ...
        ['model.d must hold one right-hand side per coupling constraint, ' ...
         'that is %d of them.'], nCtrs);
end
if ~isfield(model, 'dsense') || numel(model.dsense) ~= nCtrs
    error('couplingRows:badDsense', ...
        ['model.dsense must hold one direction per coupling constraint, ' ...
         'that is %d of them.'], nCtrs);
end

A = [sparse(C), sparse(nCtrs, nAux)];
b = reshape(model.d, [], 1);
csense = reshape(model.dsense, [], 1);
end
