function [A, b, csense] = couplingRows(model, nAux)
% Returns a model's coupling constraints as rows of a solver problem
%
% USAGE:
%
%    [A, b, csense] = couplingRows(model, nAux)
%
% INPUTS:
%    model:   COBRA model structure. Coupling constraints live in the
%             fields C (one row per constraint, one column per reaction),
%             d (the right-hand side) and dsense (its direction, as 'L',
%             'G' or 'E' per row).
%    nAux:    Number of auxiliary variables the problem carries after the
%             fluxes, so that the rows can be padded to its width.
%             (Default: 0)
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
%    added here. Coupling constraints are the case that matters: a
%    community model ties each unit's reactions to its own anchor through
%    C, and an extraction that ignored C would be free to run a unit's
%    metabolism while the unit is not growing. optimizeCbModel honours
%    them, so a model whose FBA and whose extraction disagreed would be
%    the more confusing outcome.
%
%    A model without the fields, which is every model that has never been
%    through coupleRxnList2Rxn, gives no rows and is therefore unaffected.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: coupleRxnList2Rxn, minNet, tradeOff, growthOptim

if nargin < 2 || isempty(nAux)
    nAux = 0;
end
n = numel(model.rxns);

if ~isfield(model, 'C') || isempty(model.C)
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
