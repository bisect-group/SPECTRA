function params = milpTolerances(solverTol)
% Returns the tolerances a mixed-integer formulation needs, as solver arguments
%
% USAGE:
%
%    params = milpTolerances(solverTol)
%
% OPTIONAL INPUTS:
%    solverTol: The tolerance to ask for, applied to all three of the
%               feasibility, optimality and integrality tolerances when a
%               scalar is given, or 'default' to return {} and leave the
%               solver's / COBRA session's tolerances in place.
%               (Default when omitted/empty: feasTol = 1e-8, optTol = 1e-8,
%               intTol = 1e-9)
%
% OUTPUTS:
%    params:    Cell array of name/value pairs to pass to solveCobraMILP,
%               as {'feasTol',t,'optTol',t,'intTol',t}.
%
% NOTE:
%    These are passed to the solve rather than set with
%    changeCobraSolverParams, because that is a global COBRA setting: a
%    formulation that changed it would leave every later call in the
%    session running at a tolerance nobody asked for, including plain
%    optimizeCbModel. solveCobraMILP takes the same names as arguments, so
%    the problem can carry its own tolerances and nothing outside it moves.
%
%    The integrality tolerance is the one that matters here, and the one
%    most easily overlooked. A binary is tied to its reaction by
%    lb_i*z_i <= v_i <= ub_i*z_i, so a binary sitting at the tolerance and
%    still counted as zero permits a flux of tolerance*ub_i. At a solver's
%    usual 1e-5 and the usual bound of 1000 that is a flux of 0.01, which
%    is far above inclusionCutoff, so reactions the solver switched *off*
%    read as carrying flux and are kept. The feasibility tolerance is the
%    complementary failure: too loose and a core reaction is present but
%    cannot carry the flux it was required to.
%
%    1e-9 is a floor, not a preference. Gurobi refuses any of the three
%    below it and keeps whatever it already had, which is looser than
%    either, so a smaller value is raised rather than passed through and
%    quietly dropped.
%
%    Loosening them is occasionally the right move: a large MILP that will
%    not converge at 1e-9 may solve readily at 1e-7, at the cost of having
%    to trust inclusionCutoff rather more.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: minNet, tradeOff, solveCobraMILP

FLOOR = 1e-9;

if nargin < 1 || isempty(solverTol)
    params = {'feasTol', 1e-8, 'optTol', 1e-8, 'intTol', FLOOR};
    return
end
if (ischar(solverTol) || (isstring(solverTol) && isscalar(solverTol))) && ...
        strcmpi(solverTol, 'default')
    params = {};
    return
end
if ~isnumeric(solverTol) || ~isscalar(solverTol) || ~isfinite(solverTol) || ...
        solverTol <= 0
    error('milpTolerances:badSolverTol', ...
        'solverTol must be a positive finite scalar or ''default''.');
end
if solverTol < FLOOR
    warning('milpTolerances:belowFloor', ...
        ['The solver tolerance cannot go below %.0e, so %.0e has been ' ...
         'raised to it. Asking for less would have the solver keep its ' ...
         'own tolerance instead, which is looser than either.'], ...
        FLOOR, solverTol);
    solverTol = FLOOR;
end

params = {'feasTol', solverTol, 'optTol', solverTol, 'intTol', solverTol};
end
