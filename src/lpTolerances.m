function params = lpTolerances(solverTol)
% Returns the tolerances a linear formulation needs, as solver arguments
%
% USAGE:
%
%    params = lpTolerances(solverTol)
%
% OPTIONAL INPUTS:
%    solverTol: The tolerance to ask for, applied to both the feasibility
%               and optimality tolerances, or 'default' to return {} and
%               leave the solver's / COBRA session's tolerances in place.
%               (Default: 1e-9, the tightest the solvers accept)
%
% OUTPUTS:
%    params:    Cell array of name/value pairs to pass to solveCobraLP,
%               as {'feasTol',t,'optTol',t}, or {} when solverTol is 'default'.
%
% NOTE:
%    Delegates validation and the 1e-9 solver floor to milpTolerances and
%    drops the integrality pair, leaving global CBT_LP_PARAMS untouched.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: milpTolerances, forward, forwardcc, reverse, minNet, optimBiomass

if nargin < 1
    solverTol = [];
end
milpParams = milpTolerances(solverTol);
if isempty(milpParams)
    params = {};
else
    params = milpParams(1:4);
end
end
