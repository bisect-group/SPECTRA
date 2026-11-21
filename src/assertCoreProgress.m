function stalls = assertCoreProgress(remaining, before, stalls, tol, limit)
% Stops a direction phase that is feasible but getting nowhere
%
% USAGE:
%
%    stalls = assertCoreProgress(remaining, before, stalls, tol, limit)
%
% INPUTS:
%    remaining: Logical over the reactions, true for the core reactions the
%               direction phase has still to push above the threshold.
%    before:    How many were left at the start of this round.
%    stalls:    How many rounds in a row have cleared nothing so far.
%
% OPTIONAL INPUTS:
%    tol:       The flux threshold in force, for the message. (Default: NaN)
%    limit:     How many fruitless rounds to allow before giving up.
%               (Default: 3)
%
% OUTPUTS:
%    stalls:    The updated count, zero where this round cleared something.
%
% NOTE:
%    The direction phase loops until every core reaction carries at least a
%    tenth of tol, which assumes each round clears at least one. An LP that
%    is feasible but cannot lift a reaction that far breaks the assumption,
%    and the loop then runs forever on an LP that keeps succeeding. A
%    coupling constraint is the case that brings it about in practice: a
%    unit whose anchor cannot carry flux has every one of its reactions
%    capped at couplingThreshold, so a core reaction of that unit is
%    feasible at 0.01 and unreachable at a tol above it.
%
%    An all-NaN flux vector, which is what an infeasible LP gives, is the
%    other failure and is caught where it arises rather than here.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: spectraME, spectraCCME, minimalReactome

if nargin < 4 || isempty(tol)
    tol = NaN;
end
if nargin < 5 || isempty(limit)
    limit = 3;
end

if sum(remaining) < before
    stalls = 0;
    return
end

stalls = stalls + 1;
if stalls < limit
    return
end

at = find(remaining, 3);
error('spectra:stalledCore', ...
    ['The direction phase has gone %d rounds without getting any of the ' ...
     '%d remaining core reactions above tol=%g. Their flux is feasible ' ...
     'but below the threshold, which no further round can change, so the ' ...
     'loop would run forever. Either tol is higher than these reactions ' ...
     'can carry, or something caps them: a coupling constraint pins every ' ...
     'reaction of a unit whose anchor cannot grow to couplingThreshold. ' ...
     'The first few are reactions %s. spectraCCME drops core reactions ' ...
     'that cannot carry flux and reports them as BlockedCoreRxns.'], ...
    stalls, sum(remaining), tol, mat2str(at(:)'));
end
