function [groupOf,nGroups,ungoverned] = resolveIndicatorGroups(indicatorGroups,dir0,n)
% Works out which binary governs which reaction in a MILP formulation
%
% USAGE:
%
%    [groupOf,nGroups,ungoverned] = resolveIndicatorGroups(indicatorGroups,dir0,n)
%
% INPUTS:
%    indicatorGroups: Cell array of index vectors, one per group, giving the
%                     reactions that share a single binary. Empty for the
%                     usual formulation, one binary per free reaction.
%    dir0:            Logical over the reactions, true where a reaction is
%                     free to be dropped, i.e. where direction is 0. A
%                     reaction given a direction is already forced in and
%                     cannot also be selected over.
%    n:               Number of reactions in the model.
%
% OUTPUTS:
%    groupOf:         The group index of every reaction, 0 where no binary
%                     governs it.
%    nGroups:         How many binaries there are.
%    ungoverned:      Logical over the reactions, true for the free ones no
%                     binary governs. Nothing is deciding whether to keep
%                     them, so nothing should discard them.
%
% NOTE:
%    Sharing a binary between several reactions is what expresses a
%    requirement that holds in several conditions at once: replicate the
%    network once per condition, group each reaction's copies together, and
%    the objective counts reactions while the constraints are satisfied
%    separately in every condition. That is what minimalReactome does.
%
%    Naming one reaction per group is the other use, and it is what makes a
%    selection problem over a handful of reactions tractable on a large
%    model: one binary per organism in a community, say, rather than one per
%    reaction. Every free reaction left out of the groups is kept.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: minNet, tradeOff, minimalReactome

groupOf = zeros(n,1);

if isempty(indicatorGroups)
    free = find(dir0);
    groupOf(free) = 1:numel(free);
    nGroups = numel(free);
    ungoverned = false(n,1);
    return
end

if ~iscell(indicatorGroups)
    error('resolveIndicatorGroups:badIndicatorGroups', ...
        'indicatorGroups must be a cell array of reaction index vectors.');
end
nGroups = numel(indicatorGroups);
for g = 1:nGroups
    members = indicatorGroups{g}(:);
    if isempty(members)
        error('resolveIndicatorGroups:emptyIndicatorGroup', ...
            ['Indicator group %d is empty, so its binary would govern ' ...
             'nothing while still being counted.'], g);
    end
    if any(members < 1 | members > n | members ~= round(members))
        error('resolveIndicatorGroups:badIndicatorGroup', ...
            'Indicator group %d names something that is not a reaction index.', g);
    end
    if any(~dir0(members))
        error('resolveIndicatorGroups:directedIndicatorGroup', ...
            ['Indicator group %d names a reaction that has been given a ' ...
             'direction. Such a reaction is already forced in and cannot ' ...
             'also be selected over.'], g);
    end
    clash = members(groupOf(members) > 0);
    if ~isempty(clash)
        error('resolveIndicatorGroups:overlappingIndicatorGroups', ...
            ['Reaction %d is in both group %d and group %d, so which binary ' ...
             'governs it is ambiguous. Every reaction belongs to at most ' ...
             'one.'], clash(1), groupOf(clash(1)), g);
    end
    groupOf(members) = g;
end
ungoverned = dir0 & groupOf == 0;
end
