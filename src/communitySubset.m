function reduced = communitySubset(community, keep)
% Restricts a community to some of its units
%
% USAGE:
%
%    reduced = communitySubset(community, keep)
%
% INPUTS:
%    community: Structure returned by buildCommunityModel.
%    keep:      Cell array of the unit names to retain, or a logical or
%               numeric index into community.organisms.
%
% OUTPUTS:
%    reduced:   A community holding only those units. The pools and their
%               exchanges are kept, so the survivors still trade with one
%               another and still face the same environment: dropping a
%               unit is meant to answer what the rest can do without it,
%               and rebuilding the pool around the survivors would change
%               the question.
%
% EXAMPLE:
%
%    ::
%
%        community = buildCommunityModel({a, b, c}, {'a','b','c'}, ...);
%        pair = communitySubset(community, {'a','c'});
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: buildCommunityModel, communityWithModel, minimalMicrobiome

if ~isstruct(community) || ~isfield(community, 'units')
    error('communitySubset:badCommunity', ...
        'community must be the structure buildCommunityModel returns.');
end

organisms = community.organisms(:);
if islogical(keep) || isnumeric(keep)
    wanted = organisms(keep);
elseif ischar(keep)
    wanted = {keep};
elseif iscell(keep)
    wanted = keep(:);
else
    error('communitySubset:badKeep', ...
        'keep must be unit names, or an index into community.organisms.');
end
wanted = unique(wanted, 'stable');

unknown = setdiff(wanted, organisms);
if ~isempty(unknown)
    error('communitySubset:unknownUnit', ...
        'Not a unit of this community: %s.', strjoin(unknown(:)', ', '));
end
if isempty(wanted)
    error('communitySubset:emptySubset', ...
        'Keeping no units at all would leave nothing to subset.');
end

alive = ismember(organisms, wanted);
leaving = vertcat(community.units(~alive).rxns);
if isempty(leaving)
    model = community.model;
else
    model = removeRxns(community.model, leaving);
end

reduced = community;
reduced.organisms = organisms(alive);
reduced.biomassReactions = community.biomassReactions(alive);
reduced.units = community.units(alive);
reduced.draftReactions = community.draftReactions(alive);
reduced.databaseReactions = community.databaseReactions(alive);
for p = 1:numel(reduced.pools)
    reduced.pools(p).units = ...
        reduced.pools(p).units(ismember(reduced.pools(p).units, wanted));
end
reduced = communityWithModel(reduced, model);

% The objective follows the survivors: a dropped unit's anchor is gone, and
% leaving a stale c would make the reduced model optimise nothing at all.
reduced.model.c(:) = 0;
at = findRxnIDs(reduced.model, reduced.biomassReactions);
reduced.model.c(at(at > 0)) = 1;
end
