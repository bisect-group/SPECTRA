function organism = communityOrganismOf(community, identifiers)
% Says which unit of a community a reaction or metabolite belongs to
%
% USAGE:
%
%    organism = communityOrganismOf(community, identifiers)
%
% INPUTS:
%    community:   Structure returned by buildCommunityModel.
%    identifiers: A reaction or metabolite identifier as char, or a cell
%                 array of them.
%
% OUTPUTS:
%    organism:    The unit's name, or '' for something belonging to no unit:
%                 a pool exchange, a pool metabolite, or a shared
%                 metabolite. A cell array of the same shape when a cell
%                 array was given.
%
% NOTE:
%    The tag is read off the identifier rather than looked up in the unit
%    lists, so this answers for a reaction the community model no longer
%    holds as readily as for one it does. A name is only accepted as a unit
%    if the community has a unit by that name, so a reaction that happens
%    to end in '__something' is not mistaken for one.
%
% EXAMPLE:
%
%    ::
%
%        communityOrganismOf(community, 'biomass__org1')   % 'org1'
%        communityOrganismOf(community, 'EX_glc[u]')       % ''
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: buildCommunityModel, communitySubset, communityDecompose

if ischar(identifiers) || (isstring(identifiers) && isscalar(identifiers))
    organism = lookup(char(identifiers), community.organisms);
    return
end
if ~iscell(identifiers)
    error('communityOrganismOf:badIdentifier', ...
        'identifiers must be a char or a cell array of char.');
end
organism = cell(size(identifiers));
for i = 1:numel(identifiers)
    organism{i} = lookup(char(identifiers{i}), community.organisms);
end
end


function name = lookup(identifier, organisms)
% Returns the unit named by an identifier's tag, or '' where there is none.
name = '';

% A metabolite wears its tag inside the compartment, a reaction at the end.
token = regexp(identifier, '^.*\[(\w+)\]$', 'tokens', 'once');
if ~isempty(token)
    identifier = token{1};
end

at = strfind(identifier, '__');
if isempty(at)
    return
end
candidate = identifier(at(end)+2:end);
if ismember(candidate, organisms)
    name = candidate;
end
end
