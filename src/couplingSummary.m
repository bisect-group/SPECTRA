function text = couplingSummary(community)
% Reports what each unit of a community is coupled to, and how hard
%
% USAGE:
%
%    text = couplingSummary(community)
%    couplingSummary(community)
%
% INPUTS:
%    community: Structure returned by buildCommunityModel.
%
% OUTPUTS:
%    text:      One line per unit, as char, giving its anchor, how the
%               anchor was decided, the most flux the anchor can carry and
%               the cap that imposes on the rest of the unit. Printed
%               instead when no output is asked for.
%
% NOTE:
%    The cap is ``couplingFactor * capacity + couplingThreshold``: the most
%    flux any reaction of the unit may carry when its anchor runs as fast
%    as it can. A cap far below the unit's own bounds means the coupling is
%    throttling it rather than merely tying it down, and the factor wants
%    raising. A capacity of zero means the unit cannot grow at all in this
%    community, which coupling then pins every one of its reactions to
%    zero: look at the medium before looking at the factor.
%
%    The capacity is only measured when buildCommunityModel was asked to
%    check, which is its default. A dash stands where it was not.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: buildCommunityModel, communityDecompose

if ~isstruct(community) || ~isfield(community, 'units')
    error('couplingSummary:badCommunity', ...
        'community must be the structure buildCommunityModel returns.');
end

if isempty(community.couplingFactor)
    lines = {'not coupled'};
else
    width = max([4, cellfun(@numel, community.organisms(:)')]);
    anchorWidth = max([6, cellfun(@numel, community.biomassReactions(:)')]);
    header = sprintf('%-*s  %-*s  %-10s  %10s  %10s', ...
        width, 'unit', anchorWidth, 'anchor', 'how', 'capacity', 'cap');
    lines = cell(numel(community.units) + 1, 1);
    lines{1} = header;
    for k = 1:numel(community.units)
        capacity = community.units(k).anchorCapacity;
        if isempty(capacity) || isnan(capacity)
            shown = '-';
            cap = '-';
        else
            shown = sprintf('%.4g', capacity);
            cap = sprintf('%.4g', ...
                community.couplingFactor * capacity + community.couplingThreshold);
        end
        lines{k+1} = sprintf('%-*s  %-*s  %-10s  %10s  %10s', ...
            width, community.units(k).name, ...
            anchorWidth, community.units(k).anchor, ...
            community.units(k).anchorSource, shown, cap);
    end
end

joined = strjoin(lines(:)', newline);
if nargout == 0
    fprintf('%s\n', joined);
    clear text
else
    text = joined;
end
end
