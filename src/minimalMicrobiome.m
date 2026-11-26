function result = minimalMicrobiome(community, varargin)
% Finds the smallest subset of a community that still does the job
%
% USAGE:
%
%    result = minimalMicrobiome(community, varargin)
%
% INPUTS:
%    community:        Structure returned by buildCommunityModel. Its
%                      biomass coupling is what makes an absent organism
%                      truly absent, so a community built with 'couple'
%                      false will keep organisms that contribute nothing
%                      but phantom flux.
%
% OPTIONAL INPUTS (name/value pairs):
%    products:         The community exchange reactions whose production
%                      must be kept, as a cell array of reaction IDs or a
%                      struct array with fields rxn and weight, for when
%                      one product matters more than another. With none
%                      given only the growth requirement applies.
%                      (Default: none)
%    growthFraction:   Fraction of the full community's growth the minimal
%                      one must keep. (Default: 0.8)
%    productFraction:  Fraction of its production the minimal one must
%                      keep. (Default: 0.8)
%    growthOptimumFraction: While measuring how much the full community can
%                      produce, hold each organism to this fraction of its
%                      own growth. Sub-optimal growth is the more realistic
%                      setting and gives higher production. (Default: 0.99)
%    required:         Cell array of organisms that must be kept whatever
%                      the objective says. (Default: none)
%    tol:              Flux threshold handed to the extraction, and the
%                      flux a required organism's anchor must carry.
%                      (Default: 1e-4)
%    solveTime:        Maximum runtime for the mixed-integer solve, in
%                      seconds. (Default: 300)
%    seed:             Seed for the extraction's random coefficients. There
%                      is nothing random to seed here, since those
%                      coefficients belong to the direction phase and with
%                      no core reactions that phase has nothing to settle;
%                      it is passed through for the sake of consistency.
%                      (Default: none)
%
% OUTPUTS:
%    result:           Structure with fields:
%
%                        * community - the minimal community, with the
%                          absent organisms removed
%                        * present, absent - the organisms kept and dropped
%                        * membership - struct array of name and present
%                        * growth, production - what the minimal community
%                          achieves
%                        * requiredGrowth, requiredProduction - what it had
%                          to reach
%                        * fullGrowth, fullProduction - what the whole
%                          community achieved
%                        * summary - a readable multi-line summary, as char
%
% EXAMPLE:
%
%    ::
%
%        community = buildCommunityModel(models, names, 'poolMedium', medium);
%        result = minimalMicrobiome(community, ...
%            'products', {'EX_but[u]'}, 'productFraction', 0.9);
%        disp(result.summary)
%
% NOTE:
%    This is a selection problem over organisms, and it maps onto the
%    extraction machinery directly. Give each organism's anchor a binary
%    through indicatorGroups and weight it 1, weight everything else 0, and
%    minimising the weighted count becomes minimising the number of
%    organisms present. The biomass coupling does the rest: an organism
%    whose binary is off cannot grow, and so carries no flux anywhere.
%
%    The requirements are measured before they are imposed. The full
%    community is solved for its growth, then for how much it can produce
%    while each organism keeps growthOptimumFraction of its own growth, and
%    the minimal community is held to a fraction of each. Asking for a
%    fraction of something unmeasured would make the result depend on the
%    units of the biomass reaction.
%
%    Membership is read off the binaries rather than off the fluxes of the
%    result, because the two are not the same thing: an organism whose
%    binary is off can still carry a trace of flux below the threshold, and
%    one whose binary is on need not carry any.
%
%    Follows Raghu, A. K., Palanikumar, I. and Raman, K. (2024), Designing
%    function-specific minimal microbiomes from large microbial
%    communities, npj Systems Biology and Applications 10, 46.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: buildCommunityModel, communitySubset, gapfillCommunity

if ~isstruct(community) || ~isfield(community, 'units')
    error('minimalMicrobiome:badCommunity', ...
        'community must be the structure buildCommunityModel returns.');
end

parser = inputParser();
parser.addParameter('products', {}, @(x) iscell(x) || ischar(x) || isstruct(x));
parser.addParameter('growthFraction', 0.8, @isnumeric);
parser.addParameter('productFraction', 0.8, @isnumeric);
parser.addParameter('growthOptimumFraction', 0.99, @isnumeric);
parser.addParameter('required', {}, @(x) iscell(x) || ischar(x));
parser.addParameter('tol', 1e-4, @isnumeric);
parser.addParameter('solveTime', 300, @isnumeric);
parser.addParameter('seed', [], @isnumeric);
parser.parse(varargin{:});
opt = parser.Results;

if isempty(community.couplingFactor)
    warning('minimalMicrobiome:notCoupled', ...
        ['This community was built without coupling, so an organism can ' ...
         'run its metabolism while not growing. Switching its binary off ' ...
         'then stops nothing, and the result will claim organisms are ' ...
         'absent that are still feeding the rest. Build it with ''couple'' ' ...
         'true.']);
end

tol = opt.tol;
organisms = community.organisms(:);
anchors = community.biomassReactions(:);
model = community.model;

required = opt.required;
if ischar(required)
    required = {required};
end
required = unique(required(:), 'stable');
unknown = setdiff(required, organisms);
if ~isempty(unknown)
    error('minimalMicrobiome:unknownOrganism', ...
        'Not an organism of this community: %s.', strjoin(unknown(:)', ', '));
end

[products, productWeights] = resolveProducts(model, opt.products);

%% 1. what the whole community manages
full = optimizeCbModel(changeObjective(model, anchors));
if full.stat ~= 1 || full.f <= 0
    error('minimalMicrobiome:noGrowth', ...
        ['The full community cannot grow (status %d, growth %g), so no ' ...
         'subset of it can either. Check the medium before the ' ...
         'membership.'], full.stat, full.f);
end
fullGrowth = full.f;
perOrganism = full.x(findRxnIDs(model, anchors));

%% 2. how much it can produce while still growing
fullProduction = 0;
if ~isempty(products)
    held = model;
    at = findRxnIDs(held, anchors);
    held.lb(at) = max(held.lb(at), opt.growthOptimumFraction * perOrganism);
    held = changeObjective(held, products, productWeights);
    made = optimizeCbModel(held, 'max');
    if made.stat == 1
        fullProduction = made.f;
    else
        warning('minimalMicrobiome:noProduction', ...
            ['The full community cannot produce the requested metabolites ' ...
             'while every organism keeps %g of its growth (status %d). The ' ...
             'production requirement is therefore zero.'], ...
            opt.growthOptimumFraction, made.stat);
    end
end

requiredGrowth = opt.growthFraction * fullGrowth;
requiredProduction = opt.productFraction * fullProduction;
fprintf(['minimalMicrobiome: the full community grows at %.4g and produces ' ...
         '%.4g; requiring %.4g and %.4g.\n'], fullGrowth, fullProduction, ...
        requiredGrowth, requiredProduction);

%% 3. the smallest membership that still meets both
constrained = addCOBRAConstraints(model, anchors, requiredGrowth, ...
    'c', ones(1, numel(anchors)), 'dsense', 'G', ...
    'ConstraintID', 'spectraCommunityGrowth');
if ~isempty(products) && requiredProduction > 0
    constrained = addCOBRAConstraints(constrained, products, requiredProduction, ...
        'c', productWeights(:)', 'dsense', 'G', ...
        'ConstraintID', 'spectraCommunityProduction');
end

isRequired = ismember(organisms, required);
at = findRxnIDs(constrained, anchors);
constrained.lb(at(isRequired)) = max(constrained.lb(at(isRequired)), tol);

weights = zeros(numel(constrained.rxns), 1);
weights(at(~isRequired)) = 1;
candidates = anchors(~isRequired);

if isempty(candidates)
    selectedRxns = anchors;
else
    % With coupling, and it is not optional here: an organism whose binary
    % is off must be unable to carry flux anywhere, and the coupling is the
    % only thing that says so. The growth and production requirements this
    % call has to respect are coupling constraints too.
    [~, ~, selection] = spectraME(constrained, [], tol, 'stoichiometry', ...
        weights, 1, [], 'minNetMILP', opt.solveTime, false, ...
        {}, opt.seed, candidates, [], [], true);
    if isempty(selection)
        error('minimalMicrobiome:noSelection', ...
            ['The mixed-integer solve returned no membership. It was most ' ...
             'likely infeasible: the growth and production required of the ' ...
             'minimal community cannot both be met.']);
    end
    selectedRxns = selection.selectedRxns;
end

present = organisms(isRequired | ismember(anchors, selectedRxns));
absent = setdiff(organisms, present, 'stable');

%% 4. what the minimal community actually achieves
minimal = communitySubset(community, present);
kept = minimal.biomassReactions(:);
grown = optimizeCbModel(changeObjective(minimal.model, kept));
growth = 0;
if grown.stat == 1
    growth = grown.f;
end

production = 0;
keptProducts = products(ismember(products, minimal.model.rxns));
if ~isempty(keptProducts)
    % Measured the way the requirement was: the most it can produce while
    % still growing. Reading production off the growth optimum instead
    % would report whatever the solver happened to do with the slack,
    % which is usually nothing.
    checking = addCOBRAConstraints(minimal.model, kept, requiredGrowth, ...
        'c', ones(1, numel(kept)), 'dsense', 'G', ...
        'ConstraintID', 'spectraCheckGrowth');
    keptWeights = productWeights(ismember(products, minimal.model.rxns));
    made = optimizeCbModel(changeObjective(checking, keptProducts, keptWeights), 'max');
    if made.stat == 1
        production = made.f;
    end
end

result.community = minimal;
result.present = present;
result.absent = absent;
result.membership = struct('name', organisms, ...
    'present', num2cell(ismember(organisms, present)));
result.growth = growth;
result.production = production;
result.requiredGrowth = requiredGrowth;
result.requiredProduction = requiredProduction;
result.fullGrowth = fullGrowth;
result.fullProduction = fullProduction;
result.summary = sprintf(['%d of %d organisms; growth %.4g (needed %.4g), ' ...
    'production %.4g (needed %.4g)\n  kept:    %s\n  dropped: %s'], ...
    numel(present), numel(organisms), growth, requiredGrowth, ...
    production, requiredProduction, listOf(present), listOf(absent));
end


function [products, weights] = resolveProducts(model, given)
% Returns the product exchanges and their weights, checked against the model.
if isempty(given)
    products = cell(0, 1);
    weights = zeros(0, 1);
    return
end
if ischar(given)
    given = {given};
end
if isstruct(given)
    products = {given.rxn}';
    weights = reshape([given.weight], [], 1);
else
    products = given(:);
    weights = ones(numel(products), 1);
end
missing = products(findRxnIDs(model, products) == 0);
if ~isempty(missing)
    error('minimalMicrobiome:unknownProduct', ...
        ['No reaction %s in the community model. A product is named by its ' ...
         'community exchange reaction, which is the only place anything ' ...
         'leaves the system.'], missing{1});
end
end


function text = listOf(names)
% Returns a readable list, or 'none' for an empty one.
if isempty(names)
    text = 'none';
else
    text = strjoin(names(:)', ', ');
end
end
