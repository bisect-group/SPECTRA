% Tests for reducing a community to the organisms that still do the job.
%
% Run through runSpectraTests, which sets the solver and reports.

p1 = microbiomeOrganism('p1', 'glc', 'ac');
p2 = microbiomeOrganism('p2', 'glc', 'ac');
cleaner = microbiomeOrganism('d', 'ac', 'co2');
glucose = struct('met', {'glc[e]'}, 'lb', {-10}, 'ub', {1000});
community = buildCommunityModel({p1, p2, cleaner}, {'p1', 'p2', 'd'}, ...
    'poolMedium', glucose);

%% the whole community grows, so there is something to reduce
whole = optimizeCbModel(community.model);
assert(whole.stat == 1 && whole.f > 1e-6, ...
    'the fixture community cannot grow (f = %g)', whole.f);

%% a redundant organism is dropped
% Glucose is the limit, not the number of mouths, so one glucose eater
% achieves what three organisms do.
result = minimalMicrobiome(community, 'seed', 1);
assert(numel(result.present) < 3, ...
    'nothing was dropped although the glucose supply caps the growth');
assert(numel(result.present) + numel(result.absent) == 3, ...
    'the membership does not account for every organism');
assert(result.growth >= result.requiredGrowth - 1e-6, ...
    'the minimal community grows at %g, under the %g it had to reach', ...
    result.growth, result.requiredGrowth);

%% the requirement is a fraction of what was measured, not a bare number
assert(abs(result.requiredGrowth - 0.8 * result.fullGrowth) < 1e-9, ...
    'the growth requirement is not the stated fraction of what was measured');

%% the result is a community in its own right
assert(isequal(result.community.organisms, result.present(:)), ...
    'the reduced community holds different organisms from the ones reported');
grown = optimizeCbModel(result.community.model);
assert(grown.stat == 1 && grown.f > 1e-6, ...
    'the minimal community cannot grow on its own (f = %g)', grown.f);
assert(~any(contains(result.community.model.rxns, ...
    ['__' result.absent{1}])), ...
    'a dropped organism''s reactions are still in the reduced model');

%% the membership covers every organism, present or not
assert(numel(result.membership) == 3, ...
    'the membership does not have one entry per organism');
assert(isequal({result.membership.name}, {'p1', 'p2', 'd'}), ...
    'the membership is not in the order the community reports');
assert(sum([result.membership.present]) == numel(result.present), ...
    'the membership and the present list disagree');

%% a product requirement keeps the organism that makes it
% Only the cleaner makes carbon dioxide, so preserving it forces it in,
% and something to feed it has to stay too.
withProduct = minimalMicrobiome(community, 'products', {'EX_co2[u]'}, ...
    'seed', 1);
assert(ismember('d', withProduct.present), ...
    'the only producer of the metabolite to be preserved was dropped');
assert(withProduct.fullProduction > 0, ...
    'the full community was measured as producing nothing to preserve');
assert(withProduct.production >= withProduct.requiredProduction - 1e-6, ...
    'the minimal community produces %g, under the %g it had to reach', ...
    withProduct.production, withProduct.requiredProduction);

%% a weighted product is accepted in the same place
weighted = minimalMicrobiome(community, ...
    'products', struct('rxn', {'EX_co2[u]'}, 'weight', {2}), 'seed', 1);
assert(abs(weighted.fullProduction - 2 * withProduct.fullProduction) < 1e-6, ...
    'doubling a product''s weight did not double the measured production');

%% an organism named as required is kept whatever the objective says
held = minimalMicrobiome(community, 'required', {'p2'}, 'seed', 1);
assert(ismember('p2', held.present), 'a required organism was dropped');

%% requiring everything leaves everything
all3 = minimalMicrobiome(community, 'required', {'p1', 'p2', 'd'}, 'seed', 1);
assert(numel(all3.present) == 3, ...
    'requiring every organism did not keep every organism');

%% an organism the community does not have is refused
assertMicrobiomeFails(@() minimalMicrobiome(community, 'required', {'nobody'}), ...
    'Not an organism', 'an organism that does not exist was required');

%% a product that is not a reaction of the community is refused
assertMicrobiomeFails(@() minimalMicrobiome(community, 'products', {'EX_nothing[u]'}), ...
    'community exchange', 'a product that does not exist was accepted');

%% a community that cannot grow says so rather than returning nothing
starved = buildCommunityModel({p1, p2, cleaner}, {'p1', 'p2', 'd'}, ...
    'check', false);
assertMicrobiomeFails(@() minimalMicrobiome(starved), 'cannot grow', ...
    'a community with no medium was reduced rather than refused');

%% an uncoupled community is flagged, since absence then means nothing
loose = buildCommunityModel({p1, p2, cleaner}, {'p1', 'p2', 'd'}, ...
    'poolMedium', glucose, 'couple', false, 'check', false);
previous = warning('error', 'minimalMicrobiome:notCoupled');
cleanup = onCleanup(@() warning(previous));
assertMicrobiomeFails(@() minimalMicrobiome(loose, 'seed', 1), 'coupling', ...
    'an uncoupled community was reduced without comment');
clear cleanup


function model = microbiomeOrganism(name, eats, makes)
% Returns a one-step organism: eats one metabolite, secretes another.
model = createModel();
model = addReaction(model, 'UP', 'reactionFormula', ...
    sprintf('%s[e] -> x[c]', eats), ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, 'CONV', 'reactionFormula', ...
    sprintf('x[c] -> %s[e]', makes), ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, 'biomass', 'reactionFormula', 'x[c] -> ', ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, ['EX_' eats], 'reactionFormula', ...
    sprintf('%s[e] <=> ', eats), ...
    'lowerBound', -1000, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, ['EX_' makes], 'reactionFormula', ...
    sprintf('%s[e] <=> ', makes), ...
    'lowerBound', -1000, 'upperBound', 1000, 'printLevel', 0);
model = changeObjective(model, 'biomass');
model.modelID = name;
model.description = name;
end


function assertMicrobiomeFails(call, fragment, message)
% Asserts that a call fails with a message holding the given fragment.
try
    call();
catch err
    assert(contains(err.message, fragment), ...
        '%s: failed with ''%s'', which does not mention ''%s''', ...
        message, err.message, fragment);
    return
end
error('%s: the call was expected to fail and did not', message);
end
