% Tests for reading and reshaping a community after it has been built.
%
% Run through runSpectraTests, which sets the solver and reports.

producer = helperOrganism('producer', 'glc', 'ac');
consumer = helperOrganism('consumer', 'ac', 'co2');
glucose = struct('met', {'glc[e]'}, 'lb', {-10}, 'ub', {1000});
community = buildCommunityModel({producer, consumer}, {'P', 'C'}, ...
    'poolMedium', glucose);

%% a reaction is traced back to the unit it belongs to
assert(strcmp(communityOrganismOf(community, 'biomass__P'), 'P'), ...
    'a tagged reaction was not traced to its unit');
assert(isempty(communityOrganismOf(community, 'EX_glc[u]')), ...
    'a pool exchange was claimed by a unit');
assert(strcmp(communityOrganismOf(community, 'TR_ac_u__P'), 'P'), ...
    'a transport was not traced to its unit');

%% a metabolite wears its tag inside the compartment and is still traced
assert(strcmp(communityOrganismOf(community, 'ac[e__C]'), 'C'), ...
    'a tagged metabolite was not traced to its unit');
assert(isempty(communityOrganismOf(community, 'ac[u]')), ...
    'a pool metabolite was claimed by a unit');

%% a name that merely ends in a tag is not mistaken for one
assert(isempty(communityOrganismOf(community, 'something__elsewhere')), ...
    'a reaction tagged by a name that is not a unit was claimed by one');

%% a cell array is answered elementwise
traced = communityOrganismOf(community, {'UP__C', 'EX_glc[u]'});
assert(isequal(traced, {'C', ''}), ...
    'a cell array of identifiers was not answered elementwise');

%% the coupling summary names every unit and its anchor
text = couplingSummary(community);
assert(contains(text, 'biomass__P') && contains(text, 'biomass__C'), ...
    'the coupling summary does not name both anchors');
assert(contains(text, 'objective'), ...
    'the coupling summary does not say how the anchor was decided');

%% an uncoupled community says so rather than reporting a cap
loose = buildCommunityModel({producer, consumer}, {'P', 'C'}, ...
    'poolMedium', glucose, 'couple', false, 'check', false);
assert(strcmp(couplingSummary(loose), 'not coupled'), ...
    'an uncoupled community was given a coupling table');

%% subsetting keeps the pool, so the survivors still face the environment
alone = communitySubset(community, {'P'});
assert(isequal(alone.organisms, {'P'}), 'the subset kept the wrong units');
assert(~any(contains(alone.model.rxns, '__C')), ...
    'the dropped unit''s reactions are still in the model');
assert(ismember('EX_glc[u]', alone.model.rxns), ...
    'the subset lost the pool exchange it needs to eat');
fed = optimizeCbModel(alone.model);
assert(fed.stat == 1 && fed.f > 1e-6, ...
    'the surviving unit cannot grow on its own (f = %g)', fed.f);

%% subsetting moves the objective onto the survivors
assert(sum(alone.model.c ~= 0) == 1, ...
    'the subset did not leave exactly one objective reaction');
assert(strcmp(alone.model.rxns{alone.model.c ~= 0}, 'biomass__P'), ...
    'the subset optimises something other than the survivor''s anchor');

%% a unit the community does not have is refused
assertCommunityFails(@() communitySubset(community, {'P', 'nobody'}), ...
    'Not a unit', 'a unit that does not exist was accepted');
assertCommunityFails(@() communitySubset(community, {}), ...
    'no units at all', 'an empty subset was accepted');

%% narrowing to a reduced model keeps the bookkeeping honest
cut = removeRxns(community.model, {'CONV__C'});
narrowed = communityWithModel(community, cut);
assert(~ismember('CONV__C', narrowed.units(2).rxns), ...
    'a reaction the reduction dropped is still listed against its unit');
assert(numel(narrowed.units(1).rxns) == numel(community.units(1).rxns), ...
    'narrowing touched a unit the reduction left alone');
assert(all(ismember(narrowed.communityExchanges, cut.rxns)), ...
    'the narrowed exchanges name reactions the model no longer has');

%% narrowing keeps the anchors, even where one did not survive
beheaded = communityWithModel(community, ...
    removeRxns(community.model, {'biomass__C'}));
assert(isequal(beheaded.biomassReactions, community.biomassReactions), ...
    ['narrowing forgot which reaction anchors a unit, which hides the ' ...
     'news that the reduction dropped it']);

%% decomposing gives back one model per unit, untagged
units = communityDecompose(community);
assert(numel(units) == 2, 'decompose did not return one model per unit');
assert(isequal({units.name}, {'P', 'C'}), ...
    'the units came back in the wrong order');
assert(all(ismember({'UP', 'CONV', 'biomass'}, units(1).model.rxns)), ...
    'the decomposed unit kept its tags');
assert(strcmp(units(1).anchor, 'biomass'), ...
    'the anchor was not reported untagged');
assert(~any(contains(units(1).model.rxns, '__P')), ...
    'a reaction of the decomposed unit is still tagged');
assert(~any(contains(units(1).model.mets, '__P')), ...
    'a metabolite of the decomposed unit is still tagged');

%% a transport becomes an exchange, carrying its own bounds
% Dropping the pool side of 'x[e] -> x[u]' leaves 'x[e] ->', so the rates
% measured for the unit in the community survive the split.
assert(ismember('TR_glc_u', units(1).model.rxns), ...
    'the transport did not survive the split as a boundary reaction');
single = sum(units(1).model.S ~= 0, 1) == 1;
at = strcmp(units(1).model.rxns, 'TR_glc_u');
assert(single(at), 'the transport was not left with a single metabolite');

%% the split models stand on their own
own = changeRxnBounds(units(1).model, 'TR_glc_u', -10, 'l');
grown = optimizeCbModel(own);
assert(grown.stat == 1 && grown.f > 1e-6, ...
    'the decomposed unit cannot grow on its own (f = %g)', grown.f);

%% decomposing a shared-compartment community adds the exchanges instead
% There are no transports to turn into exchanges there, so one is added
% for every pool metabolite the unit touches.
flat = buildCommunityModel({producer, consumer}, {'P', 'C'}, ...
    'mode', 'shared', 'poolMedium', glucose);
flatUnits = communityDecompose(flat);
assert(any(contains(flatUnits(1).model.rxns, 'EX_glc')), ...
    'the shared-compartment unit was given no exchange for what it eats');
assert(sum(sum(flatUnits(1).model.S ~= 0, 1) == 1) >= 2, ...
    'the shared-compartment unit has nothing entering or leaving it');

%% the helpers refuse anything that is not a community
assertCommunityFails(@() couplingSummary(struct('model', 1)), ...
    'buildCommunityModel returns', 'a bare struct was accepted');
assertCommunityFails(@() communityDecompose(producer), ...
    'buildCommunityModel returns', 'a plain model was accepted');


function model = helperOrganism(name, eats, makes)
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


function assertCommunityFails(call, fragment, message)
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
