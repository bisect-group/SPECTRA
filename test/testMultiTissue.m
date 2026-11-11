% Tests for replicating one model across tissues and joining them.
%
% Run through runSpectraTests, which sets the solver and reports.

gem = tissueModel('gem', 'glc', 'ac');
glucose = struct('met', {'glc[e]'}, 'lb', {-10}, 'ub', {1000});

%% one model is replicated across the tissues
% Tissues differ by their bounds and their pool, not by their network.
pools = struct('name', {'Bl', 'tis2_tis3'}, ...
    'units', {{'tis1', 'tis2'}, {'tis2', 'tis3'}});
body = buildMultiTissueModel(gem, {'tis1', 'tis2', 'tis3'}, pools, {'Bl'}, ...
    'anchorReactions', 'biomass', 'poolMedium', glucose, 'check', false);

assert(isequal(body.organisms, {'tis1'; 'tis2'; 'tis3'}), ...
    'the tissues were not reported in the order given');
assert(all(strcmp({body.units.anchorSource}, 'given')), ...
    'a stated anchor was not reported as given');
assert(isequal(body.biomassReactions, ...
    {'biomass__tis1'; 'biomass__tis2'; 'biomass__tis3'}), ...
    'the anchors were not tagged by their tissue');

%% only the pool that reaches the environment carries exchanges
assert(ismember('EX_glc[Bl]', body.model.rxns), ...
    'the blood has no exchange for glucose');
assert(~ismember('EX_glc[tis2_tis3]', body.model.rxns), ...
    'the interface pool was given an exchange, letting it bypass the blood');

%% a tissue off the blood reaches only the pool it belongs to
% tis3 never touches the blood, so it eats only what tis2 passes on.
third = body.units(3).rxns;
assert(~any(contains(third, '_Bl__tis3')), ...
    'the tissue off the blood was given a transport into it');
assert(any(contains(third, '_tis2_tis3__tis3')), ...
    'the tissue off the blood was given no transport into its own pool');

%% one model per tissue, when the tissues are different networks
liver = tissueModel('liver', 'glc', 'lac');
muscle = tissueModel('muscle', 'lac', 'co2');
twoPools = struct('name', {'Bl'}, 'units', {{'tis1', 'tis2'}});
pairBody = buildMultiTissueModel({liver, muscle}, {'tis1', 'tis2'}, ...
    twoPools, {'Bl'}, 'anchorReactions', {'biomass', 'biomass'}, ...
    'poolMedium', glucose, 'check', false);

assert(all(ismember({'glc[e__tis1]', 'lac[e__tis2]'}, pairBody.model.mets)), ...
    'each tissue did not keep its own metabolites');
% Each tissue kept its own network, not a copy of the first one.
assert(~ismember('glc[e__tis2]', pairBody.model.mets), ...
    'the second tissue was given the first one''s metabolites');
assert(~ismember('co2[e__tis1]', pairBody.model.mets), ...
    'the first tissue was given the second one''s metabolites');

%% the liver feeds the muscle through the blood
fed = optimizeCbModel(changeObjective(pairBody.model, 'biomass__tis2'));
assert(fed.stat == 1 && fed.f > 1e-6, ...
    'the muscle cannot grow on what the liver secretes (f = %g)', fed.f);

%% a model count that is not the tissue count is refused
try
    buildMultiTissueModel({gem, gem}, {'tis1', 'tis2', 'tis3'}, ...
        struct('name', {'Bl'}, 'units', {{'tis1', 'tis2', 'tis3'}}), {'Bl'});
    error('testMultiTissue:notRefused', ...
        'two models for three tissues were accepted');
catch err
    assert(contains(err.message, '2 models for 3 tissues'), ...
        'the wrong model count failed with ''%s''', err.message);
end

%% the pools and the environment must be given
try
    buildMultiTissueModel(gem, {'tis1', 'tis2'});
    error('testMultiTissue:notRefused', ...
        'a body was built without saying how its tissues meet');
catch err
    assert(contains(err.message, 'must all be given'), ...
        'the missing pools failed with ''%s''', err.message);
end

%% anchorReactions is checked against the tissue count
try
    buildMultiTissueModel(gem, {'tis1', 'tis2'}, twoPools, {'Bl'}, ...
        'anchorReactions', {'biomass'});
    error('testMultiTissue:notRefused', ...
        'one anchor for two tissues was accepted');
catch err
    assert(contains(err.message, '1 anchorReactions for 2 tissues'), ...
        'the wrong anchor count failed with ''%s''', err.message);
end

%% a tissue is anchored on maintenance rather than on growth
% Which is the usual case: a tissue is rarely growing.
upkeep = tissueModel('upkeep', 'glc', 'ac');
upkeep = addReaction(upkeep, 'ATPM', 'reactionFormula', 'x[c] -> ', ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
held = buildMultiTissueModel(upkeep, {'tis1', 'tis2'}, twoPools, {'Bl'}, ...
    'anchorReactions', 'ATPM', 'poolMedium', glucose, 'check', false);
assert(isequal(held.biomassReactions, {'ATPM__tis1'; 'ATPM__tis2'}), ...
    'the maintenance anchor was not used');

%% everything else is passed through to buildCommunityModel
loose = buildMultiTissueModel(gem, {'tis1', 'tis2'}, twoPools, {'Bl'}, ...
    'anchorReactions', 'biomass', 'couple', false, 'check', false, ...
    'modelId', 'body');
assert(isempty(loose.couplingFactor), 'couple was not passed through');
assert(strcmp(loose.model.modelID, 'body'), 'modelId was not passed through');


function model = tissueModel(name, eats, makes)
% Returns a one-step tissue: eats one metabolite, secretes another.
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
