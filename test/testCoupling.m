% Tests that the formulations honour a model's coupling constraints.
%
% Every formulation builds its own problem out of model.S, so a constraint
% that lives in model.C is not in the problem unless it is put there. A
% community model ties each unit's reactions to its own anchor through C,
% and an extraction that ignored C would be free to run a unit's metabolism
% while the unit is not growing -- which is exactly what the coupling is
% there to forbid.
%
% Run through runSpectraTests, which sets the solver and reports.

producer = couplingOrganism('producer', 'glc', 'ac');
consumer = couplingOrganism('consumer', 'ac', 'co2');
glucose = struct('met', {'glc[e]'}, 'lb', {-10}, 'ub', {1000});
community = buildCommunityModel({producer, consumer}, {'P', 'C'}, ...
    'poolMedium', glucose);

shut = community.model;
shut.ub(findRxnIDs(shut, 'biomass__C')) = 0;
loose = rmfield(shut, {'C', 'd', 'dsense', 'ctrs'});
n = numel(shut.rxns);

%% the model under test really does carry coupling constraints
assert(isfield(shut, 'C') && size(shut.C, 1) > 0, ...
    'the community was built without coupling, so there is nothing to test');
assert(size(shut.C, 2) == n, 'model.C does not have one column per reaction');

%% couplingRows returns the constraints, padded to the problem's width
[A, b, csense] = couplingRows(shut, 7);
assert(size(A, 1) == size(shut.C, 1), 'couplingRows lost a constraint');
assert(size(A, 2) == n + 7, 'couplingRows did not pad to the problem width');
assert(all(all(A(:, 1:n) == shut.C)), 'couplingRows changed the constraints');
assert(all(A(:, n+1:end) == 0, 'all'), ...
    'couplingRows put something in the auxiliary columns');
assert(numel(b) == size(A, 1) && numel(csense) == size(A, 1), ...
    'couplingRows returned a right-hand side of the wrong length');

%% a model with no coupling gives no rows, so nothing else changes
[none, noneB] = couplingRows(loose, 3);
assert(isempty(none) && size(none, 2) == n + 3, ...
    'a model without coupling was given constraints');
assert(isempty(noneB), 'a model without coupling was given a right-hand side');

%% an unusable C is refused rather than silently misread
broken = shut;
broken.C = broken.C(:, 1:end-1);
assertCouplingFails(@() couplingRows(broken, 0), 'columns', ...
    'a C of the wrong width was accepted');
broken = shut;
broken = rmfield(broken, 'dsense');
assertCouplingFails(@() couplingRows(broken, 0), 'dsense', ...
    'a C without a dsense was accepted');

%% the extraction refuses to force flux through a unit that cannot grow
% With biomass__C held at zero, coupling caps every reaction of C at
% couplingThreshold, which is 0.01. Asking a reaction of C to carry 0.1
% is therefore impossible, and the direction phase must say so rather
% than keep trying.
core = findRxnIDs(shut, 'TR_ac_u__C');
assertCouplingFails(@() spectraME(shut, core, 0.1, 'stoichiometry', ...
    ones(n, 1), 1, 'coreDirection', 'minNetMILP', 30, false, {}, 1), ...
    'core reactions', ...
    'the extraction forced flux through a unit the coupling has pinned');

%% the same problem without the coupling is solvable, so the test bites
% If this failed for some other reason the check above would prove nothing.
extracted = spectraME(loose, findRxnIDs(loose, 'TR_ac_u__C'), 0.1, ...
    'stoichiometry', ones(n, 1), 1, 'coreDirection', 'minNetMILP', 30, ...
    false, {}, 1);
assert(ismember('TR_ac_u__C', extracted.rxns), ...
    ['without the coupling constraints the same core reaction is ' ...
     'reachable, so the refusal above is the coupling and not the setup']);

%% the consistency check honours them too
% TR_ac_u__C can carry 0.01, which is under the threshold, so with the
% coupling in force it reads as blocked and without it does not.
withCoupling = spectraCC(shut, 0.1, 'stoichiometry');
without = spectraCC(loose, 0.1, 'stoichiometry');
assert(~ismember(core, withCoupling), ...
    'spectraCC called a coupling-pinned reaction consistent');
assert(ismember(findRxnIDs(loose, 'TR_ac_u__C'), without), ...
    ['the same reaction is consistent without the coupling, so the check ' ...
     'above is the coupling and not the setup']);

%% reorienting a reaction flips it in C as well as in S
% spectraME and spectraCC reorient every reaction whose upper bound is at
% or below zero, flipping its column in S. A constraint in C is written
% over the same flux, so leaving C alone would silently reverse what the
% constraint says about that reaction. The coupling a community builds is
% symmetric and would not notice; a constraint a caller adds need not be,
% and minimalMicrobiome adds exactly such a one.
ecoli = readCbModel('ecoli_core_model.mat');
reversible = find(ecoli.lb < 0 & ecoli.ub > 0, 1);
assert(~isempty(reversible), 'the core model has no reversible reaction');
oneWay = ecoli;
oneWay.ub(reversible) = 0;      % reverse only, so the extraction reorients it

forced = addCOBRAConstraints(oneWay, ecoli.rxns(reversible), 0.5, ...
    'c', 1, 'dsense', 'G', 'ConstraintID', 'forceForward');
assert(isempty(spectraCC(forced, 1e-4, 'stoichiometry')), ...
    ['a constraint written over a reaction the check reorients was ' ...
     'reversed, so a flux the bounds forbid looked achievable']);

% The same model without the constraint is perfectly workable, so the
% check above is the constraint and not the bounds.
assert(~isempty(spectraCC(oneWay, 1e-4, 'stoichiometry')), ...
    'the fixture is unworkable for some reason other than the constraint');


function model = couplingOrganism(name, eats, makes)
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


function assertCouplingFails(call, fragment, message)
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
