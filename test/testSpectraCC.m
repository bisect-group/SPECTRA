% Tests for the SPECTRA consistency check.
%
% Run through runSpectraTests, which sets the solver and reports.

model = readCbModel('ecoli_core_model.mat');

%% the textbook model has 87 consistent reactions of 95
ids = spectraCC(model, 1e-4);
assert(numel(ids) == 87, ...
    'expected 87 consistent reactions, got %d', numel(ids));

%% calls written against the earlier signatures still work
% detectionCutoff was added after seed, and seed after consistencyType,
% precisely so that these keep working.
ids2 = spectraCC(model, 1e-4);
ids3 = spectraCC(model, 1e-4, 'stoichiometry');
assert(isequal(ids2, ids3), ...
    'the two- and three-argument calls disagree');

%% a seed makes the result repeat
a = spectraCC(model, 1e-4, 'stoichiometry', 0);
b = spectraCC(model, 1e-4, 'stoichiometry', 0);
assert(isequal(a, b), 'the same seed gave two different answers');

%% the detection cutoff is honoured
% Passing the default explicitly must not change anything.
withDefault = spectraCC(model, 1e-4, 'stoichiometry', 0, 0.99*1e-4);
assert(isequal(withDefault, a), ...
    'passing the default detectionCutoff changed the answer');

%% topology never finds fewer reactions than stoichiometry
steady = spectraCC(model, 1e-4, 'stoichiometry', 0);
topo = spectraCC(model, 1e-4, 'topology', 0);
assert(numel(topo) >= numel(steady), ...
    'topology found %d consistent reactions against stoichiometry''s %d; it is the weaker requirement and cannot find fewer', ...
    numel(topo), numel(steady));

%% an unrecognised consistencyType is refused
caught = false;
try
    spectraCC(model, 1e-4, 'stoichometry');   % deliberate typo
catch err
    caught = strcmp(err.identifier, 'spectraCC:badConsistencyType');
end
assert(caught, 'a misspelt consistencyType was not refused');

%% the returned model is consistent in its own right
consistent = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), ids)), ...
    'metFlag', false);
again = spectraCC(consistent, 1e-4, 'stoichiometry', 1);
assert(numel(again) == numel(consistent.rxns), ...
    'the consistent model is not a fixed point: %d of %d survive a second check', ...
    numel(again), numel(consistent.rxns));
