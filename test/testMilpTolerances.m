% Tests that the mixed-integer solves ask for their own tolerances.
%
% The formulations used to call changeCobraSolverParams, which is a global
% COBRA setting: the tolerance stayed changed for everything the session
% did afterwards, plain optimizeCbModel included. They now hand the
% tolerances to the solve instead, so nothing outside the call moves.
%
% Run through runSpectraTests, which sets the solver and reports.

model = readCbModel('ecoli_core_model.mat');
n = numel(model.rxns);
core = findRxnIDs(model, {'PGI','PFK','FBA','TPI','GAPD','PGK','ENO','PYK'});

%% the three tolerances that matter are all asked for
params = milpTolerances();
names = params(1:2:end);
assert(numel(params) == 6, 'expected three name/value pairs, got %d entries', ...
    numel(params));
assert(all(ismember({'feasTol','optTol','intTol'}, names)), ...
    'milpTolerances does not set all three of feasTol, optTol and intTol');

%% the integrality tolerance is set, which is the one that is easily missed
% A binary is tied to its reaction by lb*z <= v <= ub*z, so a binary left
% at the tolerance and still read as zero permits a flux of tol*ub. At a
% solver's usual 1e-5 and the usual bound of 1000 that is 0.01, which is
% far above inclusionCutoff, so reactions the solver switched off would be
% counted as carrying flux and kept.
intAt = find(strcmp(names, 'intTol'));
assert(params{2*intAt} <= 1e-9, ...
    'the integrality tolerance is %g, which is loose enough to leak binaries', ...
    params{2*intAt});

%% the default is the tightest the solvers accept
assert(all(cell2mat(params(2:2:end)) == 1e-9), ...
    'the default tolerance is not 1e-9');

%% a tolerance below the floor is raised, not passed through and dropped
% Gurobi refuses any of the three below 1e-9 and keeps whatever it had,
% which is looser than either, so asking for less must not be obeyed.
previous = warning('off', 'milpTolerances:belowFloor');
restore = onCleanup(@() warning(previous));
lastwarn('');
tight = milpTolerances(1e-14);
[message, identifier] = lastwarn;
assert(strcmp(identifier, 'milpTolerances:belowFloor'), ...
    'asking for 1e-14 did not warn; it said ''%s''', message);
assert(all(cell2mat(tight(2:2:end)) == 1e-9), ...
    'a tolerance under the floor was passed through rather than raised');
clear restore

%% a looser tolerance is taken as given
loose = milpTolerances(1e-6);
assert(all(cell2mat(loose(2:2:end)) == 1e-6), ...
    'a tolerance above the floor was not used as given');

%% 'default' returns an empty cell array so session/solver defaults stay in place
assert(isempty(milpTolerances('default')), ...
    'milpTolerances(''default'') should return {}');
assert(isempty(lpTolerances('default')), ...
    'lpTolerances(''default'') should return {}');
defaultIDs = spectraCC(model, 1e-4, 'stoichiometry', 0, [], false, 'default');
assert(numel(defaultIDs) == 87, ...
    'spectraCC with solverTol=''default'' did not find 87 consistent reactions');

%% something that is not a tolerance is refused
assertTolFails(@() milpTolerances(0), 'positive', 'zero was accepted');
assertTolFails(@() milpTolerances(-1e-9), 'positive', 'a negative was accepted');
assertTolFails(@() milpTolerances([1e-9 1e-9]), 'scalar', 'a vector was accepted');
assertTolFails(@() milpTolerances('tight'), 'scalar', 'a char was accepted');

%% minNet leaves the global solver settings exactly as it found them
before = readMilpParams();
prepared = prepare(model);
growing = changeRxnBounds(prepared, prepared.rxns{prepared.c ~= 0}, 0.1, 'l');
[~, ~, stat] = minNet(growing, zeros(n,1), ones(n,1), 1e-4, 1, 'minNetMILP', 60);
assert(stat == 1, 'the MILP did not solve (stat %d)', stat);
after = readMilpParams();
assert(isequal(before, after), ...
    ['minNet changed the session''s MILP solver settings: feasTol %g -> %g, ' ...
     'optTol %g -> %g, intTol %g -> %g'], before(1), after(1), before(2), ...
    after(2), before(3), after(3));

%% tradeOff leaves them alone too
before = readMilpParams();
[~, ~, tradeStat] = tradeOff(growing, zeros(n,1), ones(n,1), 1e-4, 1, 60);
assert(tradeStat == 1, 'the tradeOff MILP did not solve (stat %d)', tradeStat);
assert(isequal(readMilpParams(), before), ...
    'tradeOff changed the session''s MILP solver settings');

%% so does a whole extraction, which is what a caller actually runs
before = readMilpParams();
spectraME(model, core, 1e-4, 'stoichiometry', ones(n,1), 1, 'coreDirection', ...
    'minNetMILP', 60, false, {}, 1);
assert(isequal(readMilpParams(), before), ...
    'spectraME changed the session''s MILP solver settings');

%% the tolerance reaches the formulation, rather than being accepted and ignored
% A value the helper refuses must still be refused when it arrives through
% minNet and through spectraME, which is what shows it is threaded.
assertTolFails(@() minNet(growing, zeros(n,1), ones(n,1), 1e-4, 1, ...
    'minNetMILP', 60, [], {}, {}, [], -1), 'positive', ...
    'minNet accepted a tolerance it should have refused');
assertTolFails(@() tradeOff(growing, zeros(n,1), ones(n,1), 1e-4, 1, 60, ...
    [], {}, {}, [], -1), 'positive', ...
    'tradeOff accepted a tolerance it should have refused');
assertTolFails(@() spectraME(model, core, 1e-4, 'stoichiometry', ones(n,1), ...
    1, 'coreDirection', 'minNetMILP', 60, false, {}, 1, {}, [], -1), ...
    'positive', 'spectraME accepted a tolerance it should have refused');

%% a loosened tolerance still solves
[kept, ~, looseStat] = minNet(growing, zeros(n,1), ones(n,1), 1e-4, 1, ...
    'minNetMILP', 60, [], {}, {}, [], 1e-7);
assert(looseStat == 1, 'the MILP did not solve at 1e-7 (stat %d)', looseStat);
assert(sum(kept) > 0, 'the loosened solve kept nothing');

%% lpTolerances asks for feasTol and optTol without mutating global LP settings
lpParams = lpTolerances();
assert(isequal(lpParams, {'feasTol', 1e-9, 'optTol', 1e-9}), ...
    'lpTolerances did not return default 1e-9 feasTol/optTol');
[lpFeasBefore, lpOptBefore] = getCobraSolverParams('LP', {'feasTol','optTol'}, struct());
spectraCC(model, 1e-4, 'stoichiometry', 0);
spectraME(model, core, 1e-4, 'stoichiometry', ones(n,1), 1, 'coreDirection', ...
    'minNetLP', 60, false, {}, 1);
[lpFeasAfter, lpOptAfter] = getCobraSolverParams('LP', {'feasTol','optTol'}, struct());
assert(isequal([lpFeasBefore, lpOptBefore], [lpFeasAfter, lpOptAfter]), ...
    'LP solves mutated global CBT_LP_PARAMS');

%% MILP warm starts run automatically across empty-core, one-way-core, and tradeOff
growingZeroC = growing;
growingZeroC.c(:) = 0;
[keptZeroC, lpsZeroC] = spectraME(growingZeroC, [], 1e-4, 'stoichiometry', ...
    ones(n, 1), 1, [], 'minNetMILP', 60, false, {}, 1);
assert(lpsZeroC > 0 && numel(keptZeroC.rxns) > 0, ...
    'empty-core minNetMILP with forced growth and c=0 should run LP warm start');

irrevCore = findRxnIDs(model, {'PFK'});
[keptIrrev, lpsIrrev] = spectraME(model, irrevCore, 1e-4, 'stoichiometry', ...
    ones(n, 1), 1, [], 'minNetMILP', 60, false, {}, 1);
assert(lpsIrrev > 0 && ismember('PFK', keptIrrev.rxns), ...
    'irreversible-core minNetMILP should run sparse LP warm start');

signedW = -ones(n, 1);
signedW(irrevCore) = 0;
[keptTrade, lpsTrade] = spectraME(model, irrevCore, 1e-4, 'stoichiometry', ...
    signedW, 1, [], 'tradeOff', 60, false, {}, 1);
assert(lpsTrade > 0 && ismember('PFK', keptTrade.rxns), ...
    'signed tradeOff should run sparse LP warm start and polish');


function values = readMilpParams()
% Returns the session's MILP tolerances, globals included.
%
% The struct argument is what makes getCobraSolverParams consult the global
% CBT_MILP_PARAMS at all; without it the defaults come back whatever the
% session has done, which would make this check pass vacuously.
[feasTol, optTol, intTol] = getCobraSolverParams('MILP', ...
    {'feasTol','optTol','intTol'}, struct());
values = [feasTol, optTol, intTol];
end


function model = prepare(model)
% Reorients the model the way spectraME hands it to a formulation.
n = numel(model.rxns);
backwards = model.ub <= 0;
swap = model.ub(backwards);
model.S(:, backwards) = -model.S(:, backwards);
model.ub(backwards) = -1 * model.lb(backwards);
model.lb(backwards) = -1 * swap;
model.rev = ones(n, 1);
model.rev(model.lb >= 0) = 0;
end


function assertTolFails(call, fragment, message)
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
