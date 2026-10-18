% Tests for gap-filling towards growth and towards metabolic tasks.
%
% Run through runSpectraTests, which sets the solver and reports.

universal = readCbModel('ecoli_core_model.mat');
biomass = universal.rxns{universal.c ~= 0};

SALTS  = {'EX_pi(e)','EX_h2o(e)','EX_h(e)','EX_nh4(e)','EX_o2(e)','EX_co2(e)'};
CARBON = {'EX_glc(e)','EX_ac(e)','EX_succ(e)'};

media = struct('label',{},'rxns',{},'lb',{});
for i = 1:numel(CARBON)
    media(i).label = CARBON{i};
    media(i).rxns  = [SALTS, CARBON(i)];
    media(i).lb    = [-1000*ones(1,numel(SALTS)), -10];
end

%% a draft missing four reactions is repaired, and grows everywhere
knock = {'ENO','ACONTb','SUCDi','NH4t'};
draft = removeRxns(universal, knock, 'metFlag', false);

result = gapfillForGrowth(draft, universal, media, 0.1, biomass);
assert(isempty(result.unsatisfied), ...
    'the gap-filled model still fails in %d medium/media', numel(result.unsatisfied));
assert(numel(result.satisfied) == numel(media), ...
    'expected to grow in all %d media', numel(media));
assert(all(ismember(result.added, knock)), ...
    'reactions were added that had not been removed: %s', ...
    strjoin(setdiff(result.added, knock)', ', '));

%% what it adds is in the result, and the draft is still there
assert(all(ismember(draft.rxns, result.model.rxns)), ...
    'the gap-filled model lost reactions the draft had');
assert(all(ismember(result.added, result.model.rxns)), ...
    'a reaction reported as added is not in the model');

%% the answer does not move with tol
% Growth is demanded through a bound on biomass, not through tol, which is
% what keeps the two independent.
counts = zeros(1,3);
tolerances = [1e-3 1e-4 1e-5];
for k = 1:numel(tolerances)
    r = gapfillForGrowth(draft, universal, media, 0.1, biomass, tolerances(k));
    counts(k) = numel(r.added);
end
assert(numel(unique(counts)) == 1, ...
    'the number of additions moved with tol: %s', mat2str(counts));

%% tradeOff is refused for growth
caught = false;
try
    gapfillForGrowth(draft, universal, media, 0.1, biomass, 1e-4, 'tradeOff');
catch err
    caught = strcmp(err.identifier, 'gapfillForGrowth:badProbType');
end
assert(caught, 'tradeOff should be refused when gap-filling towards growth');

%% a medium naming an unknown exchange is refused
bad = media(1);
bad.rxns{end+1} = 'EX_unobtainium(e)';
bad.lb(end+1) = -10;
caught = false;
try
    gapfillForGrowth(draft, universal, bad, 0.1, biomass);
catch err
    caught = strcmp(err.identifier, 'gapfillForGrowth:unknownExchange');
end
assert(caught, 'a medium naming a reaction the model lacks was accepted');

%% gap-filling towards tasks repairs what the tasks need
equation = taskEquation({'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'}, ...
                        [-1 -1 1 1 1], [1 1000]);
AER = {'glc-D[e]','o2[e]','pi[e]','h2o[e]','h[e]'};
OUT = {'co2[e]','h2o[e]','h[e]','ac[e]','etoh[e]','for[e]'};
OUT = OUT(ismember(OUT, universal.mets));

t1 = metabolicTask('atp_aerobic', 'inputs', AER, 'outputs', OUT, ...
                   'equations', equation);
t2 = metabolicTask('atp_anaerobic', 'inputs', AER(~strcmp(AER,'o2[e]')), ...
                   'outputs', OUT, 'equations', equation);
tasks = [t1; t2];

taskKnock = {'CYTBD','ACALD','GLUSy','AKGDH'};
taskDraft = removeRxns(universal, taskKnock, 'metFlag', false);
assert(~all([checkTasks(taskDraft, tasks, 0).ok]), ...
    'the draft was supposed to fail at least one task');

filled = gapfillForTasks(taskDraft, universal, tasks, 1e-3);
assert(isempty(filled.unsatisfied), ...
    'the gap-filled model still fails: %s', strjoin(filled.unsatisfied', ', '));
assert(all(ismember(filled.added, universal.rxns)), ...
    'a reaction was added that is not in the universal model; task scaffolding leaked');

%% a task the universal model cannot perform is skipped, not attempted
impossible = metabolicTask('impossible', 'outputs', OUT, ...
    'equations', equation);   % nothing may be taken up
withImpossible = [tasks; impossible];
filled2 = gapfillForTasks(taskDraft, universal, withImpossible, 1e-3);
assert(ismember('impossible', filled2.skipped), ...
    'a task the universal model cannot perform should have been skipped');
