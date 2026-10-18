% Tests for checkExtraction, which reports what an extraction did and did
% not deliver.
%
% Run through runSpectraTests, which sets the solver and reports.

model = readCbModel('ecoli_core_model.mat');
consistentIDs = spectraCC(model, 1e-4, 'stoichiometry', 0);
consistent = removeRxns(model, ...
    model.rxns(setdiff(1:numel(model.rxns), consistentIDs)), 'metFlag', false);

core = {'PGI','PFK','FBA'};
core = core(ismember(core, consistent.rxns));
extracted = spectraME(consistent, findRxnIDs(consistent, core), 1e-3, ...
    'stoichiometry', [], 1, 'coreDirection', 'minNetLP', [], 0, {}, 0);

%% a sound extraction reports valid
report = checkExtraction(extracted, core, 1e-3);
assert(isempty(report.missingCore), ...
    'core reactions went missing: %s', strjoin(report.missingCore', ', '));
assert(report.isValid, 'the extraction should have been valid: %s', report.summary);
assert(report.isConsistent, 'the extracted model should be flux consistent');
assert(report.nReactions == numel(extracted.rxns), ...
    'the reaction count does not match the model');

%% a missing core reaction is noticed
% Ask after a reaction the extraction was never given.
absent = setdiff(consistent.rxns, extracted.rxns);
if ~isempty(absent)
    report = checkExtraction(extracted, [core, absent(1)'], 1e-3);
    assert(ismember(absent{1}, report.missingCore), ...
        'a core reaction absent from the model was not reported missing');
    assert(~report.isValid, 'a report with a missing core reaction is not valid');
end

%% indices are refused, since they do not survive the extraction
% The extracted model is smaller, so indices into the model it came from
% point at different reactions.
caught = false;
try
    checkExtraction(extracted, [1 2 3], 1e-3);
catch err
    caught = strcmp(err.identifier, 'checkExtraction:numericCore');
end
assert(caught, 'numeric core reactions should be refused');

%% an unknown method is refused
caught = false;
try
    checkExtraction(extracted, core, 1e-3, 'stoichiometry', [], [], 'guesswork');
catch err
    caught = strcmp(err.identifier, 'checkExtraction:badMethod');
end
assert(caught, 'an unknown method should be refused');

%% fva cannot answer a topology question, and says so
% find_blocked-style FVA always assumes a steady state.
caught = false;
try
    checkExtraction(extracted, core, 1e-3, 'topology', [], [], 'fva');
catch err
    caught = strcmp(err.identifier, 'checkExtraction:fvaCannotDoTopology');
end
assert(caught, 'fva should refuse a topology check rather than answer the wrong question');

%% the spectra method runs and agrees that a sound extraction is valid
report = checkExtraction(extracted, core, 1e-3, 'stoichiometry', [], 0, 'spectra');
assert(report.isValid, ...
    'the spectra method disagreed with fva on a sound extraction: %s', report.summary);
