function report = checkExtraction(extracted,coreRxns,tol,consType,detectionCutoff,seed,method)
% Check an extracted model against what the extraction promised.
%
% USAGE:
%   report = checkExtraction(extracted,coreRxns,tol,consType,detectionCutoff,seed,method)
%
% INPUTS:
%   extracted: COBRA model structure returned by spectraME or spectraCCME.
%   coreRxns:  Reaction IDs that were asked for, as a cell array of char.
%              spectraME takes core reactions as *indices into the model it
%              was given*, and the extracted model is smaller, so those
%              indices no longer point at the same reactions. Pass
%              model.rxns(coreRxns) rather than coreRxns itself. For
%              spectraCCME, leave out the reactions it reported as blocked;
%              it was never going to include those.
%
% OPTIONAL INPUTS:
%   tol:             Flux threshold the extraction was run with
%                    (Default: 1e-4). Use the same value, or the answer does
%                    not describe the model you have.
%   consType:        Consistency type the extraction was run with, either
%                    'stoichiometry' or 'topology'. Only used by
%                    method 'spectra'. (Default: 'stoichiometry')
%   detectionCutoff: Absolute flux at which a reaction counts as carrying
%                    flux. (Default: the solver's feasTol for method 'fva',
%                    0.99*tol for method 'spectra')
%   seed:            Seed for the consistency check's random coefficients.
%                    Only used by method 'spectra'. (Default: none)
%   method:          How to decide whether a reaction is blocked, either
%                    'fva' or 'spectra'. (Default: 'fva')
%                    'fva' maximises and minimises every reaction in turn.
%                    It is slow but authoritative. 'spectra' reuses this
%                    package's own consistency check, which is far quicker
%                    but *can report a dead reaction as consistent* on a
%                    model holding numerically marginal reactions - exactly
%                    the models an LP extraction produces. Use it only when
%                    speed matters more than the answer.
%
% OUTPUTS:
%   report: Structure with fields
%             nReactions   - how many reactions the extracted model has
%             missingCore  - core reactions absent from it. Should always be
%                            empty; an extraction that drops one is a bug.
%             blocked      - reactions in it that cannot carry flux in it
%             blockedCore  - the core reactions among blocked. These are the
%                            damaging ones: the reaction is present, so the
%                            model looks right, but it cannot carry flux and
%                            so cannot play the role it was chosen for.
%             isConsistent - true when nothing in the model is blocked
%             isValid      - true when every core reaction is present and
%                            able to carry flux, and nothing else is blocked
%             summary      - a readable multi-line summary, as char
%
% NOTE:
%   An extraction promises two things: that every core reaction is in the
%   result, and that the result is flux consistent. The first always holds.
%   The second does not always hold for the LP formulations. minNetLP and
%   growthOptim read their answer off a flux vector, keeping any reaction
%   whose flux clears a cutoff near the solver's own tolerance. A few
%   reactions end up carrying flux just under that cutoff while doing
%   load-bearing balancing work, and dropping them leaves a mass-balance
%   residual of the same order. That residual is above the solver's
%   feasibility tolerance, so the kept flux vector is not actually feasible
%   in the extracted model, and chains that relied on those sub-cutoff
%   reactions can then carry no flux at all - a core reaction included.
%   The mixed-integer formulations should not have this problem, since a
%   binary at zero forces its reaction's flux to exactly zero and so
%   discards nothing.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end
if ~exist('consType', 'var') || isempty(consType)
    consType = 'stoichiometry';
end
if ~exist('method', 'var') || isempty(method)
    method = 'fva';
end
if ~exist('seed', 'var')
    seed = [];
end

if ~any(strcmp(method,{'fva','spectra'}))
    error('checkExtraction:badMethod', ...
        'method must be ''fva'' or ''spectra'', not ''%s''.', method);
end
if strcmp(method,'fva') && ~strcmp(consType,'stoichiometry')
    % FVA always assumes a steady state, so it cannot answer the question
    % the accumulation condition asks. Saying so beats quietly checking
    % something other than what was requested.
    error('checkExtraction:fvaCannotDoTopology', ...
        ['method ''fva'' cannot check consType ''%s'', because flux ' ...
         'variability always assumes a steady state. Pass method ' ...
         '''spectra'' to check the accumulation condition.'], consType);
end

if ischar(coreRxns)
    coreRxns = {coreRxns};
end
if isnumeric(coreRxns)
    error('checkExtraction:numericCore', ...
        ['coreRxns must be reaction IDs, not indices. The extracted model ' ...
         'is smaller than the one you extracted from, so indices into it ' ...
         'no longer point at the same reactions. Pass model.rxns(coreRxns).']);
end
coreRxns = coreRxns(:);

present = ismember(coreRxns, extracted.rxns);
missingCore = coreRxns(~present);

if strcmp(method,'fva')
    if ~exist('detectionCutoff', 'var') || isempty(detectionCutoff)
        detectionCutoff = getCobraSolverParams('LP', 'feasTol');
    end
    % optPercentage 0 leaves the objective unconstrained, so this asks what
    % each reaction can carry at all rather than what it can carry at
    % optimal growth.
    [minF,maxF] = fluxVariability(extracted, 0, 'max', extracted.rxns, 0);
    isBlocked = (abs(minF) < detectionCutoff) & (abs(maxF) < detectionCutoff);
    blocked = extracted.rxns(isBlocked);
else
    if ~exist('detectionCutoff', 'var') || isempty(detectionCutoff)
        detectionCutoff = 0.99*tol;
    end
    consIDs = spectraCC(extracted, tol, consType, seed, detectionCutoff);
    isBlocked = true(numel(extracted.rxns),1);
    isBlocked(consIDs) = false;
    blocked = extracted.rxns(isBlocked);
end
blocked = blocked(:);
blockedCore = coreRxns(ismember(coreRxns, blocked));

report.nReactions   = numel(extracted.rxns);
report.missingCore  = missingCore;
report.blocked      = blocked;
report.blockedCore  = blockedCore;
report.isConsistent = isempty(blocked);
report.isValid      = isempty(missingCore) && isempty(blocked);

lines = {sprintf('%d reactions in the extracted model', report.nReactions)};
if isempty(missingCore)
    lines{end+1} = '  every core reaction is present';
else
    lines{end+1} = sprintf('  %d core reaction(s) MISSING: %s', ...
        numel(missingCore), strjoin(missingCore(1:min(5,end))', ', '));
end
if isempty(blocked)
    lines{end+1} = '  flux consistent: every reaction can carry flux';
else
    lines{end+1} = sprintf('  %d reaction(s) cannot carry flux', numel(blocked));
    if ~isempty(blockedCore)
        lines{end+1} = sprintf('  %d of those are CORE reactions: %s', ...
            numel(blockedCore), strjoin(blockedCore(1:min(5,end))', ', '));
    end
end
report.summary = strjoin(lines, newline);

if ~report.isValid
    warning('checkExtraction:notValid', ...
        ['The extracted model is not valid: %d core reaction(s) missing, ' ...
         '%d reaction(s) blocked, %d of those core.'], ...
        numel(missingCore), numel(blocked), numel(blockedCore));
end
end
