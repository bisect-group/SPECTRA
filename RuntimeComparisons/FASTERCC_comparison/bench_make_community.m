function [model, info] = bench_make_community(nSpecies, rep, cfg)
% bench_make_community  One AGORA2 community instance under random media.
%
%   [model, info] = bench_make_community(200, 1, cfg)
%
% Builds a community of nSpecies randomly chosen AGORA2 reconstructions with
% a shared [u] compartment, then imposes a RANDOM medium on the community
% exchanges. Everything is a pure function of (nSpecies, rep, cfg.baseSeed),
% so an instance is rebuildable from its name alone - the models are ~200 KB
% each and a 200-member community is ~210k reactions, far too large to store
% per instance.
%
% WHY MEDIA AND NOT REACTION REMOVAL
% The other two experiments make reactions inconsistent by deleting them,
% which is surgical and unlike anything a modeller does. Closing uptake
% reactions is how blocked content actually arises in practice: a reaction is
% consistent in isolation but dead under the medium you gave it. It also
% blocks in correlated clusters - a whole pathway goes dark when its one
% carbon source closes - rather than at random, which is the regime neither
% method has been tested in.
%
% THE MEDIUM
%   - each community uptake is opened with probability cfg.mediaOpenProb
%   - open uptakes get lb = -10^U(log10(minInflow), log10(maxInflow)),
%     log-uniform so the draw spans orders of magnitude rather than
%     clustering near the top of the range
%   - everything else is closed (lb = 0), which is what makes reactions
%     blocked
%   - secretion (ub) is left untouched; constraining only uptake is the
%     standard way a diet is applied
%
% The log-uniform floor matters. At cfg.mediaMinInflow = 1e-6 an open uptake
% can cap a pathway's flux at 1e-6, which is BELOW spectraCC's default
% detection cutoff of 0.99*epsilon. With epsilon = 1e-4 such a reaction is
% reported blocked though FVA calls it consistent - a tolerance artefact, not
% an algorithmic difference. That is why the community runs use epsilon=1e-8.
% info.minOpenInflow records the smallest inflow actually drawn so the margin
% is auditable per instance rather than assumed.

if nargin < 3 || isempty(cfg), cfg = bench_config('smoke'); end

t0 = tic;
info.seed = cfg.baseSeed + 7919*nSpecies + rep;
rng(info.seed, 'twister');

%% pick the species ------------------------------------------------------
% Sample from a FROZEN pool file, never from dir() of the model folder.
%
% dir() would make an instance depend on what had finished downloading when
% it ran: the same seed would build a different community before and after
% the download completed, and com200_rep01 would silently not be com200_rep01
% any more. The pool is written once by bench_freeze_agora_pool and is the
% definition of the universe these samples are drawn from.
poolFile = fullfile(fileparts(cfg.agoraDir), 'pool.txt');
assert(exist(poolFile,'file') > 0, 'bench_make_community:noPool', ...
    ['Model pool %s does not exist. Run bench_freeze_agora_pool once the ' ...
     'download has finished; sampling from a changing directory would make ' ...
     'instances irreproducible.'], poolFile);

pool = strtrim(string(readlines(poolFile)));
pool = pool(pool ~= "");
assert(numel(pool) >= nSpecies, 'bench_make_community:tooFew', ...
    'Asked for %d species but the pool holds only %d.', nSpecies, numel(pool));

pick  = sort(randsample(numel(pool), nSpecies));
names = cellstr(pool(pick));
info.poolSize = numel(pool);

modelCell = cell(1, nSpecies);
spBm      = cell(nSpecies, 1);
spATPM    = cell(nSpecies, 1);
for i = 1:nSpecies
    S = load(fullfile(cfg.agoraDir, names{i}));
    f = fieldnames(S);
    m = S.(f{1});

    bm = agora_biomass_rxn(m);
    assert(~isempty(bm), 'bench_make_community:noBiomass', ...
        '%s has no identifiable biomass reaction.', names{i});

    modelCell{i} = m;
    spBm{i}      = bm;
    if any(strcmp(m.rxns, 'DM_atp_c_')), spATPM{i} = 'DM_atp_c_';
    else,                                spATPM{i} = ''; end
end

%% join ------------------------------------------------------------------
% couple = FALSE. createCommModel defaults it true, matching
% build_community_model in spectra-cobra, which bounds every reaction of an
% organism by its own biomass flux. That is right for growth simulation and
% wrong here: it adds constraints that can make a reaction blocked which the
% uncoupled network carries fine, so a consistency benchmark run with
% coupling on would be measuring the coupling, not the network.
%
% check = FALSE for the same reason - it only reports whether the coupling
% factor would throttle an organism, and there is no coupling to report on.
% It also costs one LP per organism, which would be charged to no arm.
%
% The medium is applied below rather than through options.poolMedium: the
% draw is random over whatever community exchanges the join produces, and
% those are not known until it has run.
options = struct('spBm', {spBm}, 'spATPM', {spATPM}, 'sepUtEx', false, ...
                 'couple', false, 'check', false, 'verbFlag', false);
[model, infoCom] = createCommModel(modelCell, options); %#ok<ASGLU>

%% random medium ---------------------------------------------------------
% Community uptakes are the [u] exchanges createCommModel adds. Match on the
% compartment tag rather than on an organism suffix, so this does not depend
% on how the joiner happens to name per-organism reactions.
isComEx = startsWith(model.rxns, 'EX_') & contains(model.rxns, '(u)');
if ~any(isComEx)
    isComEx = startsWith(model.rxns, 'EX_') & contains(model.rxns, '[u]');
end
assert(any(isComEx), 'bench_make_community:noComEx', ...
    'No community [u] exchange reactions found after joining.');

idx  = find(isComEx);
open = rand(numel(idx), 1) < cfg.mediaOpenProb;

lo = log10(cfg.mediaMinInflow);
hi = log10(cfg.mediaMaxInflow);
inflow = 10.^(lo + (hi - lo)*rand(numel(idx), 1));

model.lb(idx)         = 0;                      % closed by default
model.lb(idx(open))   = -inflow(open);          % opened uptakes

info.nSpecies      = nSpecies;
info.rep           = rep;
info.species       = names(:);
info.nComEx        = numel(idx);
info.nOpen         = nnz(open);
info.minOpenInflow = min(inflow(open));
info.maxOpenInflow = max(inflow(open));
info.nRxns         = numel(model.rxns);
info.nMets         = numel(model.mets);
info.name          = sprintf('com%03d_rep%02d', nSpecies, rep);
info.time          = toc(t0);

fprintf(['  community %s: %d species, %d rxns, %d mets, %d/%d uptakes open, ' ...
         'inflow [%.2g %.2g], built in %.1f s\n'], ...
        info.name, nSpecies, info.nRxns, info.nMets, info.nOpen, info.nComEx, ...
        info.minOpenInflow, info.maxOpenInflow, info.time);
end

% -------------------------------------------------------------------------
function bm = agora_biomass_rxn(m)
% Name of a model's biomass reaction, or '' if it has none.
%
% Identify it by the OBJECTIVE VECTOR, not by its name. AGORA2 is not
% consistent about naming: a few models use biomass<NNN> (e.g.
% Abiotrophia_defectiva uses biomass087) but most use bio1, and every model
% also carries an EX_biomass(e) exchange that must never be mistaken for the
% reaction itself. Matching on the string 'biomass' therefore rejected 87% of
% the set. model.c is how the reconstruction declares which reaction is
% growth, so that is what we read.
bm = '';
if isfield(m,'c')
    oi = find(m.c ~= 0 & ~startsWith(m.rxns,'EX_'));
    if ~isempty(oi), bm = m.rxns{oi(1)}; return; end
end
% Fallback for a model shipped with a zeroed objective.
cand = m.rxns(~startsWith(m.rxns,'EX_') & ...
              (contains(lower(m.rxns),'biomass') | strncmpi(m.rxns,'bio',3)));
if ~isempty(cand), bm = cand{1}; end
end
