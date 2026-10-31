function n = bench_freeze_agora_pool(cfg)
% bench_freeze_agora_pool  Fix the AGORA2 sampling universe, once.
%
%   n = bench_freeze_agora_pool(cfg)
%
% Writes agora2/pool.txt: the sorted list of model files the community
% experiment may sample from. bench_make_community reads only this file,
% never the directory, so an instance means the same thing whenever it is
% rebuilt - before, during or after any further downloading.
%
% Run this ONCE, after the download settles. Re-running it after more models
% arrive changes the universe and therefore changes every community instance,
% so it refuses to overwrite an existing pool unless explicitly forced.
%
% Each file is opened and checked, not just listed. Drive serves a small HTML
% error page when it rate-limits, and a truncated or non-model .mat left in
% the folder would otherwise enter the pool and fail much later, mid-run,
% after hours of queueing.

if nargin < 1 || isempty(cfg), cfg = bench_config('smoke'); end

poolFile = fullfile(fileparts(cfg.agoraDir), 'pool.txt');
if exist(poolFile, 'file')
    error('bench_freeze_agora_pool:exists', ...
        ['%s already exists. The pool defines which community every seed ' ...
         'builds, so overwriting it invalidates every result produced so ' ...
         'far. Delete it deliberately if that is what you want.'], poolFile);
end

d = dir(fullfile(cfg.agoraDir, '*.mat'));
assert(~isempty(d), 'bench_freeze_agora_pool:empty', ...
    'No .mat files in %s', cfg.agoraDir);

keep = false(numel(d), 1);
why  = strings(numel(d), 1);
for i = 1:numel(d)
    f = fullfile(cfg.agoraDir, d(i).name);
    try
        S = load(f);
        fn = fieldnames(S);
        m  = S.(fn{1});
        hasFields = all(isfield(m, {'S','rxns','mets','lb','ub'}));
        hasBm     = ~isempty(agora_biomass_rxn(m));
        if hasFields && hasBm
            keep(i) = true;
        else
            why(i) = "missing fields or biomass reaction";
        end
    catch ME
        why(i) = string(ME.message);
    end
    if mod(i, 250) == 0
        fprintf('  checked %d/%d, %d usable\n', i, numel(d), nnz(keep));
    end
end

names = sort({d(keep).name}');
fid = fopen(poolFile, 'w');
fprintf(fid, '%s\n', names{:});
fclose(fid);

n = numel(names);
fprintf('\nPool frozen: %d usable of %d files -> %s\n', n, numel(d), poolFile);
bad = find(~keep);
if ~isempty(bad)
    fprintf('Rejected %d:\n', numel(bad));
    for i = bad(1:min(10,numel(bad)))'
        fprintf('  %-55s %s\n', d(i).name, why(i));
    end
end
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
