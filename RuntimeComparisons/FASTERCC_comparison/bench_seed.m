function s = bench_seed(k, nRemove, rep, cfg)
% bench_seed  Deterministic seed for one benchmark condition.
%
% The whole test set is reproducible from cfg.baseSeed alone, which is what
% lets instances be regenerated on demand instead of stored. Any change to
% this function invalidates every cached result.

if nargin < 4 || isempty(cfg), cfg = bench_config(); end

% Index into the CANONICAL removal list, never the stage's subset: otherwise
% k07_rm500_rep01 would get one seed under 'pilot' (where 500 is the 2nd
% level) and a different one under 'full' (where it is the 6th), and the same
% named instance would not be the same model across stages.
r = find(cfg.allRemovalCounts == nRemove, 1);
assert(~isempty(r), 'bench_seed:removal', ...
       '%d is not one of cfg.allRemovalCounts', nRemove);
s = cfg.baseSeed + 1000000*k + 1000*r + rep;
end
