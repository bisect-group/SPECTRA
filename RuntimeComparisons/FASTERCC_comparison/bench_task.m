function bench_task(taskId, stage)
% bench_task  Array-job entry point: run one instance by its manifest index.
%
%   bench_task(7)
%   bench_task(7, 'pilot')
%
% Called from the submit script as
%   matlab -batch "bench_task($SLURM_ARRAY_TASK_ID, 'pilot')"

if nargin < 2 || isempty(stage), stage = 'pilot'; end
if ischar(taskId) || isstring(taskId), taskId = str2double(taskId); end

cfg = bench_config(stage);
bench_addpaths(cfg);

M = load(fullfile(cfg.modelDir, sprintf('manifest_%s.mat', cfg.stage)));
manifest = M.manifest;
assert(taskId >= 1 && taskId <= height(manifest), ...
    'Task id %d outside manifest range 1..%d', taskId, height(manifest));

r = manifest(taskId, :);
fprintf('task %d -> %s (~%d rxns)\n', taskId, r.name{1}, r.nRxns);
bench_run_one(r.clusters, r.nRemoved, r.rep, cfg);
end
