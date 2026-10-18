% Exports every live script in this directory to HTML, for the docs build.
%
% USAGE:
%
%    run mlx2html.m
%
% NOTE:
%    The .mlx files are the source of truth; edit those in MATLAB and re-run
%    this to refresh the HTML alongside them.

here = fileparts(mfilename('fullpath'));
scripts = dir(fullfile(here,'*.mlx'));

for k = 1:numel(scripts)
    name = scripts(k).name;
    fprintf('exporting %s\n', name);
    export(fullfile(here,name), Format='html');
end

fprintf('%d live script(s) exported\n', numel(scripts));
