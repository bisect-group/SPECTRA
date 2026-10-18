function forced = taskIsForced(task)
% Return whether anything in a task requires nonzero flux.
%
% USAGE:
%
%    forced = taskIsForced(task)
%
% INPUTS:
%    task:               Task structure from metabolicTask.
%
% OUTPUTS:
%    forced:             true when some output or equation has a positive
%                        lower bound. A task for which this is false is
%                        satisfied by the all-zero flux vector, so it
%                        passes trivially.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

forcedOut = ~isempty(task.outputBounds) && any(task.outputBounds(:,1) > 0);

forcedEqu = false;
for k = 1:numel(task.equations)
    if task.equations(k).bounds(1) > 0
        forcedEqu = true;
        break
    end
end

forced = forcedOut || forcedEqu;
end
