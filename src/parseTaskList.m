function [tasks,unresolved] = parseTaskList(path,model,skipUnresolved,freeInputs,freeOutputs)
% Reads a published metabolic task list and binds it to a model
%
% USAGE:
%
%    [tasks,unresolved] = parseTaskList(path,model,skipUnresolved,freeInputs,freeOutputs)
%
% INPUTS:
%    path:               Tab-separated task file. The columns follow the
%                        convention of the published lists: ID,
%                        DESCRIPTION, SHOULD FAIL, IN, IN LB, IN UB, OUT,
%                        OUT LB, OUT UB, EQU, EQU LB, EQU UB. Metabolites
%                        are written as name[compartment] and several are
%                        separated by semicolons.
%    model:              COBRA model structure the names are resolved
%                        against, by metabolite name and compartment rather
%                        than by identifier.
%
% OPTIONAL INPUTS:
%    skipUnresolved:     Whether to drop tasks naming metabolites the model
%                        does not have (Default: true). They can never be
%                        performed, so keeping them only produces failures
%                        no gap-fill could repair.
%    freeInputs:         References, written as name[compartment], that
%                        every task may take up whether or not it lists
%                        them (Default: none). Published lists are often
%                        written against a convention where something
%                        common may always be taken up, most often oxygen,
%                        and omit it.
%    freeOutputs:        References that every task may produce whether or
%                        not it lists them (Default: none). The same
%                        omission happens on the output side with common
%                        waste.
%
% OUTPUTS:
%    tasks:              Task structures from metabolicTask, as a struct
%                        array.
%    unresolved:         Struct array with fields taskID and references,
%                        listing what could not be resolved against the
%                        model.
%
% EXAMPLE:
%
%    [tasks, unresolved] = parseTaskList('metabolicTasks_Essential.txt', model);
%    % loosen every task by letting oxygen in and water out:
%    [tasks, ~] = parseTaskList('tasks.txt', model, true, {'O2[s]'}, {'H2O[s]'});
%
% NOTE:
%    References are resolved by metabolite *name* and compartment, not by
%    identifier, because the published lists are written against no
%    particular namespace. Where several metabolites share a name and
%    compartment the first is taken.
%
%    freeInputs and freeOutputs both loosen every task, and a loose enough
%    task passes for the wrong reason. Check them against the list's SHOULD
%    FAIL controls before trusting the result.
%
%    A metabolite written in both IN and OUT has two readings of its
%    bounds. The rule used here is the published format's: the output's
%    upper bound always wins, the output's lower bound wins only when it is
%    positive, and the input's lower bound is discarded either way. So an
%    output that merely permits production leaves the input's allowance to
%    consume intact, and the metabolite may still flow both ways; an output
%    that requires production replaces it. Requiring both at once is
%    contradictory and raises an error.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('skipUnresolved', 'var') || isempty(skipUnresolved)
    skipUnresolved = true;
end
if ~exist('freeInputs', 'var') || isempty(freeInputs)
    freeInputs = {};
end
if ~exist('freeOutputs', 'var') || isempty(freeOutputs)
    freeOutputs = {};
end

DEFAULT_BOUNDS = [0 1000];
DEFAULT_EQU_BOUNDS = [-1000 1000];

index = metaboliteIndex(model);
rows = readRows(path);

tasks = [];
unresolved = struct('taskID',{},'references',{});
seen = containers.Map('KeyType','char','ValueType','double');

% A list can hold hundreds of tasks, many of which are permissive by
% design. metabolicTask warns about each one; count them instead.
warnState = warning('off','metabolicTask:notForced');
cleanup = onCleanup(@() warning(warnState));
nUnforced = 0;

for r = 1:numel(rows)
    row = rows(r);
    identifier = strtrim(field(row,'ID'));
    description = strtrim(field(row,'DESCRIPTION'));
    if isempty(identifier) || startsWith(identifier,'#') || isempty(description)
        continue
    end

    % The published lists reuse a category as the ID, so make it unique.
    if isKey(seen, identifier)
        seen(identifier) = seen(identifier) + 1;
    else
        seen(identifier) = 1;
    end
    taskID = sprintf('%s_%03d', identifier, seen(identifier));

    missing = {};
    [inMets,inBounds,missing] = readSide(row,'IN',index,DEFAULT_BOUNDS,missing);
    [outMets,outBounds,missing] = readSide(row,'OUT',index,DEFAULT_BOUNDS,missing);

    [inMets,inBounds] = reconcileSides(taskID,inMets,inBounds,outMets,outBounds);

    [inMets,inBounds,missing] = addFree(freeInputs,index,inMets,inBounds, ...
        DEFAULT_BOUNDS,missing);
    [outMets,outBounds,missing] = addFree(freeOutputs,index,outMets,outBounds, ...
        DEFAULT_BOUNDS,missing);

    equations = [];
    equBounds = readBounds(row,'EQU',DEFAULT_EQU_BOUNDS);
    rawEquations = strsplit(field(row,'EQU'), ';');
    for e = 1:numel(rawEquations)
        raw = strtrim(rawEquations{e});
        if isempty(raw)
            continue
        end
        [mets,coeffs,bad] = parseEquation(raw,index);
        missing = [missing; bad]; %#ok<AGROW>
        if ~isempty(mets)
            equations = [equations; taskEquation(mets,coeffs,equBounds)]; %#ok<AGROW>
        end
    end

    if ~isempty(missing)
        unresolved(end+1).taskID = taskID; %#ok<AGROW>
        unresolved(end).references = unique(missing,'stable');
        if skipUnresolved
            continue
        end
    end

    task = metabolicTask(taskID, ...
        'description', description, ...
        'inputs', inMets, 'inputBounds', inBounds, ...
        'outputs', outMets, 'outputBounds', outBounds, ...
        'equations', equations, ...
        'shouldFail', ~isempty(strtrim(field(row,'SHOULD FAIL'))));
    if ~taskIsForced(task) && ~task.shouldFail
        nUnforced = nUnforced + 1;
    end
    tasks = [tasks; task]; %#ok<AGROW>
end

fprintf(['parseTaskList: %d task(s) parsed; %d named metabolites the model ' ...
         'does not have\n'], numel(tasks), numel(unresolved));
if nUnforced > 0
    warning('parseTaskList:unforcedTasks', ...
        ['%d parsed task(s) have no positive lower bound anywhere, so the ' ...
         'all-zero flux vector satisfies them and they pass trivially.'], ...
        nUnforced);
end
end


function index = metaboliteIndex(model)
% Indexes a model's metabolites by name and compartment
%
% USAGE:
%
%    index = metaboliteIndex(model)
%
% INPUTS:
%    model:              COBRA model structure.
%
% OUTPUTS:
%    index:              containers.Map from 'name|compartment' to the
%                        metabolite ID, with a lowercased key alongside
%                        each, since published lists are not always careful
%                        about case. Where several share both, the first is
%                        kept.

index = containers.Map('KeyType','char','ValueType','char');
if isfield(model,'metNames')
    names = model.metNames;
else
    names = model.mets;
end
for k = 1:numel(model.mets)
    comp = compartmentOf(model.mets{k});
    if isempty(comp)
        continue
    end
    name = strtrim(names{k});
    for key = {[name '|' comp], [lower(name) '|' comp]}
        if ~isKey(index, key{1})
            index(key{1}) = model.mets{k};
        end
    end
end
end


function comp = compartmentOf(metID)
% Returns the compartment symbol of a metabolite identifier
%
% USAGE:
%
%    comp = compartmentOf(metID)
%
% INPUTS:
%    metID:              A metabolite identifier, written either as
%                        'name[c]' or as 'name_c'.
%
% OUTPUTS:
%    comp:               The compartment symbol, or '' if there is none.

token = regexp(metID, '\[(\w+)\]\s*$', 'tokens', 'once');
if isempty(token)
    token = regexp(metID, '_(\w)\s*$', 'tokens', 'once');
end
if isempty(token)
    comp = '';
else
    comp = token{1};
end
end


function metID = resolveReference(reference,index)
% Returns the metabolite a 'name[compartment]' reference points to
%
% USAGE:
%
%    metID = resolveReference(reference,index)
%
% INPUTS:
%    reference:          The reference, as written in the task file.
%    index:              The map from metaboliteIndex.
%
% OUTPUTS:
%    metID:              The metabolite ID, or '' if the model lacks it.

metID = '';
token = regexp(strtrim(reference), '^(.+?)\s*\[\s*([^\]]+)\s*\]$', 'tokens', 'once');
if isempty(token)
    return
end
name = strtrim(token{1});
comp = strtrim(token{2});
for key = {[name '|' comp], [lower(name) '|' comp]}
    if isKey(index, key{1})
        metID = index(key{1});
        return
    end
end
end


function [mets,coeffs,unresolved] = parseEquation(equation,index)
% Turns a task equation into stoichiometry
%
% USAGE:
%
%    [mets,coeffs,unresolved] = parseEquation(equation,index)
%
% INPUTS:
%    equation:           The equation, such as 'ATP[c] + H2O[c] => ADP[c]'.
%    index:              The map from metaboliteIndex.
%
% OUTPUTS:
%    mets:               Metabolite IDs, as a cell array of char.
%    coeffs:             The coefficient of each, consuming on the left and
%                        producing on the right.
%    unresolved:         References the model does not have.
%
% NOTE:
%    Terms are split on ' + ' rather than on '+', because the published
%    lists write the proton as 'H+[c]' and a bare split shears it into 'H'
%    and '[c]'.

ARROWS = {'<==>','<=>','-->','=>','->'};

mets = {};
coeffs = [];
unresolved = {};

arrow = '';
for a = 1:numel(ARROWS)
    if contains(equation, ARROWS{a})
        arrow = ARROWS{a};
        break
    end
end
if isempty(arrow)
    unresolved = {equation};
    return
end

parts = strsplit(equation, arrow);
sides = {parts{1}, strjoin(parts(2:end), arrow)};
signs = [-1 1];

acc = containers.Map('KeyType','char','ValueType','double');
order = {};
for s = 1:2
    terms = regexp(sides{s}, '\s\+\s', 'split');
    for t = 1:numel(terms)
        term = strtrim(terms{t});
        if isempty(term)
            continue
        end
        [first,rest] = strtok(term);
        factor = str2double(first);
        if ~isnan(factor) && ~isempty(strtrim(rest))
            reference = strtrim(rest);
        else
            factor = 1;
            reference = term;
        end
        metID = resolveReference(reference,index);
        if isempty(metID)
            unresolved{end+1,1} = reference; %#ok<AGROW>
            continue
        end
        if isKey(acc,metID)
            acc(metID) = acc(metID) + signs(s)*factor;
        else
            acc(metID) = signs(s)*factor;
            order{end+1,1} = metID; %#ok<AGROW>
        end
    end
end

mets = order;
coeffs = zeros(numel(order),1);
for k = 1:numel(order)
    coeffs(k) = acc(order{k});
end
end


function bounds = readBounds(row,prefix,default)
% Reads a row's lower and upper bound for one column group
%
% USAGE:
%
%    bounds = readBounds(row,prefix,default)
%
% INPUTS:
%    row:                One parsed row of the task file.
%    prefix:             The column prefix, 'IN', 'OUT' or 'EQU'.
%    default:            What to use where the row is blank.
%
% OUTPUTS:
%    bounds:             [lower upper].

bounds = default;
rawLower = strtrim(field(row,[prefix ' LB']));
rawUpper = strtrim(field(row,[prefix ' UB']));
if ~isempty(rawLower) && ~isnan(str2double(rawLower))
    bounds(1) = str2double(rawLower);
end
if ~isempty(rawUpper) && ~isnan(str2double(rawUpper))
    bounds(2) = str2double(rawUpper);
end
end


function [mets,bounds,missing] = readSide(row,prefix,index,default,missing)
% Reads one side of a task row into metabolite IDs and bounds
%
% USAGE:
%
%    [mets,bounds,missing] = readSide(row,prefix,index,default,missing)
%
% INPUTS:
%    row:                One parsed row of the task file.
%    prefix:             'IN' or 'OUT'.
%    index:              The map from metaboliteIndex.
%    default:            Bounds to use where the row is blank.
%    missing:            References already found unresolvable.
%
% OUTPUTS:
%    mets:               Metabolite IDs on that side.
%    bounds:             An n-by-2 matrix of bounds, one row per metabolite.
%    missing:            The input list, with anything new appended.

sideBounds = readBounds(row,prefix,default);
references = strsplit(field(row,prefix), ';');

mets = {};
bounds = zeros(0,2);
for k = 1:numel(references)
    reference = strtrim(references{k});
    if isempty(reference)
        continue
    end
    metID = resolveReference(reference,index);
    if isempty(metID)
        missing{end+1,1} = reference; %#ok<AGROW>
        continue
    end
    hit = strcmp(mets, metID);
    if any(hit)
        bounds(hit,:) = sideBounds;
    else
        mets{end+1,1} = metID; %#ok<AGROW>
        bounds(end+1,:) = sideBounds; %#ok<AGROW>
    end
end
end


function [inMets,inBounds] = reconcileSides(taskID,inMets,inBounds,outMets,outBounds)
% Settles the bounds of a metabolite written on both sides of a task
%
% USAGE:
%
%    [inMets,inBounds] = reconcileSides(taskID,inMets,inBounds,outMets,outBounds)
%
% INPUTS:
%    taskID:             The task's identifier, for the error message.
%    inMets, inBounds:   The input side.
%    outMets, outBounds: The output side.
%
% OUTPUTS:
%    inMets, inBounds:   The input side, with overlaps settled.
%
% NOTE:
%    The output's lower bound wins when it is positive, in which case the
%    input is dropped entirely; otherwise the input keeps its allowance to
%    consume but loses any requirement. Requiring uptake and production of
%    the same metabolite at once is contradictory and is an error.

drop = false(numel(inMets),1);
for k = 1:numel(inMets)
    hit = strcmp(outMets, inMets{k});
    if ~any(hit)
        continue
    end
    inLower = inBounds(k,1);
    outLower = outBounds(hit,1);
    outLower = outLower(1);
    if inLower > 0 && outLower > 0
        error('parseTaskList:contradictorySides', ...
            ['Task ''%s'' requires ''%s'' to be both taken up at %g or ' ...
             'more and produced at %g or more, which cannot both hold.'], ...
            taskID, inMets{k}, inLower, outLower);
    end
    if outLower > 0
        drop(k) = true;
    else
        inBounds(k,1) = 0;
    end
end
inMets = inMets(~drop);
inBounds = inBounds(~drop,:);
end


function [mets,bounds,missing] = addFree(references,index,mets,bounds,default,missing)
% Adds the references every task is allowed, where it does not list them
%
% USAGE:
%
%    [mets,bounds,missing] = addFree(references,index,mets,bounds,default,missing)
%
% INPUTS:
%    references:         References to add, as name[compartment].
%    index:              The map from metaboliteIndex.
%    mets, bounds:       The side to add to.
%    default:            Bounds to give them.
%    missing:            References already found unresolvable.
%
% OUTPUTS:
%    mets, bounds:       The side, with the new entries appended.
%    missing:            The input list, with anything new appended.

for k = 1:numel(references)
    metID = resolveReference(references{k},index);
    if isempty(metID)
        missing{end+1,1} = references{k}; %#ok<AGROW>
        continue
    end
    if ~any(strcmp(mets, metID))
        mets{end+1,1} = metID; %#ok<AGROW>
        bounds(end+1,:) = default; %#ok<AGROW>
    end
end
end


function rows = readRows(path)
% Reads a tab-separated task file into a struct array of columns
%
% USAGE:
%
%    rows = readRows(path)
%
% INPUTS:
%    path:               The task file.
%
% OUTPUTS:
%    rows:               Struct array with fields headers and values, one
%                        entry per data row. Read through field() rather
%                        than directly, since the published headers contain
%                        spaces and are not valid field names.

fid = fopen(path,'r');
if fid < 0
    error('parseTaskList:cannotOpen','Cannot open ''%s''.', path);
end
cleanup = onCleanup(@() fclose(fid));

TAB = sprintf('\t');
headerLine = fgetl(fid);
if ~ischar(headerLine)
    error('parseTaskList:emptyFile','''%s'' is empty.', path);
end
headers = strtrim(strsplit(headerLine, TAB, 'CollapseDelimiters', false));

rows = struct('headers',{},'values',{});
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    if isempty(strtrim(line))
        continue
    end
    values = strsplit(line, TAB, 'CollapseDelimiters', false);
    if numel(values) < numel(headers)
        values(end+1:numel(headers)) = {''};
    end
    rows(end+1).headers = headers; %#ok<AGROW>
    rows(end).values = values(1:numel(headers));
end
end


function value = field(row,name)
% Returns one column of a row by its header, or '' if absent
%
% USAGE:
%
%    value = field(row,name)
%
% INPUTS:
%    row:                One entry of the struct array from readRows.
%    name:               The column header.
%
% OUTPUTS:
%    value:              The cell's contents, as char.

hit = strcmpi(row.headers, name);
if any(hit)
    value = row.values{find(hit,1)};
else
    value = '';
end
end
