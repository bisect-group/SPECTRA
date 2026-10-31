function bench_addpaths(cfg)
% bench_addpaths  Put SPECTRA, FASTERCC, scFASTCORMICS and COBRA on the path
% and initialise COBRA.
%
% Paths come from environment variables (see env.sh), so nothing here is
% machine specific. Adding a path is not a modification: all four codebases
% stay exactly as checked out.
%
% PATH PRECEDENCE MATTERS HERE. rFASTCORMICS ships several folders holding
% files with identical names, and COBRA is added with genpath, which PREPENDS.
% If the wrong copy wins, the benchmark silently measures a different
% implementation and still produces plausible timings. The checks at the end
% of this function are the safeguard, and they run LAST, after every addpath,
% so nothing added later can slip past them.

if nargin < 1 || isempty(cfg), cfg = bench_config(); end

addpath(cfg.root);

% Gurobi's MATLAB interface. Both methods solve through it now: SPECTRA_CC via
% solveCobraLP as always, FASTERCC via the LP backend swap in
% fastercc_lp_backend/. See bench_config for why CPLEX was abandoned.
gurobiMatlab = getenv('GUROBI_MATLAB_DIR');
if ~isempty(gurobiMatlab) && exist(gurobiMatlab, 'dir')
    addpath(gurobiMatlab);
end

% The IBM CPLEX MATLAB API. Only added if explicitly requested, so a stray
% cplexlp on the path cannot quietly pull FASTERCC back onto the unstable MEX.
if ~isempty(getenv('CPLEX_MATLAB_DIR')) && ...
        exist(getenv('CPLEX_MATLAB_DIR'), 'dir') && cfg.requireCplexlp
    addpath(getenv('CPLEX_MATLAB_DIR'));
end

% genpath for SPECTRA/scFASTCORMICS (they have subfolders); FASTERCC is flat.
for p = {cfg.spectraSrc, cfg.scfastSrc}
    if ~isempty(p{1}) && exist(p{1}, 'dir'), addpath(genpath(p{1})); end
end

if ~isempty(cfg.fasterccSrc) && exist(cfg.fasterccSrc, 'dir')
    addpath(cfg.fasterccSrc);
    rRoot = fileparts(cfg.fasterccSrc);

    % build_expanded_input_model calls fastcc_fastcore, in a sibling folder.
    % Safe to prepend: every file there is *_fastcore suffixed.
    addIfPresent(fullfile(rRoot, 'FASTCORE', 'Scripts'), false);

    % simplifyDuplicatedGenesFastbox needs GPRrulesMapper_rFASTCORMICS, and
    % call_build_expanded_input_model needs findEX_Rxns. Both live only in the
    % RNA-seq tree, which ALSO ships its own LP7/LP3/LP9_rFASTCORMICS.m and
    % fixIrr_rFASTCORMICS.m - the same names fasterCC calls. LP7/LP3/LP9 are
    % byte-identical to the FASTERCC copies, but fixIrr_rFASTCORMICS is NOT:
    % that file declares "function[model] = fixIrr(model)" and is a different
    % implementation. So both folders go on with '-end', behind FASTERCC.
    addIfPresent(fullfile(rRoot, 'rFASTCORMICS for RNA-seq data', 'Scripts'), true);
    addIfPresent(fullfile(rRoot, 'rFASTCORMICS for RNA-seq data', ...
                          'rFASTCORMICS Website Version 2019', 'Scripts'), true);
end

%% COBRA ------------------------------------------------------------------
if ~isempty(cfg.cobraSrc) && exist(cfg.cobraSrc, 'dir')
    addpath(genpath(cfg.cobraSrc));
    if isempty(getenv('COBRA_INITIALISED'))
        initCobraSafely();
    end
end

%% Precedence checks - LAST, so they see the final path -------------------
if ~isempty(cfg.fasterccSrc) && exist(cfg.fasterccSrc, 'dir')
    for f = {'fasterCC','fixIrr_rFASTCORMICS','LP7_rFASTCORMICS', ...
             'LP3_rFASTCORMICS','findDeadEndsFastbox','structureAnalyseFastbox'}
        assertResolvesIn(f{1}, cfg.fasterccSrc);
    end
end
if ~isempty(cfg.spectraSrc) && exist(cfg.spectraSrc, 'dir')
    for f = {'spectraCC','forwardcc','reverse'}
        assertResolvesIn(f{1}, cfg.spectraSrc);
    end
end
end

% -------------------------------------------------------------------------
function addIfPresent(d, atEnd)
if exist(d, 'dir')
    if atEnd, addpath(d, '-end'); else, addpath(d); end
end
end

% -------------------------------------------------------------------------
function assertResolvesIn(fname, wantDir)
% Guard against a shadowed function silently changing what is benchmarked.
p = which(fname);
assert(~isempty(p), 'bench_addpaths:missing', '%s not found on the path.', fname);

% Canonicalise both sides: on fir, /home/pavan/scratch is a symlink to
% /scratch/pavan, so the same file has two spellings and a plain string
% compare gives a false positive.
gotDir = canonical(fileparts(p));
wantC  = canonical(wantDir);
if ~strcmp(gotDir, wantC)
    error('bench_addpaths:shadowed', ...
        ['%s resolves to\n    %s\nbut must come from\n    %s\n' ...
         'A shadowing copy is ahead on the path; the benchmark would ' ...
         'measure a different implementation.'], fname, gotDir, wantC);
end
end

% -------------------------------------------------------------------------
function c = canonical(p)
% Resolve symlinks and relative segments. Falls back to the input if the JVM
% is unavailable (matlab -nojvm).
try
    c = char(java.io.File(p).getCanonicalPath());
catch
    c = p;
end
end

% -------------------------------------------------------------------------
function initCobraSafely()
% Work around a COBRA/Lmod interaction that breaks initCobraToolbox on fir.
%
% Lmod is wired into every subshell MATLAB spawns and invokes flatpak, which
% cannot load libXau.so.6 under MATLAB's own LD_LIBRARY_PATH. So inside
% MATLAB, system('git --version') returns 132 characters, not 18:
%
%   flatpak: error while loading shared libraries: libXau.so.6: ELF load
%   command address/offset not properly aligned
%   git version 2.41.0
%
% COBRA's checkGit() locates 'git version' with strfind (position
% independent) but then slices from a HARDCODED offset of 12, so it reads
% "ror whil" out of "error while loading", str2num returns [], and
% initCobraToolbox dies at
%
%   if installedGit && versionGit > 2100     % [] > 2100 is empty -> error
%
% Clearing LD_LIBRARY_PATH across init alone silences flatpak and lets the
% version parse. It is restored immediately afterwards so the CPLEX MEX still
% loads - solvers are selected AFTER the restore, never during init.
%
% Verified on fir 2026-10-23 (MATLAB 2023b.2, CPLEX 12.10):
%   initCobraToolbox OK | ibm_cplex LP set: 1 | solveCobraLP stat=1
% No edit to the COBRA checkout is required.

% initCobraToolbox(updateToolbox, mode): mode 'agent' sets agentMode, which
% skips EVERY git operation inside it. That fixes two separate failures:
%   - submodule init, which collides when several SLURM jobs share one COBRA
%     checkout ("The submodules could not be initialized")
%   - checkGit, whose fixed-offset version parse breaks here (see below)
% Fall back to the LD_LIBRARY_PATH workaround on a COBRA too old for the arg.
try
    initCobraToolbox(false, 'agent');
catch ME
    if isempty(strfind(ME.message, 'Too many input arguments')) %#ok<STREMP>
        rethrow(ME);
    end
    old = getenv('LD_LIBRARY_PATH');
    cleanup = onCleanup(@() setenv('LD_LIBRARY_PATH', old));
    setenv('LD_LIBRARY_PATH', '');
    initCobraToolbox(false);
    clear cleanup   % restore before any solver is selected
end
setenv('COBRA_INITIALISED', '1');
end
