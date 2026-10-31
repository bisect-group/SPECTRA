function [class, info] = spectraCCdir(model, tol, consistencyType, seed, detectionCutoff, opts)
% spectraCCdir  Directional consistency check in the style of fasterCC revopt=1.
%
%   [class, info] = spectraCCdir(model, tol, consistencyType, seed, detectionCutoff, opts)
%
% Stock spectraCC answers one question - can this reaction carry flux at all -
% because it marks a reaction consistent on abs(flux) >= cutoff and lets the
% forward and reverse passes feed the SAME remaining-core set. Once the
% forward LP proves a reaction, the reverse LP never tests it, so the
% direction information is destroyed.
%
% This variant keeps the two directions apart and classifies every reaction:
%
%   'forward'     carries flux only in the forward direction
%   'reverse'     carries flux only in the reverse direction
%   'reversible'  carries flux in both
%   'blocked'     proven in neither
%
% HOW IT DIFFERS FROM spectraCC, concretely:
%   1. Two independent remaining-core sets, coreF and coreR, instead of one.
%   2. SIGNED flux decides membership: flux >= +cutoff proves forward,
%      flux <= -cutoff proves reverse. Never abs().
%   3. The loop runs until BOTH sets stop shrinking, not just one.
%
% INPUTS / OPTIONAL INPUTS
%   model, tol, consistencyType, seed, detectionCutoff
%       As spectraCC. Defaults match it: tol 1e-4, 'stoichiometry',
%       detectionCutoff 0.99*tol, seed unset.
%   opts    struct, all optional:
%       .crossEvidence  (default TRUE) harvest the opposite direction from
%                       each LP as well. Every LP returns a complete flux
%                       vector satisfying S*v = 0 and the bounds, so if a
%                       forward LP happens to return v_j <= -cutoff, that IS
%                       a feasible steady-state distribution in which j
%                       carries reverse flux - it proves reverse feasibility
%                       just as firmly as the reverse LP would. Free evidence,
%                       strictly fewer LPs, identical classification.
%                       Set false to keep the passes fully independent.
%       .confirmBlocked (default false) after the batch loop, retest each
%                       still-unproven reaction one at a time, as fastcc and
%                       fasterCC do. See the caveat below.
%
% OUTPUTS
%   class   n x 1 categorical, one label per reaction of the INPUT model
%   info    struct with
%       .forward/.reverse/.reversible/.blocked   reaction index lists
%       .fwdFeasible/.revFeasible                n x 1 logicals
%       .LPS                                     LPs solved
%       .consistent      indices feasible in either direction. This is the
%                        set stock spectraCC returns, so it is directly
%                        comparable to spectraCC/fastcc/fasterCC output.
%       .flipped         reactions normalised from reverse-irreversible, whose
%                        labels were mapped back to the input orientation
%
% CAVEAT ON 'blocked'
%   spectraCC's loop stops as soon as one batch LP proves nothing new; unlike
%   fastcc it has no one-by-one (LP3) fallback. A batch LP maximises a random
%   objective, so a reaction can be feasible yet never get proven. 'blocked'
%   therefore means "not proven feasible", which can over-report blocked
%   relative to FVA. Set opts.confirmBlocked = true to retest leftovers
%   individually, and validate against FVA before publishing - which is
%   exactly what the FASTERCC paper did for its own revopt=1 reversibility
%   correction (its Table 3, FVA on Recon3D, iHuman1, iMM1415, iMM1865).
%   check_spectraCCdir_fva.m does that check.

if ~exist('tol','var') || isempty(tol), tol = 1e-4; end
if ~exist('consistencyType','var') || isempty(consistencyType)
    consistencyType = 'stoichiometry';
end
if ~exist('detectionCutoff','var') || isempty(detectionCutoff)
    detectionCutoff = 0.99*tol;
end
if ~exist('opts','var') || isempty(opts), opts = struct(); end
if ~isfield(opts,'crossEvidence'),  opts.crossEvidence  = true;  end
if ~isfield(opts,'confirmBlocked'), opts.confirmBlocked = false; end
if exist('seed','var') && ~isempty(seed), rng(seed); end

switch consistencyType
    case 'stoichiometry', steadystate = 1;
    case 'topology',      steadystate = 0;
    otherwise
        error('spectraCCdir:badConsistencyType', ...
            ['consistencyType must be ''stoichiometry'' or ''topology'', ' ...
             'not ''%s''.'], consistencyType);
end

n0 = numel(model.rxns);

%% Normalisation - identical to spectraCC -------------------------------
% Reactions irreversible in the REVERSE direction are flipped so that every
% irreversible reaction points forward. Their labels are flipped back at the
% end, so the result is always stated in the INPUT model's orientation.
[~, n] = size(model.S);
IrR  = model.ub <= 0;
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR)  = -1*model.lb(IrR);
model.lb(IrR)  = -1*temp;
rev = ones(n,1);
rev(model.lb >= 0) = 0;
model.rev = rev;

%% Two independent passes ------------------------------------------------
% coreF: not yet proven forward-feasible - every reaction is a candidate.
% coreR: not yet proven reverse-feasible - only reversible reactions can be,
%        since after normalisation an irreversible reaction has lb >= 0.
coreF = true(n,1);
coreR = logical(rev);

LPS = 0;
prevF = -1; prevR = -1;

while sum(coreF) ~= prevF || sum(coreR) ~= prevR
    prevF = sum(coreF);
    prevR = sum(coreR);

    % ---- forward pass: drives core reactions to v >= +tol
    if any(coreF)
        flux = forwardcc(model, coreF, tol, steadystate);
        LPS = LPS + 1;
        coreF(flux >=  detectionCutoff) = false;
        if opts.crossEvidence
            coreR(flux <= -detectionCutoff) = false;
        end
    end

    % ---- reverse pass: drives reversible core reactions to v <= -tol
    if any(coreR)
        flux = reverse(model, coreR, tol, steadystate);
        LPS = LPS + 1;
        coreR(flux <= -detectionCutoff) = false;
        if opts.crossEvidence
            coreF(flux >=  detectionCutoff) = false;
        end
    end
end

fwdFeasible = ~coreF;
revFeasible = logical(rev) & ~coreR;

%% Optional one-by-one confirmation of the leftovers ---------------------
if opts.confirmBlocked
    for j = find(~fwdFeasible)'
        c = false(n,1); c(j) = true;
        flux = forwardcc(model, c, tol, steadystate);
        LPS = LPS + 1;
        if flux(j) >= detectionCutoff, fwdFeasible(j) = true; end
        if opts.crossEvidence
            revFeasible(flux <= -detectionCutoff & logical(rev)) = true;
        end
    end
    for j = find(logical(rev) & ~revFeasible)'
        c = false(n,1); c(j) = true;
        flux = reverse(model, c, tol, steadystate);
        LPS = LPS + 1;
        if flux(j) <= -detectionCutoff, revFeasible(j) = true; end
        if opts.crossEvidence
            fwdFeasible(flux >= detectionCutoff) = true;
        end
    end
end

%% Map labels back to the input orientation ------------------------------
% For a flipped reaction, forward in the normalised model is reverse in the
% input model. Swap before classifying so 'forward' always means forward as
% the caller supplied it.
fwdOut = fwdFeasible; revOut = revFeasible;
fwdOut(IrR) = revFeasible(IrR);
revOut(IrR) = fwdFeasible(IrR);

class = repmat(categorical({'blocked'}), n0, 1);
class(fwdOut & ~revOut) = categorical({'forward'});
class(~fwdOut &  revOut) = categorical({'reverse'});
class(fwdOut &  revOut)  = categorical({'reversible'});

info.forward     = find(class == 'forward');
info.reverse     = find(class == 'reverse');
info.reversible  = find(class == 'reversible');
info.blocked     = find(class == 'blocked');
info.fwdFeasible = fwdOut;
info.revFeasible = revOut;
info.consistent  = find(fwdOut | revOut);   % comparable to spectraCC output
info.LPS         = LPS;
info.flipped     = find(IrR);
info.opts        = opts;
end
