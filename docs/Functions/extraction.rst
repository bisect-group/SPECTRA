Model extraction
================

Builds a context-specific model around a set of core reactions.
``spectraCCME`` does the consistency check and the extraction in one pass,
for a model that is not known to be consistent, while ``spectraME2`` applies
the pure thermodynamic cone variant (ignoring finite box bounds and solving
``minNetLP``).

During the iterative direction-discovery phase, settled core reactions are
pinned directly at ``+/-tol`` (with zero convex combinations), ``forward``
(or ``forwardcc``) and ``reverse`` share a single random weight draw per
round, and Stage-1 LPs are skipped automatically when ``coreRxns`` is empty
or purely one-way. Reverse-irreversible reactions (``ub <= 0``) are oriented
forward only on an internal working copy -- the returned model is always
subsetted from the input model so original COBRA reaction stoichiometries
and bounds are preserved. ``altSolMethod`` defaults to ``'coreDirection'``
for LP formulations (``'minNetLP'``, ``'optimBiomass'``) and
``'pathwayExclusion'`` for MILP formulations (``'minNetMILP'``,
``'tradeOff'``).

Both take ``useCoupling`` as their last argument, and both default it to
false, so the problem is built from ``S``, ``lb`` and ``ub`` alone. Turn it
on for a model whose coupling constraints are part of what you are asking
about -- a community, where a unit that cannot grow must not be able to run
its metabolism either -- and whenever the extracted model is going straight
into FBA, which applies them regardless. See :doc:`consistency` for the
argument, and :doc:`formulations` for where the constraints enter the
problem.

``indicatorGroups`` is the other argument worth knowing: it says which
reactions share a binary in the MILP, named either by index or by reaction
ID. One binary per thing being selected over, rather than one per reaction,
is what makes a selection problem tractable on a large model, and the third
output reports what those binaries decided. Reading membership off the
returned model is not equivalent -- a reaction whose binary is off can
still carry a trace of flux.

.. automodule:: src
   :members: spectraME, spectraME2, spectraCCME, forward
