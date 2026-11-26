Gap-filling
===========

Adds the cheapest set of reactions from a universal model that lets a draft
meet a requirement: growth in each of a panel of media, or a list of
metabolic tasks.

What a reaction costs is what steers the result. The draft's own reactions
cost nothing, so the fill is free to keep all of them, and everything else
costs one, so it adds as little as it can. The ``weights`` argument sets
particular reactions against that default, which is where sequence evidence
goes: a reaction the organism's genome supports should cost less than one
it does not, so the fill prefers it. ``gapfillForTasks`` accepts all four
network inference formulations (``'minNetMILP'``, ``'minNetLP'``,
``'tradeOff'``, and ``'optimBiomass'``) for its core extraction step,
converting signed ``'tradeOff'`` evidence weights monotonically into
non-negative penalties for its task-by-task ``'minNetMILP'`` repair pass.

``gapfillCommunity`` does the same for every unit of a community at once,
with the biomass coupling applied -- a database reaction that can only
carry flux by ignoring the coupling is not a candidate.
Filling them together rather than one at a time lets a gap in one be closed
by another's secretion instead of by a new reaction, which is why the joint
result is the smaller one -- often by the whole of what a unit would
otherwise have had to invent.

.. automodule:: src
   :members: gapfillForGrowth, gapfillForTasks, gapfillCommunity
