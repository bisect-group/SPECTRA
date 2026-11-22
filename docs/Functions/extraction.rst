Model extraction
================

Builds a context-specific model around a set of core reactions.
``spectraCCME`` does the consistency check and the extraction in one pass,
for a model that is not known to be consistent.

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
   :members: spectraME, spectraCCME, forward
