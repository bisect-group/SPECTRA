Network inference formulations
==============================

What the extraction optimises. ``minNetLP`` minimises weighted total
absolute flux (using ``"With-Chain"`` median-normalised weights and
warm-started IRLS), ``minNetMILP`` the weighted reaction count itself,
``tradeOff`` maximises net weight so that evidence can push a reaction in or
out, and ``optimBiomass`` rewards the objective reactions and charges for
everything else.

Note that because ``tradeOff`` maximises :math:`\sum w_i z_i` subject only to
mass balance and bounds, any disconnected internal stoichiometrically
balanced cycle (Type III loop) with net positive weight can circulate flux of
at least ``tol`` independently of the core set and exchanges. In both
``minNetMILP`` and ``tradeOff``, the extracted subnetwork takes the union of
core reactions, switched-on binary indicators (``z >= 0.5``), and reactions
carrying flux at or above ``inclusionCutoff`` (``tol * 1e-7``), while
``pathwayExclusion`` cuts act strictly over the switched-on indicator sets.

Each formulation builds its own solver problem out of ``model.S``, so a
constraint that is not in ``S`` is not in the problem unless it is put
there. ``couplingRows`` is what puts it there: it reads a model's coupling
constraints -- the ``C``, ``d`` and ``dsense`` fields that
``coupleRxnList2Rxn`` writes -- and returns them as rows of the problem,
padded to its width.

It is **off by default**, and every routine here takes a ``useCoupling``
argument to turn it on. Flux consistency, and the extractions built on it,
are conventionally defined over the stoichiometry alone: that is what
fastcc, FASTCORE, FASTERCC and swiftcore all answer, and a comparison
against them has to answer the same question. A coupling constraint is a
modelling device rather than stoichiometry, so a reaction blocked only by
one is not blocked in the network -- it is blocked in that
parameterisation.

Switch it on when the parameterisation is the point. A community model ties
each unit's reactions to its own anchor through ``C``, and an extraction
ignoring that is free to run a unit's metabolism while the unit is not
growing, which is the whole thing the coupling forbids; ``gapfillCommunity``
and ``minimalMicrobiome`` therefore ask for it. Switch it on, too, whenever
the answer has to agree with ``optimizeCbModel`` on the same model, since
that honours ``C`` always. On a toy community whose consumer cannot grow,
FBA puts a transport's range at ``[-0.01, 0]`` while a check ignoring ``C``
reports it able to carry 0.1.

A model with no such fields gives no rows either way.

The linear and mixed-integer formulations hand their tolerances to each solve
rather than setting them globally, assembled by ``lpTolerances`` and
``milpTolerances``. Pass ``solverTol = 'default'`` to either helper (or any
caller) to return ``{}`` and leave the active COBRA session / solver
tolerances untouched. The integrality tolerance is the one that decides what
the MILP answer contains and
the one most easily overlooked: a binary is tied to its reaction by
:math:`lb_i y_i \le v_i \le ub_i y_i`, so a binary sitting at the
tolerance and still counted as zero permits a flux of
:math:`\text{tol} \cdot ub_i` -- at a solver's usual 1e-5 and the usual
bound of 1000, a flux of 0.01, far above ``inclusionCutoff``. Reactions the
solver switched *off* would then read as carrying flux and be kept. The
feasibility tolerance is the complementary failure: too loose, and a core
reaction is present but cannot carry the flux it was required to. Both are
silent, which is why the default is 1e-9, the floor the solvers accept,
rather than their own looser settings.

Sharing a binary between several reactions, through ``indicatorGroups``, is
how a requirement that holds in several conditions at once is expressed:
replicate the network once per condition, group each reaction's copies
together, and the objective counts reactions while the constraints are
satisfied separately in every condition. Naming one reaction per group is
the other use, and it is what makes a selection problem tractable on a
large model -- one binary per organism of a community, say, rather than one
per reaction. Every free reaction left out of the groups keeps no binary
and is therefore kept.

.. automodule:: src
   :members: minNet, tradeOff, optimBiomass, resolveIndicatorGroups, couplingRows, milpTolerances, lpTolerances, assertCoreProgress
