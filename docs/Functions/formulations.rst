Network inference formulations
==============================

What the extraction optimises. ``minNetLP`` minimises weighted total
absolute flux, ``minNetMILP`` the weighted reaction count itself,
``tradeOff`` maximises net weight so that evidence can push a reaction in or
out, and ``growthOptim`` rewards the objective reactions and charges for
everything else.

Each formulation builds its own solver problem out of ``model.S``, so a
constraint that is not in ``S`` is not in the problem unless it is put
there. ``couplingRows`` is what puts it there: it reads a model's coupling
constraints -- the ``C``, ``d`` and ``dsense`` fields that
``coupleRxnList2Rxn`` writes and ``optimizeCbModel`` honours -- and returns
them as rows of the problem, padded to its width. A community model ties
each unit's reactions to its own anchor through ``C``, so an extraction
that ignored it would be free to run a unit's metabolism while the unit is
not growing. A model with no such fields gives no rows and is unaffected.

The mixed-integer formulations hand their tolerances to the solve rather
than setting them globally, which ``milpTolerances`` assembles. The
integrality tolerance is the one that decides what the answer contains and
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
   :members: minNet, tradeOff, growthOptim, resolveIndicatorGroups, couplingRows, milpTolerances, assertCoreProgress
