Consistency checking
====================

Removes the blocked reactions of a model, leaving one in which every
reaction can carry flux.

In each iteration, ``spectraCC`` draws a single random weight vector
:math:`w \sim U(1, 1.1)` over the still-unverified reactions and shares it
between ``forwardcc`` and ``reverse`` (after internally orienting
reverse-only irreversible reactions ``ub <= 0`` forward). Sharing the weight
vector cancels anti-parallel reversible terms across the forward/reverse pair
and guarantees almost-sure progress in every active round. Each LP is solved
at ``solverTol = 1e-9`` by default via :func:`~src.lpTolerances` without
mutating global ``CBT_LP_PARAMS``; pass ``solverTol = 'default'`` for
severely ill-scaled forced-flux models such as ``Harvey`` that require the
solver's default ``1e-6`` tolerance.

Flux consistency here is defined over the stoichiometry alone -- ``S``,
``lb`` and ``ub`` -- which is what fastcc, FASTCORE, FASTERCC and swiftcore
all answer, and what a comparison against any of them has to answer too. A
model may also carry coupling constraints, the ``C``, ``d`` and ``dsense``
fields that ``coupleRxnList2Rxn`` writes and ``buildCommunityModel`` uses
to tie every reaction of a unit to that unit's anchor. Those are left out
by default: a reaction blocked only by one is not blocked in the network,
it is blocked in that parameterisation.

Pass ``useCoupling`` to bring them in, and note which answer agrees with
the rest of COBRA -- ``optimizeCbModel`` applies ``C`` whether asked to or
not. On a community whose consumer cannot grow, FBA puts a transport's
range at ``[-0.01, 0]`` while a check ignoring ``C`` certifies it able to
carry 0.1. If the result is going to be built into a model and then
optimised, that disagreement is the one to avoid.

.. automodule:: src
   :members: spectraCC, forwardcc, reverse
