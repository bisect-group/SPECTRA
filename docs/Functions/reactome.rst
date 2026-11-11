Minimal reactomes
=================

Cuts a model down to the smallest set of reactions that still meets a stated
requirement: growth in a medium, the production of a metabolite, a list of
metabolic tasks, or any combination of them.

This is the original question of Burgard, Vaidyaraman and Maranas (2001),
who asked how few reactions *Escherichia coli* needs to grow and found the
answer depends entirely on what is in the medium: 224 reactions on glucose
alone, 229 on acetate alone, 122 on a medium engineered to supply as much as
possible. A minimal reactome is minimal *for an environment*, never in the
abstract.

That is a minimum-cardinality problem, so it is a MILP. Every reaction gets
a binary :math:`y_j` tied to its flux by
:math:`lb_j y_j \le v_j \le ub_j y_j`, so that :math:`y_j = 0` pins the flux
to zero, and :math:`\sum_j w_j y_j` is minimised subject to :math:`Sv = 0`
and whatever the network has to achieve. That is ``minNetMILP`` with no core
reactions at all: the requirements enter as bounds rather than as a core set.

Growth, production and the medium are one condition, solved as one flux
vector. Tasks are not, because a task defines its own medium. The conditions
are therefore replicated -- one copy of the network each -- and tied
together by one binary per reaction, so that what is minimised is the size
of the union while every condition is satisfied separately:

.. math::

   \min_y \sum_j w_j y_j \quad\text{s.t.}\quad
   S^{(k)} v^{(k)} = 0, \;
   lb^{(k)}_j y_j \le v^{(k)}_j \le ub^{(k)}_j y_j
   \quad \forall k

This is the "protected functions" idea of NetworkReducer and MinNW: a
reduced network must keep not one phenotype but a list of them. Solving the
conditions jointly is what makes the answer minimal *across* them. Filling
them one at a time and taking the union, as ``gapfillForTasks`` does, gives
a feasible network but not the smallest one.

It is the opposite operation to :doc:`gapfilling`, and the two use the same
vocabulary of tasks: gap-filling adds reactions until the tasks pass, this
removes reactions while they still do.

.. automodule:: src
   :members: minimalReactome
