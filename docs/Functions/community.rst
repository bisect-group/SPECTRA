Communities and bodies
======================

Joins models into one network whose units trade metabolites through a pool.
The units may be organisms of a microbial community or tissues of a body;
nothing in the construction distinguishes them.

Each unit keeps its own copy of the network, tagged so that its reactions
and metabolites cannot be confused with another's, and the units meet in a
**pool**: a compartment several of them can reach. One unit's secretion
arrives in the pool and is available to the rest, which is what makes
cross-feeding possible. A pool named in ``environment`` also carries
exchange reactions, which is how anything enters or leaves the whole system.

There are two ways to build the pool, and both are in published use:

``'pooled'``
   Each unit keeps a private external compartment, and a transport reaction
   moves each exchangeable metabolite between it and the pool. Because the
   transport belongs to one unit, it can be bounded per unit -- which is
   where measured uptake and secretion rates go. Several pools can coexist,
   so the units need not all meet in the same place.

``'shared'``
   The units share the external compartment outright: one row of the
   stoichiometric matrix per exchangeable metabolite, for the whole
   community. There are no transport reactions and there is only ever one
   pool. For a hundred units across a few thousand metabolites this saves a
   hundred thousand reactions, which is the difference between a model that
   can be built and one that cannot.

The compartmentalised layout, the biomass coupling and its parameters follow
the joint FBA formulation used for minimal microbiomes in Raghu et al.
(2024), *npj Systems Biology and Applications* **10**, 46.

.. automodule:: src
   :members: buildCommunityModel, buildMultiTissueModel
