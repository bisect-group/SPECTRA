Network inference formulations
==============================

What the extraction optimises. ``minNetLP`` minimises weighted total
absolute flux, ``minNetMILP`` the weighted reaction count itself,
``tradeOff`` maximises net weight so that evidence can push a reaction in or
out, and ``growthOptim`` rewards the objective reactions and charges for
everything else.

.. automodule:: src
   :members: minNet, tradeOff, growthOptim
