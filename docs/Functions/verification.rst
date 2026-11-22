Verification
============

An extraction promises that every core reaction is present and that the
result is flux consistent. The first always holds; the second does not
always hold for the LP formulations. ``checkExtraction`` reports on a result
so the difference is visible rather than silent.

Run it with the same arguments the extraction used, ``tol`` and
``useCoupling`` included: a check that answers a different question from
the extraction is not a check of it. Both of its methods are held to that
choice, which takes some doing on the ``fva`` side -- ``fluxVariability``
applies a model's coupling constraints whether or not you want it to, so
when ``useCoupling`` is false they are dropped from the copy it is given.

.. automodule:: src
   :members: checkExtraction
