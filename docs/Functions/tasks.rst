Metabolic tasks
===============

A task states something a network must be able to do: take up these
metabolites, produce those, carry flux through this reaction.

``essentialReactionsForTasks`` is what makes tasks usable by the extraction.
Requiring a model to perform a task is difficult to express as a constraint,
but requiring it to contain a reaction is what ``spectraME`` already does, so
the reactions without which a task becomes infeasible become its core set.

.. automodule:: src
   :members: metabolicTask, taskEquation, taskIsForced, taskConstraints,
             missingMetabolites, checkTask, checkTasks,
             essentialReactionsForTask, essentialReactionsForTasks,
             parseTaskList
