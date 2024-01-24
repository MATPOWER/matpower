.. _sec_task:

Task Object
===========

The task object is the one that builds and manages the model objects in order to solve the problem of interest. The :class:`mp.task` base class implements much of the functionality, with PF, CPF and OPF subclasses, namely :class:`mp.task_pf`, :class:`mp.task_cpf`, and :class:`mp.task_opf`, respectively, specifying the model classes to use and implementing other problem-specific functionality. The typical usage pattern is simply to construct the task object for the problem of interest, then call its :meth:`run() <mp.task.run>` method, passing in a struct of input data, a |MATPOWER| options struct, and an optional cell array of |MATPOWER| extensions.

Running a Task
--------------

Most of the action related to the task object occurs in the :meth:`run() <mp.task.run>` method.
In a typical case, as illustrated in :numref:`code_task_run_eg`, it simply builds the objects for the three model layers sequentially, solves the math model, and uses the results to update the network and data models. However, the actual :meth:`run() <mp.task.run>` method also allows for each model layer to iterate with a modified instance of the model as shown in the flowchart in :numref:`fig_task_run`. This can be used, for example, to iteratively update and re-solve a power flow in order to automatically satisfy the generator reactive power limits.

.. _fig_task_run:
.. figure:: figures/task-run.*
   :alt: Flowchart of task run() method
   :align: center

   Flowchart of task :meth:`run() <mp.task.run>` method


Building Model and Converter Objects
------------------------------------

Each of the build steps, marked with the stars in :numref:`fig_task_run`, consists of the following sub-steps:

1. Determine the class for the corresponding container object. There is a default, defined by a task method, but it can be overridden by a task subclass, or modified by user options or extensions.
2. Construct the container object.
3. Determine the set of classes for the individual element objects. The container class defines the defaults, but they can also be modified by user options or extensions.
4. Call the container object's :meth:`build` method to construct the element objects and complete the build process.

Iterative Execution
-------------------

As mentioned above, the :meth:`run() <mp.task.run>` method allows for an iterative solution at any of the three modeling layers. This is accomplished by overridding the :meth:`next_dm() <mp.task.next_dm>`, :meth:`next_nm() <mp.task.next_nm>`, or :meth:`next_mm() <mp.task.next_mm>` methods, respectively, for the data, network or math models. By default, these methods return an empty matrix, indicating that iteration should terminate. On the other hand, if a modified model object is returned, it triggers a new iteration with the modified model.

This feature is used by both PF and CPF to implement enforcement of certain constraints, such as generator reactive power limits.

Other Methods
-------------

A task also has a :meth:`print_soln() <mp.task.print_soln>` method for pretty printing the solution to the console and a :meth:`save_soln() <mp.task.save_soln>` method for saving the saved case to a file.
