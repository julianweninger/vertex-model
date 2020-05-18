
``PCPTopology`` - Perform topological changes on a cell tissue
==============================================================


Model Fundamentals
------------------

This model implements topological changes and other operations for a cell tissue that is described in the PCPVertex model.

The idea is to initiate a topological change, e.g. a cell division, by changing parameters of the PCPVertex model or its objects or to call operations on the  PCPVertex model at specific times.
These specific times are specified in ``OperationBundle``, please check its documentation.

.. note::
    The currently implemented operations are available in the collections ``OperationCollection``.
    Please check their documentation for the configuration of the respective operation.

The energy is minimized by the vertex model according to ``OperationBundle::MinimizationMode`` after every or all iterations.
The parameters for minimization can be updated for every operation wrt to the defaults of ``PCPTopology``.

The order of application of the operations is maintained throughout the simulation and is according to the order of registration.
Note, that not every operation is called in every step.

Default Model Configuration
---------------------------

Below are the default configuration parameters for the ``PCPTopology`` model.

.. literalinclude:: ../../src/models/PCPTopology/PCPTopology_cfg.yml
   :language: yaml
   :start-after: ---


References
----------

* Aigouy, B. et al. Cell Flow Reorients the Axis of Planar Polarity in the Wing Epithelium of Drosophila. 14 (2010) doi:10.1016/j.cell.2010.07.042.
* Etournay, R. et al. Interplay of cell dynamics and epithelial tension during morphogenesis of the Drosophila pupal wing. eLife 4, e07090 (2015).
* Farhadifar, R., Röper, J.-C., Aigouy, B., Eaton, S. & Jülicher, F. The Influence of Cell Mechanics, Cell-Cell Interactions, and Proliferation on Epithelial Packing. Current Biology 17, 2095–2104 (2007).