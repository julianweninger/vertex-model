
``PCPTopology`` - Perform topological changes on a cell tissue
==============================================================


Model Fundamentals
------------------

This model implements topological changes for a cell tissue that is described
in the PCPVertex model.

The idea is to initiate a topological change, e.g. a cell division, by changing
parameters of the PCPVertex model or its objects or to call operations on the 
PCPVertex model at specific times.
After every operations the PCPVertex model is iterated until relative change in
energy falls below a threshold value.

Currently implemented topological changes
    * Cell division

The topological changes
^^^^^^^^^^^^^^^^^^^^^^^

from Farhadifar et al [2007] and Aigouy et al. [2010]

**Cell division**
    A cell division begins with selecting a random cell from the oldest generation:
    Cells from initialization are called generation 0, upon cell division the parent cell of generation i is divided into 2 daughter cells that are called generation i+1.
    This random cell is told to grow to the double of its usual size by setting the cell's parameter area_preferential accordingly.
    This process of area increase is done incrementally, i.e. the preferential area is increased a bit and the PCPVertex model is iterated to an equilibrium in successive steps until the target area 2 A0 is reached.
    Simultaneously to the increase of area the domain size is increased to eventually fit the new cell.
    
    At the point the cell reached the target area in an equilibrated tissue,
    the cell is divided along an axis of division through the cell's center at a random angle. 
    Both cells have the same properties with area_preferential reset to the original value.
    Thereafter, the PCPVertex model is again equilibrated.

    Cell divisions occur at a prescribed probability.


Implementation Details
----------------------

Process ordering
^^^^^^^^^^^^^^^^

#. Cell division


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