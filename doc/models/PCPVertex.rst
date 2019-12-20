
``PCPVertex`` - A vertex-model for planar cell polarity
=======================================================


Model Fundamentals
------------------

This model implements the relaxation of an energy function of a planar cell polarity system towards a local minimum.

The energy function currently includes the following terms
    * area elasticity
    * line tension, alias surface tension

The energy relaxation is performed in one of the following ways
    * along steepest descent with fixed step size

The model accounts for the following topological changes
    * T1 transition: cell intercalation changes neighbourhood of cells
        i.e. an edge shrinks to neglecting length and is replaced with an orthogonal edge connecting the two next neighbour cells. Thus the cells adjacent to the removed edge are no longer neighbours
    * T2 transition: cell extrusion when cell area shrinks below threshold value

The energy function and its gradient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

from Farhadifar et al [2007] and Aigouy et al. [2010]

**Area elasitcity**
    .. math::
        E = \sum_\alpha \frac{K_\alpha}{2}(A_\alpha - A_\alpha^0)^2

    summing over all cells, where

    .. math::
        A_\alpha = 0.5 | \sum_i x_i (y_{i+1} - y_{i-1}) |
        
    using the Shoelace formula for the area of a polygon. Note that the sum is positive if vertices are ordered anti-clockwise and negative if ordered clockwise.

    .. math::
        \Rightarrow F^alpha_i = \frac{K_\alpha}{2} (A_\alpha - A_\alpha^0) \begin{pmatrix} y_{i+1} - y_{i-1} \\ x_{i-1} - x_{i+1} \end{pmatrix}  * \mathrm{sgn}(A_\alpha)

    where :math:`\mathrm{sgn}(A_\alpha)` is 1 if the vertices are ordered anti-clockwise and -1 if ordered clockwise.


**Linetension**
    .. math::    
        E = \sum_{<i,j>} \Lambda_{ij} l_{ij}

    summing over all edges

    .. math::
        \Rightarrow F^{<i,j>}_i = \frac{\Lambda_{<i,j>}}{l_{<i,j>}} \begin{pmatrix}  x_j - x_i \\ y_j - y_i \end{pmatrix}


**Contractility of cell perimeter**
    .. math::
        E = \sum_\alpha \frac{\Gamma_\alpha}{2} L_\alpha^2

    Not implemented


**PCP protein interaction**
    Not implemented


Implementation Details
----------------------

Process ordering
^^^^^^^^^^^^^^^^

#. reset vertex forces, calculate cell area and edge length
#. perform T2 transitions on cells
#. perform T1 transitions on edges 
#. Linetension on edges
#. Area elasticity on cells
#. Update vertex positions on vertices


Default Model Configuration
---------------------------

Below are the default configuration parameters for the ``PCPVertex`` model.

.. literalinclude:: ../../src/models/PCPVertex/PCPVertex_cfg.yml
   :language: yaml
   :start-after: ---


References
----------

* Aigouy, B. et al. Cell Flow Reorients the Axis of Planar Polarity in the Wing Epithelium of Drosophila. 14 (2010) doi:10.1016/j.cell.2010.07.042.
* Etournay, R. et al. Interplay of cell dynamics and epithelial tension during morphogenesis of the Drosophila pupal wing. eLife 4, e07090 (2015).
* Farhadifar, R., Röper, J.-C., Aigouy, B., Eaton, S. & Jülicher, F. The Influence of Cell Mechanics, Cell-Cell Interactions, and Proliferation on Epithelial Packing. Current Biology 17, 2095–2104 (2007).