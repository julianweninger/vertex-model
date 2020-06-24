
``NotchDelta`` - A model for the NotchDelta pathway of cell differentiation
===========================================================================

This is a model for the differentiation of `progenitor` cells to the terminal cell types `hair` and `support` in the NotchDelta pathway.


Model Fundamentals
------------------

In the NotchDelta pathway the `progenitor` cells differentiate to `support` cells at a constant rate and to `hair` cells if the concentration of atoh1 exceeds a threshold value.
Furthermore, `hair` cells suppress the accumulation of atoh1 in neighboring cells.
Eventually, the two terminally differentiated cells interact at their cellular junctions in a way that `support`-`hair` and `support`-`support` contacts are energetically favored over `hair`-`hair` contacts.
Thinking of a real tissue of cells as a liquid on long time-scales, this disfavored `hair`-`hair` junctions lead to what is called T1 transitions -- the exchange of neighborhood, i.e. the elimination of a `hair`-`hair` contact.
In this model, we include these T1 transitions as swapping of states: Whenever a `hair` cell is in contact with an alike cell it swaps position with a third non-`hair` cell.

This leads to the following set of rules

#. atoh1 increases with exponential random distribution ($\lambda$)
#. HCs suppress atoh1 in neighbouring cells
#. Progenitors differentiate to SCs at constant rate
#. Progenitors differentiate to HCs if atoh1 accumulated to threshold value
#. Tag hair cells that have at least one hair cell neighbor.
#. if HC-HC contact tagged, T1 transitions possible, i.e. swap state with one random non HC neighbor


Implementation Details
----------------------

The accumulation of atoh1 is handled by the `Environment` model, which not only implements different types of random distributions, but also different non-uniform ways of accumulating atoh1.

The ``NotchDelta`` allows to define a costume neighborhood.
This is especially helpful if the model is coupled to a more sophisticated description of the cellular arrangement.
To this end, there is also an additional cell type `inactive`.
This allows to handle the cells on a square lattice, but selecting only a subset of the cells.


Default Model Configuration
---------------------------

Below are the default configuration parameters for the ``NotchDelta`` model.

.. literalinclude:: ../../src/models/NotchDelta/NotchDelta_cfg.yml
   :language: yaml
   :start-after: ---


Interpretation of the Model Parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The time-scale :math:`\tau` of the model is defined by the rate of support-cell differentiation ``rate_ps``.
The condition of atoh1 accumulating to a threshold value defines a delay time :math:`\tau_{delay} = \lambda` that is given by the average time to reach this value depending on the chosen mechanism of accumulation.

In the absence of suppression and T1 transitions one can calculate the steady state fraction of `hair`-cells:

.. math::
  \frac{h}{s} = \frac{r_h \exp{(-\tau_{delay}/\tau)}}{r_s + r_h(1-\exp{(-\tau_{delay}/\tau)})} \\

  \frac{h}{s}(\bar{\beta}, \bar{\lambda}) =
    \frac{\bar{\beta} \exp{(-\bar{\lambda})}}
         {1 + \bar{\beta}(1 - \exp{(-\bar{\lambda}}))}

or

.. math::
  h (\bar{\beta}, \bar{\lambda}) =
    \frac{\frac{\bar{\beta} \exp{(-\bar{\lambda})}}
               {1 + \bar{\beta}(1 - \exp{(-\bar{\lambda}}))}}
         {1 + \frac{\bar{\beta} \exp{(-\bar{\lambda})}}
                   {1 + \bar{\beta}(1 - \exp{(-\bar{\lambda}}))}}

  h (\bar{\beta}, \bar{\lambda}) = \frac{\bar{\beta}}{1 + \bar{\beta}} \exp{(-\bar{\lambda})}


References
----------

