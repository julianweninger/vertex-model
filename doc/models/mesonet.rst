
``mesonet`` — Model of Consumer-Resource Systems
================================================

``mesonet`` models the interaction of consumer species with the resources in their environment.

This implementation is preceded by a `Python implementation <https://ts-gitlab.iup.uni-heidelberg.de/yunus/mesonet>`_ of the same name, where the interactions were represented using a network.
The present Utopia implementation foregoes the network approach; while this restricts interactions to the consumer-resource scenario, the implementation is simpler, faster, and allows to explore other aspects of the model.

*This document is very much WIP*


.. contents::
   :local:
   :depth: 2

----


Default Configurations
----------------------
Below are the default model and plot configuration parameters for the ``mesonet`` model.

Default Model Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../src/models/mesonet/mesonet_cfg.yml
   :language: yaml
   :start-after: ---

Default Plot Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../src/models/mesonet/mesonet_plots.yml
   :language: yaml
   :start-after: ---

Base Plot Configuration
^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../src/models/mesonet/mesonet_base_plots.yml
   :language: yaml
   :start-after: ---

Note that some of these plots are based on the utopya base plot configurations, which are not part of this documentation.
