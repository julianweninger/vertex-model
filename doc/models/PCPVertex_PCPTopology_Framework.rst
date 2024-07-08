
``PCPVertex / PCPTopology`` - A vertex-model framework with planar cell polarity
================================================================================

Framework Fundamentals
----------------------

The ``PCPVertex`` / ``PCPTopology`` framework provides the infrastructure to obtain a minimized vertex-model (``PCPVertex``) while performing topological operations (``PCPTopology``) on it.
Such an operation could be for example a cell division or external forcing (convergent extension).
Every operation in the ``PCPTopology`` model will be followed by an instruction to the``PCPVertex`` model on how far to iterate towards the work function's minimum.
In this sense, the ``PCPVertex`` model continuously iterates towards the minimum of the work function, while ``PCPTopology`` provides at discrete times operations to the simulated tissue.


Entities - The tissue infrastructure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Vertex models are defined by a set of ``Vertices`` :math:`\{v_i\}`.
Certain of these vertices are connected to form an ``Edge`` :math:`e_<i,j>`.
When egdes together form a closed loop, this defines a ``Cell`` :math:`c_a`.
As ``Utopia::Agents``, their properties are defined by a state, the ``VertexState``, ``EdgeState``, and ``CellState``, respectively.
These 3 classes of entities are organised into 3 ``Utopia::AgentManagers``, living within the ``PCPVertex::EntitiesManager``.

This ``EntitiesManager`` provides the main interface to obtain information on the tissue's configuration.
Here, properties, such as ``position_of(vertex)``, ``length_of(edge)``, ``area_of(cell)``, and many more are defined.
The ``EntitiesManager`` also organises the relations between vertices, edges, and cells through so called ``CustomLinkContainers``, containg most importantly for edges their start and end vertices i and j, and for cells a a vector containing pairs of edges with a ``flip`` boolean.
This vector defines a cell's boundary, starting from vertex i of the first edge, j of that edge will be identical with vertex i from the consecutive edge, where the boolean provides whether the direction of the edge, i.e. :math:`i` and :math:`j` have to be ``flip``-ped.
The vertex :math:`j` of the last edge is identical to the vertex :math:`i` of the first edge, and the closed path is anti-clockwise along the boundary of the cell.
Furthermore, the ``EntitiesManager`` provides secondary topological relations, such of the ``adjoints_of(edge)``, the cells :math:`a` and :math`b` on either side of the edge.

Vertex :math:`i`:
   * ``id()``
   * ``position()`` (also via entities manager)
   * ``state``
      * ``get_force()``
      * ``add_force(force)``
      * parameter interface: Getter and Setter to map {``key`` (string): ``value`` (double)}

Edge :math:`\langle i, j \rangle`:
   * ``id()``
   * ``custom_links()``
      * vertices ``a`` and ``b``
   * ``state``
      * parameter interface: Getter and Setter to map {``key`` (string): ``value`` (double)}

Cell :math:`c_a`:
   * ``id()``
   * ``custom_links()``
      * edge container <``edge``, ``flip``>, defining a closed boundary
   * ``state``
      * ``type``
      * ``area_preferential``
      * parameter interface: Getter and Setter to map {``key`` (string): ``value`` (double)}

The ``EntitiesManager`` resides within the ``PCPVertex`` model an is initialised at initialisation of the PCPVertex model, typically as a hexagonal lattice with {Nx, Ny} cells.


PCPVertex model - The work-function minimizer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The PCPVertex model defines the most basic version of a vertex-model, therefore it contains the EntitiesManager, the work-function, and defines the algorithmic minimization of the work-function.
Every iteration of the vertex-model, moves the vertices towards the minimum of the work-function, including remodelling of junctions and cell extrusions events, so called T1 and T2 transitions, respectively.

Individual work-function terms can be registered to the vertex-model.
Each term, being derived from PCPVertex::WorkFunctionTerm, defines how to calculate the energy on a set of vertices, edges, and cells, and instructions on how to calculate the derivative dW/dx_i (called force), with x_i being the position of vertex i.
The work-function terms is the cumulative result w.r.t energy and force.

Custom work-function terms - the ``WorkFunctionTermBuilder``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Custom work-function terms can be registered to the vertex-model, either directly to the model after its initialisation or through a so-called work-function-builder.
Direct registration adds the initialised work-function to the work-function, which is automatically active.
In contrast, the work-function-builder provides instructions on how to initialise the work-function during the configuration of the PCPVertex model given a set of parameters.
The advantage of providing a builder is that an extensive collection of (custom) work-function terms can be created with unique term-identifiers.
During configuration, i.e. after compilation, the required work-functions can be registered using the model's run-configuration.
A collection of ready-to-use work-function terms can be found in the PCPVertex::WorkFunction namespace.

.. note::

   Work function terms can also be registered through the ``PCPTopology`` model, which mirrors most Getters and Setters of the ``PCPVertex`` model.


PCPTopology model - The operation manager
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the study of biological morphogenesis, developmental processes are occuring altering the internal or external condition of the tissue.
Using vertex-model, the first studies were investigating the imapct of continuous cell proliferation on tissue morphology [Farhadifar et al, Current Biology, 2007, Figure 2].
In the PCPVertex / PCPTopology framework, instructions on when to divide which cell are organised by a seperate model, the ``PCPTopology`` model in so called operations.

``Operations`` are occuring at discrete times associated with the iterative time of the PCPTopology model.
Every operation is associated with instructions when to apply the operation to an ``OperationBundel``.
These instructions also include instructions on how to minimize the work-function using the PCPVertex model.


Operation application
"""""""""""""""""""""

Every operation is therefore applied therefore as follows
   #. Apply operation
   #. Increment time in PCPVertex model
   #. Iterate ``num_steps`` in PCPvertex model
   #. If ``num_repeat > 1``, repeat 4. ``num_repeat`` times.
   #. If ``iterations > 1``, repeat 1.-4. ``iterations`` times.
   
Note, that the time of the ``PCPTopology`` typically increases much slower than that of the ``PCPVertex`` model.
Indeed, if ``num_steps`` and ``num_repeat`` are sufficiently large, the ``PCPVertex`` model reaches equilibrium after every application of the operation.
The ``PCPTopology`` model thus records only the equilibrium states after ``iterations`` applications of the operation, while the ``PCPVertex`` model records the process of how it is approaching this new equilibrium.


PCPVertex minimization
""""""""""""""""""""""

The above iteration of ``PCPVertex`` model using ``num_steps`` and ``num_repeat`` permits to include fluctuations into the minimization process, facilitating the PCPVertex model to find a ''global'' minimum, typically involving a number of T1 transitions.
Fluctuations are implemented in multiple ways, depending on the use case.

**Jiggle vertices**.
``jiggle_intensity``, displace randomly all vertices at a length scale :math:`l = I * \sqrt{\langle A \rangle}`, where :math:`I` is the `jiggle_intensity` and :math:`\langle A \rangle` the mean area of a cell.
The jiggle approach is useful to analyse a truly quasi-static development. 
To do so the ``PCPVertex`` model is iterated over ``num_steps`` sufficiently large to reach a new equilibrium after every jiggle.
Indeed, a ``jiggle_tolerance`` can be defined, when the energy in the PCPVertex model changes by less than this tolerance, the model is not iterated further but the next jiggle performed.
After ``num_repeat`` applications of the jiggling, a last minimization is performed, where energy change should fall below ``tolerance``.

**Linetension parameter fluctuations**.
Every junction receives an additional linetension parameter, which develops according to an Ornstein Uhlenbeck process, with zero mean and amplitude ``linetension_fluctuations``.
The parameter's characteristic time is given by ``linetension_fluctuations_tau``, where parameter self-correlation decay at this characteristic time.
This approach is particularly useful as it defines a timescale of the ``PCPVertex`` model.

In fact, multiple timescales can be defined with linetension parameter fluctuations.
First there is the numeric timescale, :math:`\tau_\mathrm{num} = 1 / dt`, defining the numeric step size.
Second, the ``PCPVertex`` model has an ''elastic'' timescale :math:`\tau`.
This is the timescale at which the ``PCPVertex`` model approaches an energetic minimum with no or few T1 transitions.
It capture the elastic material properties of the simulated tissue and depends on the parameters of the work-function.
A third timescale is then explicitly defined by the ``linetension_fluctuations_tau`` :math:`\tau_\Lambda`, typically being (much) larger than :math:`\tau`.
Eventually, for sufficiently large ``linetension_fluctuations`` amplitude, a diffusive timescale :math:`\tau_D` can be measured, where cells move over distances :math:`d = \sqrt{2 D t}` larger than a cell's diameter.


Details on operations
"""""""""""""""""""""

A sequence of operations can be defined within the ``PCPTopology`` model.
In this case, in every iteration of the ``PCPTopology`` model, the operations are called in sequence.
However, operations do not need to be active at every time.
A list of ``times`` is provided to every operation; it is possible to define ``times`` using a range-generator using ``start``, ``stop``, and ``step``.
Only at those times the operation is active, otherwise it is skipped.
Additionally, ``iterations`` can be specified to apply this operation multiple times before moving to the next operation in sequence.
After every application of an operation, the minimization instructions are followed, which can be customized for every operation.
In particular, minimization can be applied after ``every`` application, ``once`` after the operation instructions have been applied ``iteration`` times (without minimzing between applications), or turned ``off`` not minimizing the work function.
These can be specified using the minimization ``mode``.
Additionally, operations can be applied during the model's prolog, prior to ``t = 0``, using ``iterations_prolog``, or after ``t = t_\mathrm{max}`` using the ``iterations_epilog``; following the same logic as ``iterations``.

Operations can be specified using the run-configuration as follows
 
.. code-block:: yaml
   
   ---
   parameter_space:
      PCPTopology:
         # Defaults to minimization instructions
         minimzation:
            # numeric step size
            dt: 0.05

            # How many steps per minimization
            num_steps: 20

            # How often to repeat minimization
            num_repeat: 3

            # other parameters can go here, e.g. 
            #     jiggle_intensity, jiggle_tolerance, and tolerance
            #     linetension_fluctuations and linetension_fluctuations_tau

         operations:
            - my_first_operation:
               # The times in PCPTopology when to apply this operation
               times: [1, 10, 20, 21, 22, 23]
               
               # At every time it is applied 10 times
               iterations: 10

               # Put your parameters here
               parameter_a: 10
               parameter_b: 1

               minimization:
                  # Minimize after every iteration
                  mode: every

            - the_second_operation:
               # The times in PCPTopology when to apply this operation
               times:
                  start: 1
                  stop: 11
                  step: 2
                  # -> [1, 3, 5, 7, 9]
               
               # At every time it is applied 10 times
               iterations: 3

               # Put your parameters here
               parameter_c: 10
               parameter_d: 1

               minimization:
                  # Minimize once after all iterations
                  mode: once

                  # Howver, do 20 repeats
                  num_repeat: 20
                  
            - the_third_operation:
               # The times in PCPTopology when to apply this operation
               times: [1]
               
               # Apply 3 times during prolog, prior to t = 0
               iterations_prolog: 3

               # Put your parameters here
               parameter_e: 10
               parameter_f: 1

               minimization:
                  # Never minimize
                  mode: off

In this instruction, the following will happen:
In PCPTopology prolog, prior to t = 0, 3 applications of ``the_third_operation``, no minimization of the PCPVertex model (t = 0).
At PCPTopology t = 1, 10 iterations of ``my_first_operation`` with 3 times 20 steps performed by PCPVertex model (op1, t -> 1, t->61, op1, t-> 62, t->122, ..., op1, t-> 610), then 3 iterations of ``the_second_operation`` with 20 x 20 steps once, (op2, op2, op2, t->611, t->631, ..., t->1011), and eventually 1 iteration of ``the_third_operation``, no iteration of the PCPVertex model, (op3, t -> 1011).
At PCPTopology t = 2, only ``my_first_operation`` is active, with 10 iterations including 3 times 20 steps performed by PCPVertex model after each application.
Eventually, at PCPTopology t = 3, 10 iterations of ``my_first_operation`` with 3 times 20 steps performed by PCPVertex model, followed by 3 iterations of ``the_second_operation`` with 20 x 20 steps once.


Custom operations - The ``OperationBundelBuilder``
""""""""""""""""""""""""""""""""""""""""""""""""""

A collection of operations is provided by the ``PCPTopology`` model.
Furthermore, it is possible to register custom operations.
Similar to the registration of work-function terms with the PCPVertex model, ``OperationBundel`` can be registered with the PCPTopology model directly or via an builder.
Direct registration adds the initialised operation at the end of the operation sequence.
In contrast, the operation-builder provides instructions on how to initialise the operation during the configuration of the PCPTopology model given a set of parameters.


The ``CopyMeVertex`` model - Building a custom model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you are familiar with running and adjusting operations and work-function terms of the `PCPVertex` and `PCPTopology` model, you can write your own customized operations and terms.
In that case, I propose to use the `register_operation` and `register_work_function_builder` interface of the PCPTopology model respectively.
An example of how this can be done is given by the `CopyMeVertex` model.
Copy this model and rename `CopyMeVertex` -> `YourModelName`.

The `CopyMeVertex` model is just a `PCPTopololgy` model with customized name.
However, it permits to also define custom work-fucntion terms and operations. Their implementation and intergration into respective Builders are exemplarily shown in the `CopyMeVertex` model.
So, copy-paste and modify to obtain your own first custom terms and operations.

Enjoy!


Default Model Configuration
---------------------------

Below are the default configuration parameters for the ``PCPTopology`` model.

.. literalinclude:: ../../src/models/PCPTopology/PCPTopology_cfg.yml
   :language: yaml
   :start-after: ---
