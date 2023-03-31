// /// Minimize using alternating brownian motion and minimization
// /** Iterates the vertex model using a brownian motion update scheme,
//  *  then minimizes energy using the usual configuration.
//  * 
//  *  The following parameter are extracted from cfg 
//  *  (besides those passed to `OperationParams`):
//  *      - `num_repeat` (std::size_t): How often to repeat the cycle
//  *      - `brownian_iteration` (minimization Config): The configuration for
//  *              brownian motion minimization update.
//  *      - `minimization` (minimization Config): The update to the default
//  *              minimization config. Used to minimize after brownian motion.
//  */
// template<typename Logger>
// OperationBundle build_brownian_noise (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params,
//         std::shared_ptr<Logger> logger = nullptr,
//         std::function<void()> monitor = [](){ })
// {
//     using MinimizationMode = OperationParams::MinimizationMode;

//     OperationParams params(name, cfg, default_minim_params);

//     MinimizationParams brownian_iteration(
//         get_as<Config>("brownian_iteration", cfg),
//         params.minimization_params);
//     auto minimization = params.minimization_params;
//     bool minimize = (    params.minimization_mode == MinimizationMode::Every
//                      and not params.disable);

//     auto num_repeat = get_as<std::size_t>("num_repeat", cfg);
//     params.add_num_minimisations = 2 * num_repeat;

//     Operation operation = [num_repeat, brownian_iteration, minimization,
//                            minimize,
//                            monitor, logger]
//                           (PCPVertex& vertex_model)
//     {
//         logger->debug("Repeating cycle of minimization and brownian motion {} "
//                       "times ..", num_repeat);
//         for (std::size_t i = 0; i < num_repeat; i++) {
//             logger->trace("   Minimizing energy {} / {} ..", i, num_repeat);
//             if (minimize) {
//                 vertex_model.minimize_energy(minimization, monitor);
//             }
            
//             logger->trace("   Iterating brownian motion {} / {} ..",
//                           i, num_repeat);
//             vertex_model.minimize_energy(brownian_iteration, monitor);
//         }
//     };

//     return std::make_pair(operation, params);
// }




// /// Enable or disable topological transitions
// /** The following parameter are extracted from cfg 
//  *  (besides those passed to `OperationParams`):
//  *      - `enable_T1s` (bool): Whether to enable T1 transitions
//  *      - `enable_T2s` (bool): Whether to enable T2 transitions
//  */
// OperationBundle build_enable_transitions (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     auto enable_T1s = get_as<bool>("enable_T1s", cfg);
//     auto enable_T2s = get_as<bool>("enable_T2s", cfg);

//     Operation operation = [enable_T1s, enable_T2s](PCPVertex& vertex_model)
//     {
//         vertex_model.enable_transitions(enable_T1s, enable_T2s);
//     };

//     return std::make_pair(operation, params);
// }

// /// An operation to fix boundary vertices in space
// /** Sets the fix_boundary in the vertex model and is updated continuously
//  *  in update
//  * 
//  *  The following parameter are extracted from cfg 
//  *  (besides those passed to `OperationParams`):
//  *      - `fix_boundary` (bool, default: true): Whether to fix the boundary
//  *                                              in space
//  */
// OperationBundle build_fix_boundary (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     bool fix_boundary = get_as<bool>("fix_boundary", cfg, true);

//     Operation operation = [fix_boundary] (PCPVertex& vertex_model)
//     {
//         vertex_model.fix_boundary(fix_boundary);
//     };

//     return std::make_pair(operation, params);
// }

// /// The operation to increment preferential area
// /** \details The following parameter are extracted from cfg 
//  *           (besides those passed to `OperationParams`):
//  *               - `progenitor` (double, default: 0.): incremental value for
//  *                      cells of type progenitor
//  *               - `hair` (double, default: 0.): incremental value for cells
//  *                      of type hair
//  *               - `support` (double, default: 0.): incremental value for cells
//  *                      of type support
//  *               - `adapt_support` (bool, default: false): If true, the total
//  *                      area of hair and support cells remains constant, hence
//  *                      support cells compensate area changes from hair cells.
//  *               - `relax_domain` (bool, default: false): If true, the domain
//  *                      size is adapted to fit the cells as per preferential 
//  *                      area.
//  *               - `gradient_hair_max` and `gradient_hair_min`
//  *                  (double, default: 0.): The difference between
//  *                      the cells that grow the most and the least. Using a
//  *                      negative slope for non-periodic and a cos(x)**2 in 
//  *                      periodic BC along the horizontal (P-D) axis.
//  *                      If `adapt_support`, the total A0 of all cells remains
//  *                      constant (uniform substraction of the change)
//  *               - `hair_gradient_center` (double, default: 0.5): The relative 
//  *                      position of dA = 0 in non-periodic BC and the relative 
//  *                      position of dA = max in periodic BC.
//  * 
//  *  \note This affects the current entities properties, it does not overwrite
//  *        changes in the past.
//  */
// OperationBundle build_increment_area (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double prog(get_as<double>("progenitor", cfg, 0.));
//     double hair(get_as<double>("hair", cfg, 0.));
//     double support(get_as<double>("support", cfg, 0.));
//     bool adapt_support(get_as<bool>("adapt_support", cfg, false));
//     bool relax_domain(get_as<bool>("relax_domain", cfg, false));
//     bool relax_domain_PD_axis(get_as<bool>("relax_domain_PD_axis",
//                                            cfg, false));
//     double gradient_hair_max(get_as<double>("gradient_hair_max", cfg, 0.));
//     double gradient_hair_min(get_as<double>("gradient_hair_min", cfg, 0.));
//     double hair_gradient_center(get_as<double>("hair_gradient_center",
//                                                cfg, 0.5));
    
//     if (adapt_support and fabs(support) > 1.e-12) {
//         throw std::invalid_argument(fmt::format(
//             "In cfg {}: If `adapt_support = True`, "
//             "then the increment value for the support cells "
//             "`support` must be zero, but was {}!", name, support));
//     }

//     if (adapt_support + relax_domain + relax_domain_PD_axis > 1) {
//         throw std::invalid_argument(fmt::format(
//             "In cfg {}: `adapt_support` and `relax_domain` cannot both be "
//             "true!", name));
//     }

//     Operation operation = [prog, hair, support, adapt_support,
//                            gradient_hair_max, gradient_hair_min,
//                            hair_gradient_center,
//                            relax_domain, relax_domain_PD_axis]
//             (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;
//         using SpaceVec = PCPVertex::SpaceVec;

//         const auto& am = vertex_model.get_am();
//         const auto& cells = am.cells();

//         double _support = support;

//         if (adapt_support) {
//             std::size_t num_hcs = std::count_if(
//                 cells.begin(), cells.end(), 
//                 [](const auto& c) {
//                     return c->state.type == CellType::hair;
//                 });

//             if (num_hcs < cells.size()) {
//                 _support = -1. * hair * num_hcs / (cells.size() - num_hcs);
//             }
//         }

//         PCPVertex::RuleFuncCell update = [vertex_model, prog, hair, _support]
//                 (const auto& cell)
//         {
//             auto state = cell->state;
            
//             double dA0 = 0.;
//             if (state.type == CellType::progenitor) {
//                 dA0 = prog;
//             }
//             else if (state.type == CellType::hair) {
//                 dA0 = hair;
//             }
//             else if (state.type == CellType::support) {
//                 dA0 = _support;
//             }

//             state._area_preferential += dA0;
//             return state;
//         };
        
//         apply_rule<Update::sync>(update, cells);


//         // add spatial gradient in HC species
//         if (fabs(gradient_hair_max) > 1.e-12) {
//             std::function<double(const SpaceVec&)> map;
//             if (not vertex_model.get_space()->periodic) {
//                 auto [min, max] = am.get_extent();

//                 double L = (max - min)[0];
//                 double x_min = min[0];

//                 map = 
//                     [gradient_hair_max, gradient_hair_min,
//                      hair_gradient_center, L, x_min]
//                     (const SpaceVec& pos)
//                     {
//                         double hair_gradient = (  gradient_hair_max
//                                                 - gradient_hair_min);
//                         double rel_x = (pos[0] - x_min) / L;
//                         double fac = (hair_gradient_center - rel_x);
//                         return hair_gradient * fac + gradient_hair_min;
//                     };
//             }
//             else {
//                 SpaceVec domain = am.get_space()->get_domain_size();
//                 double L = domain[0];
//                 map = 
//                     [gradient_hair_max, gradient_hair_min,
//                      hair_gradient_center, L]
//                     (const SpaceVec& pos)
//                     {
//                         double hair_gradient = (  gradient_hair_max
//                                                 - gradient_hair_min);
//                         double rel_x = (pos[0] / L - hair_gradient_center);
//                         double fac = std::pow(cos(rel_x * M_PI), 2);
//                         return hair_gradient * fac + gradient_hair_min;
//                     };
//             }

//             double tissue_area = std::accumulate(
//                 cells.begin(), cells.end(), 0.,
//                 [](const double& val, const auto& cell){
//                     return val + cell->state._area_preferential;
//                 });

//             std::size_t num_SC = 0;
//             for (const auto& cell : cells) {
//                 if (cell->state.type == CellType::hair) {
//                     SpaceVec pos = am.barycenter_of(cell);
//                     double dA = map(pos);
//                     cell->state._area_preferential += dA;
//                 }
//                 else {
//                     ++num_SC;
//                 }
//             }
//             if (adapt_support and num_SC > 0) {
//                 double new_tissue_area = std::accumulate(
//                     cells.begin(), cells.end(), 0.,
//                     [](const double& val, const auto& cell){
//                         return val + cell->state._area_preferential;
//                     });

//                 for (const auto& cell : cells) {
//                     cell->state._area_preferential -= ( 
//                         (new_tissue_area - tissue_area) / cells.size()
//                     );
//                 }
//             }
//         }

//         if (relax_domain) {
//             double area = std::accumulate(cells.begin(), cells.end(), 0.,
//                                 [](const double& val, const auto& cell) {
//                                     return val + cell->state._area_preferential;
//                                 });

//             PCPVertex::SpaceVec domain_size = 
//                     vertex_model.get_space()->get_domain_size();
//             vertex_model.increase_domain_size(
//                 area - domain_size[0]*domain_size[1]);
//         }
//         else if (relax_domain_PD_axis) {
//             double area = std::accumulate(cells.begin(), cells.end(), 0.,
//                             [](const double& val, const auto& cell) {
//                                 return val + cell->state._area_preferential; });

//             PCPVertex::SpaceVec domain_size = 
//                     vertex_model.get_space()->get_domain_size();
//             double dA = area - domain_size[0]*domain_size[1];

//             vertex_model.stretch_domain({dA/domain_size[1], 0.},
//                                         false, false, false);
//         }

//         return;
//     };

//     return std::make_pair(operation, params);
// }


// /// The operation to locally relax SC to the available area constraint
// /** SC target area is adapted such that the local mean (of cells within a 
//  *  given Manhatten-distance) of target area is relaxed to a given area.
//  *  Furthermore, relax target area to the average of local mean in SC 
//  *  target area:
//  * 
//  *  \f$ dA_\alpha^(i+1) = A_\alpha^i
//  *                        - (<A>_N / A_0 - 1) * A_\alpha^i / \tau
//  *                        + (<A>_{S \in nbs} - A_\alpha^i) / \tau
//  *  \f$ 
//  * 
//  *  Extracts the following cfgs:
//  *      - `adapt_SC_tau` (double): the characteristic time of relaxation in 
//  *              units of iterations.
//  *      - `use_area_fluctuations` (bool): The target area change is 
//  *              in area_fluctuations. Requires a minimization with 
//  *              Ornstein-Uhlenbeck fluctuations of target area
//  *      - `target_area` (double): The average target area of a cell in a 
//  *              patch. The value to which SC area is adapted
//  *      - `manhatten_distance` (uint): The distance to the furthest cell to 
//  *              include in neighborhood. Distance in unit of cell count.
//  *      - `min_area` (double): A minimum target area.
//  *          
//  */
// OperationBundle build_relax_SC_area (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double adapt_SC_tau(get_as<double>("adapt_SC_tau", cfg));
//     bool use_fluctuations(get_as<bool>("use_area_fluctuations", cfg));
//     double target_area(get_as<double>("target_area", cfg));
//     std::size_t distance(get_as<std::size_t>("manhatten_distance", cfg));
//     double min(get_as<double>("min_area", cfg));
//     // double center(get_as<double>("gradient_center", cfg, 0.5));

//     Operation operation = [adapt_SC_tau, target_area, distance, min,
//                            use_fluctuations]
//             (PCPVertex& vertex_model)
//     {
//         using Cell = PCPVertex::Cell;
//         using CellType = PCPVertex::CellType;

//         const auto& am = vertex_model.get_am();
//         const auto& cells = am.cells();

//         std::map<std::shared_ptr<Cell>, double> targets({});
//         for (const auto& cell : cells) {
//             if (cell->state.type == CellType::hair) {
//                 continue;
//             }

//             const auto& neighbors = am.neighbors_of(cell, distance);
            
//             std::size_t N = neighbors.size() + 1;
//             std::size_t N_SC = 1;
//             double average_area = cell->state._area_preferential;
//             double average_area_sc = cell->state._area_preferential;
//             for (const auto& n : neighbors) {
//                 average_area += n->state._area_preferential;
                
//                 if (n->state.type != CellType::hair) {
//                     average_area_sc += n->state._area_preferential;
//                     N_SC++;
//                 }
//             }
//             average_area /= static_cast<double>(N);
//             average_area_sc /= static_cast<double>(N_SC);

//             auto A = cell->state._area_preferential;
//             double A_target = A - (average_area / target_area - 1.) * A / adapt_SC_tau;
//             A_target += (average_area_sc - A) / adapt_SC_tau;

//             targets[cell] = std::max(A_target, min);
//         }

//         for (const auto& [cell, A_target] : targets) {
//             double dA = A_target - cell->state._area_preferential;
//             cell->state._area_preferential += dA;
//             if (use_fluctuations) {
//                 cell->state._area_preferential_fluctuations -= dA;
//             }
//         }

//         return;
//     };

//     return std::make_pair(operation, params);
// }

// OperationBundle build_increment_area_gradient (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double init(get_as<double>("initialise_homogeneous_area", cfg));
//     double incr_hom(get_as<double>("increment_homogeneous_HC", cfg));
//     double incr_grad(get_as<double>("increment_HC_gradient", cfg));
    
//     bool use_fluctuations(get_as<bool>("use_area_fluctuations", cfg));
//     double min(get_as<double>("min_area", cfg));

//     auto homogeneous = std::make_shared<double>(init);
//     auto gradient = std::make_shared<double>(0.);

//     auto __relax_SC_bundle = build_relax_SC_area("relax_SC", cfg,
//                                                  default_minim_params);
//     auto relax_SC = std::get<0>(__relax_SC_bundle);

//     Operation operation = [homogeneous, gradient, incr_hom, incr_grad,
//                            relax_SC, min, use_fluctuations]
//             (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;
//         using SpaceVec = PCPVertex::SpaceVec;

//         const auto& am = vertex_model.get_am();
//         const auto& cells = am.cells();

//         *homogeneous += incr_hom;
//         *gradient += incr_grad;

//         // add spatial gradient in HC species
//         std::function<double(const SpaceVec&)> map;
//         if (not vertex_model.get_space()->periodic) {
//             auto [min, max] = am.get_extent();

//             double L = (max - min)[0];
//             double x_min = min[0];

//             map = 
//                 [homogeneous, gradient, L, x_min]
//                 (const SpaceVec& pos)
//                 {
//                     double rel_x = (pos[0] - x_min) / L;

//                     return *gradient * (1. - rel_x) + *homogeneous;
//                 };
//         }
//         else {
//             SpaceVec domain = am.get_space()->get_domain_size();
//             double L = domain[0];
//             map = 
//                 [homogeneous, gradient, L]
//                 (const SpaceVec& pos)
//                 {
//                     double rel_x = pos[0] / L;

//                     return (  *gradient * std::pow(cos(rel_x * M_PI), 2)
//                             + *homogeneous);
//                 };
//         }

//         // fix the area of HC according to position
//         for (const auto& cell : cells) {
//             if (cell->state.type != CellType::hair) {
//                 continue;
//             }

//             SpaceVec pos = am.barycenter_of(cell);
            
//             double A = std::max(map(pos), min);
            
//             double dA = A - cell->state._area_preferential;
//             cell->state._area_preferential += dA;
//             if (use_fluctuations) {
//                 cell->state._area_preferential_fluctuations -= dA;
//             }
//         }

//         relax_SC(vertex_model);

//         return;
//     };

//     return std::make_pair(operation, params);
// }


// /// The operation to increment cell contractility
// /** \details The following parameter are extracted from cfg 
//  *           (besides those passed to `OperationParams`):
//  *               - `progenitor` (double, default: 0.): incremental value for
//  *                      cells of type progenitor
//  *               - `hair` (double, default: 0.): incremental value for cells
//  *                      of type hair
//  *               - `support` (double, default: 0.): incremental value for cells
//  *                      of type support
//  * 
//  *  \note This affects the current entities properties, it does not overwrite
//  *        changes in the past.
//  */
// OperationBundle build_increment_cell_contractility (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double incr_prog(get_as<double>("progenitor", cfg, 0.));
//     double incr_hair(get_as<double>("hair", cfg, 0.));
//     double incr_support(get_as<double>("support", cfg, 0.));

//     Operation operation = [incr_prog, incr_hair, incr_support]
//             (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;

//         const auto& am = vertex_model.get_am();
//         const auto& cells = am.cells();
        
//         PCPVertex::RuleFuncCell update = [incr_prog, incr_hair, incr_support]
//                 (const auto& cell)
//         {
//             auto state = cell->state;

//             if (state.type == CellType::progenitor) {
//                 state.contractility += incr_prog;
//             }
//             else if (state.type == CellType::hair) {
//                 state.contractility += incr_hair;
//             }
//             else if (state.type == CellType::support) {
//                 state.contractility += incr_support;
//             }

//             return state;
//         };
        
//         apply_rule<Update::sync>(update, cells);

//         return;
//     };

//     return std::make_pair(operation, params);
// }

// /// The operation to increment the domain size
// /** The following parameter are extracted from cfg 
//  *  (besides those passed to `OperationParams`):
//  *      - `value` (SpaceVec<2> or float): The incremental increase in x and y,
//  *              or the incremental area (ratio Lx to Ly constant)
//  *      - `compensate` (bool): Whether to compensate the increase in area in the
//  *              preferential area of the cells. If true, keeps the mechanical
//  *              properties constant
//  *      - `fix_hair_cell_area` (bool, default: false): Whether to fix the area
//  *              of cells of type hair.
//  */
// OperationBundle build_increment_domain (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     using SpaceVec = PCPVertex::SpaceVec;

//     OperationParams params(name, cfg, default_minim_params);

//     if (cfg["value"] and cfg["value"].IsSequence()) {
//         SpaceVec increment(get_as_SpaceVec<2>("value", cfg));
//         bool compensate(get_as<bool>("compensate", cfg));
//         bool fix_hc_area(get_as<bool>("fix_hair_cell_area", cfg, false));
//         bool fix_sc_area(get_as<bool>("fix_support_cell_area", cfg, false));

//         if (compensate and fix_hc_area and fix_sc_area) {
//             throw std::invalid_argument("Cannot compensate area while fixing "
//                 "hair and support cell area in operation 'increment_domain'!");
//         }
        
//         Operation operation = [increment, compensate, fix_hc_area,
//                                fix_sc_area] (PCPVertex& vertex_model)
//         {
//             vertex_model.stretch_domain(increment, compensate,
//                                         fix_hc_area, fix_sc_area);
//         };

//         return std::make_pair(operation, params);
//     }
//     else {
//         double increment(get_as<double>("value", cfg));

//         Operation operation = [increment] (PCPVertex& vertex_model)
//         {
//             vertex_model.increase_domain_size(increment);
//         };

//         return std::make_pair(operation, params);
//     }
// }

// /// The operation to increment edge contractility
// /** The following parameter are extracted from cfg 
//  *  (besides those passed to `OperationParams`):
//  *      - `progenitor_progenitor` (double, default: 0.): Incremental value for 
//  *              edges between cells of types progenitor and progenitor.
//  *      - `progenitor_hair` (double, default: 0.): Incremental value for 
//  *              edges between cells of types progenitor and hair.
//  *      - `progenitor_support` (double, default: 0.): Incremental value for 
//  *              edges between cells of types progenitor and support.
//  *      - `hair_hair` (double, default: 0.): Incremental value for 
//  *              edges between cells of types hair and hair.
//  *      - `hair_support` (double, default: 0.): Incremental value for 
//  *              edges between cells of types hair and support.
//  *      - `support_support` (double, default: 0.): Incremental value for 
//  *              edges between cells of types support and support.
//  * 
//  *  The following paramter are excluded and passed on to the Vertex model
//  *      - `set_ppMLC` (double, default: 0.): The contractility of SC-SC
//  *              junctions with a specified orientation.
//  *              See PCPVertex::set_ppMLC_contractility.
//  *      - `set_ppMLC_angle` (double, default: 0.): The angle of the ppMLC
//  *              contractility orientation, where 0 is the (1, 0) axis and
//  *              Pi the (-1, 0) axis.
//  *              See PCPVertex::set_ppMLC_contractility.
//  *      - `set_pMLC` (double, default: 0.): The contractility of HC-SC junctions
//  *              where the orientation is specified by the HC's polarity vector.
//  *              See PCPVertex::set_pMLC_contractility.
//  *              
//  * 
//  *  \note This increments the matrix of contractility, hence affects current and
//  *        future edges between cells of corresponding type.
//  */
// OperationBundle build_increment_edge_contractility (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     using SpaceVec = PCPVertex::SpaceVec;

//     OperationParams params(name, cfg, default_minim_params);
    
//     double incr_prog_prog(get_as<double>("progenitor_progenitor", cfg, 0.));
//     double incr_prog_hair(get_as<double>("progenitor_hair", cfg, 0.));
//     double incr_prog_supp(get_as<double>("progenitor_support", cfg, 0.));
//     double incr_prog_bnd(get_as<double>("progenitor_boundary", cfg, 0.));
//     double incr_hair_hair(get_as<double>("hair_hair", cfg, 0.));
//     double incr_hair_supp(get_as<double>("hair_support", cfg, 0.));
//     double incr_hair_bnd(get_as<double>("hair_boundary", cfg, 0.));
//     double incr_supp_supp(get_as<double>("support_support", cfg, 0.));
//     double incr_supp_bnd(get_as<double>("support_boundary", cfg, 0.));

//     double ppMLC(get_as<double>("set_ppMLC", cfg, 0.));
//     double ppMLC_angle(get_as<double>("set_ppMLC_angle", cfg, 0.));
//     auto ppMLC_axis = SpaceVec({cos(ppMLC_angle), sin(ppMLC_angle)});
//     bool ppMLC_curved_axis(get_as<bool>("ppMLC_curved_axis", cfg, false));

//     double pMLC(get_as<double>("set_pMLC", cfg, 0.));
//     double adaptivity(get_as<double>("polarity_adaptivity", cfg, 1));
//     auto reset_HC_polarity(get_as<std::pair<bool, double>>(
//         "reset_HC_polarity", cfg, std::make_pair(false, -M_PI_2)));
//     std::size_t reset_HC_polarity_antispin_nshift = get_as<std::size_t>(
//         "reset_HC_polarity_antispin_nshift", cfg, 0);
//     std::size_t reset_HC_polarity_antispin_n = get_as<std::size_t>(
//         "reset_HC_polarity_antispin_n", cfg, 0);
//     auto reset_HC_polarity_shuffle(get_as<double>(
//         "reset_HC_polarity_shuffle", cfg, 0.));

//     if (cfg["hair_progenitor"]) {
//         double incr_hair_prog = get_as<double>("hair_progenitor", cfg);
//         if (cfg["progenitor_hair"] and incr_hair_prog != incr_prog_hair) {
//             throw KeyError("", cfg, fmt::format("Increment has to be "
//                 "symmetric matrix."));
//         }
//         incr_prog_hair = incr_hair_prog;
//     }
//     if (cfg["support_progenitor"]) {
//         double incr_supp_prog = get_as<double>("hair_progenitor", cfg);
//         if (cfg["progenitor_support"] and incr_supp_prog != incr_prog_supp) {
//             throw KeyError("", cfg, fmt::format("Increment has to be "
//                 "symmetric matrix."));
//         }
//         incr_prog_supp = incr_supp_prog;
//     }
//     if (cfg["support_hair"]) {
//         double incr_supp_hair = get_as<double>("hair_progenitor", cfg);
//         if (cfg["hair_support"] and incr_supp_hair != incr_hair_supp) {
//             throw KeyError("", cfg, fmt::format("Increment has to be "
//                 "symmetric matrix."));
//         }
//         incr_hair_supp = incr_supp_hair;
//     }

//     Operation operation = [incr_prog_prog, incr_prog_hair, incr_prog_supp,
//                            incr_hair_hair, incr_hair_supp, incr_supp_supp,
//                            incr_prog_bnd, incr_hair_bnd, incr_supp_bnd,
//                            ppMLC, ppMLC_axis, ppMLC_curved_axis, pMLC,
//                            adaptivity, reset_HC_polarity, 
//                            reset_HC_polarity_antispin_nshift,
//                            reset_HC_polarity_antispin_n,
//                            reset_HC_polarity_shuffle]
//             (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;

//         auto contractility = vertex_model.get_edge_contractility();

//         contractility(CellType::progenitor,
//                       CellType::progenitor) += incr_prog_prog;
//         contractility(CellType::progenitor, CellType::hair) += incr_prog_hair;
//         contractility(CellType::progenitor,
//                       CellType::support) += incr_prog_supp;
//         contractility(CellType::progenitor,
//                       CellType::num_cell_types) += incr_prog_bnd;
        
//         contractility(CellType::hair, CellType::progenitor) += incr_prog_hair;
//         contractility(CellType::hair, CellType::hair) += incr_hair_hair;
//         contractility(CellType::hair, CellType::support) += incr_hair_supp;
//         contractility(CellType::hair,
//                       CellType::num_cell_types) += incr_hair_bnd;
        
//         contractility(CellType::support,
//                       CellType::progenitor) += incr_prog_supp;
//         contractility(CellType::support, CellType::hair) += incr_hair_supp;
//         contractility(CellType::support, CellType::support) += incr_supp_supp;
//         contractility(CellType::support,
//                       CellType::num_cell_types) += incr_supp_bnd;
                      
//         contractility(CellType::num_cell_types,
//                       CellType::progenitor) += incr_prog_bnd;
//         contractility(CellType::num_cell_types,
//                       CellType::hair) += incr_hair_bnd;
//         contractility(CellType::num_cell_types,
//                       CellType::support) += incr_supp_bnd;

//         // set the contractility and update the edge properties
//         vertex_model.set_edge_contractility(contractility, true);

//         vertex_model.set_ppMLC_contractility(ppMLC, ppMLC_axis,
//                                              ppMLC_curved_axis);
//         vertex_model.set_pMLC_contractility(pMLC, adaptivity);

//         if (std::get<bool>(reset_HC_polarity)) {
//             const auto& cells = vertex_model.get_am().cells();

//             std::size_t nshift = reset_HC_polarity_antispin_nshift;
//             std::size_t n = reset_HC_polarity_antispin_n;

//             std::size_t cell_cnt = 0;
//             for (const auto& cell : cells) {
//                 if (cell->state.type == CellType::support) {
//                     continue;
//                 }
//                 auto state = cell->state;
//                 state._polarity = std::get<double>(reset_HC_polarity);
//                 state._polarity_fluctuations = 0.;

//                 std::size_t flip = 0;
//                 // in a flip row
//                 if (nshift > 0 and (cell_cnt / nshift) % 2 == 1) {
//                     flip++;
//                 }
//                 // in a row flip each nth
//                 if (n > 0 and cell_cnt % n == n - 1) {
//                     flip++;
//                 }
//                 if (flip % 2 == 1) {
//                     state._polarity = std::fmod(state._polarity + M_PI, 
//                                                 2 * M_PI);
//                 }
//                 cell->state = state;

//                 cell_cnt += 1;
//             }
//         }

//         if (fabs(reset_HC_polarity_shuffle) > 1.e-12) {
//             auto distr = std::uniform_real_distribution<double>(
//                 -reset_HC_polarity_shuffle, reset_HC_polarity_shuffle);
//             const auto& cells = vertex_model.get_am().cells();
//             const auto& rng = vertex_model.get_rng();
//             for (const auto& cell : cells) {
//                 cell->state._polarity = std::fmod(
//                     cell->state._polarity + distr(*rng) + 2 * M_PI, 
//                     2 * M_PI);
//                 cell->state._polarity_fluctuations = 0.;
//             }
//         }


//     };

//     return std::make_pair(operation, params);
// }

// /// The operation to increment linetension
// /** The following parameter are extracted from cfg 
//  *  (besides those passed to `OperationParams`):
//  *      - `progenitor_progenitor` (double, default: 0.): Incremental value for 
//  *              edges between cells of types progenitor and progenitor.
//  *      - `progenitor_hair` (double, default: 0.): Incremental value for 
//  *              edges between cells of types progenitor and hair.
//  *      - `progenitor_support` (double, default: 0.): Incremental value for 
//  *              edges between cells of types progenitor and support.
//  *      - `hair_hair` (double, default: 0.): Incremental value for 
//  *              edges between cells of types hair and hair.
//  *      - `hair_support` (double, default: 0.): Incremental value for 
//  *              edges between cells of types hair and support.
//  *      - `support_support` (double, default: 0.): Incremental value for 
//  *              edges between cells of types support and support.
//  * 
//  *  \note This increments the matrix of linetension, hence affects current and
//  *        future edges between cells of corresponding type.
//  */
// OperationBundle build_increment_linetension (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);
    
//     double incr_prog_prog(get_as<double>("progenitor_progenitor", cfg, 0.));
//     double incr_prog_hair(get_as<double>("progenitor_hair", cfg, 0.));
//     double incr_prog_supp(get_as<double>("progenitor_support", cfg, 0.));
//     double incr_prog_bnd(get_as<double>("progenitor_boundary", cfg, 0.));
//     double incr_hair_hair(get_as<double>("hair_hair", cfg, 0.));
//     double incr_hair_supp(get_as<double>("hair_support", cfg, 0.));
//     double incr_hair_bnd(get_as<double>("hair_boundary", cfg, 0.));
//     double incr_supp_supp(get_as<double>("support_support", cfg, 0.));
//     double incr_supp_bnd(get_as<double>("support_boundary", cfg, 0.));

//     if (cfg["hair_progenitor"]) {
//         double incr_hair_prog = get_as<double>("hair_progenitor", cfg);
//         if (cfg["progenitor_hair"] and incr_hair_prog != incr_prog_hair) {
//             throw KeyError("", cfg, fmt::format("Increment has to be "
//                 "symmetric matrix."));
//         }
//         incr_prog_hair = incr_hair_prog;
//     }
//     if (cfg["support_progenitor"]) {
//         double incr_supp_prog = get_as<double>("hair_progenitor", cfg);
//         if (cfg["progenitor_support"] and incr_supp_prog != incr_prog_supp) {
//             throw KeyError("", cfg, fmt::format("Increment has to be "
//                 "symmetric matrix."));
//         }
//         incr_prog_supp = incr_supp_prog;
//     }
//     if (cfg["support_hair"]) {
//         double incr_supp_hair = get_as<double>("hair_progenitor", cfg);
//         if (cfg["hair_support"] and incr_supp_hair != incr_hair_supp) {
//             throw KeyError("", cfg, fmt::format("Increment has to be "
//                 "symmetric matrix."));
//         }
//         incr_hair_supp = incr_supp_hair;
//     }

//     Operation operation = [incr_prog_prog, incr_prog_hair, incr_prog_supp,
//                            incr_hair_hair, incr_hair_supp, incr_supp_supp,
//                            incr_prog_bnd, incr_hair_bnd, incr_supp_bnd]
//             (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;

//         auto linetension = vertex_model.get_linetension();

//         linetension(CellType::progenitor,
//                     CellType::progenitor) += incr_prog_prog;
//         linetension(CellType::progenitor, CellType::hair) += incr_prog_hair;
//         linetension(CellType::progenitor, CellType::support) += incr_prog_supp;
//         linetension(CellType::progenitor,
//                     CellType::num_cell_types) += incr_prog_bnd;
        
//         linetension(CellType::hair, CellType::progenitor) += incr_prog_hair;
//         linetension(CellType::hair, CellType::hair) += incr_hair_hair;
//         linetension(CellType::hair, CellType::support) += incr_hair_supp;
//         linetension(CellType::hair,
//                     CellType::num_cell_types) += incr_hair_bnd;
        
//         linetension(CellType::support, CellType::progenitor) += incr_prog_supp;
//         linetension(CellType::support, CellType::hair) += incr_hair_supp;
//         linetension(CellType::support, CellType::support) += incr_supp_supp;
//         linetension(CellType::support,
//                     CellType::num_cell_types) += incr_supp_bnd;
                    
//         linetension(CellType::num_cell_types,
//                     CellType::progenitor) += incr_prog_bnd;
//         linetension(CellType::num_cell_types,
//                     CellType::hair) += incr_hair_bnd;
//         linetension(CellType::num_cell_types,
//                     CellType::support) += incr_supp_bnd;

//         // set the contractility and update the edge properties
//         vertex_model.set_linetension(linetension, true);
//     };

//     return std::make_pair(operation, params);
// }

// /// The operation to increment the shape index
// /** The following parameter are extracted from cfg 
//  *  (besides those passed to `OperationParams`):
//  *      - `progenitor` (double, default: 0.): Incremental value for cells of
//  *              type progenitor.
//  *      - `hair` (double, default: 0.): Incremental value for cells of
//  *              type hair.
//  *      - `support` (double, default: 0.): Incremental value for cells of
//  *              type support.
//  * 
//  *  \note This affects the current entities properties, it does not overwrite
//  *        changes in the past.
//  */
// OperationBundle build_increment_shape_index (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);
    
//     /// Incremental value for progenitor cells.
//     /** \details default: 0 . pass as 'progenitor'. */
//     double incr_prog(get_as<double>("progenitor", cfg, 0.));

//     /// Incremental value for hair cells.
//     /** \details default: 0 . pass as 'hair'. */
//     double incr_hair(get_as<double>("hair", cfg, 0.));

//     /// Incremental value for support cells.
//     /** \details default: 0 . pass as 'support'. */
//     double incr_supp(get_as<double>("support", cfg, 0.));

//     Operation operation = [incr_prog, incr_hair, incr_supp]
//             (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;

//         PCPVertex::RuleFuncCell update = [incr_prog, incr_hair, incr_supp]
//                 (const auto& cell)
//         {
//             auto state = cell->state;
//             if (state.type == CellType::progenitor) {
//                 state.shape_index_preferential += incr_prog;
//             }
//             else if (state.type == CellType::hair) {
//                 state.shape_index_preferential += incr_hair;
//             }
//             else if (state.type == CellType::support) {
//                 state.shape_index_preferential += incr_supp;
//             }
//             return state;
//         };

//         apply_rule<Update::sync>(update, vertex_model.get_am().cells());
//     };

//     return std::make_pair(operation, params);
// }

// /// Initialises a stripe that fits all vertices
// OperationBundle build_initialise_stripe_boundary(
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     Operation operation = [cfg](PCPVertex& vertex_model) {
//         vertex_model.initialise_stripe_boundary(cfg);
//     };

//     return std::make_pair(operation, params);
// }

// /// Increment the width of the stripe
// /** Stripe needs to be initialised!
//  */
// template<typename Logger>
// OperationBundle build_increment_stripe_width (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params,
//         std::shared_ptr<Logger> logger = nullptr,
//         std::function<void()> monitor = [](){ })
// {
//     OperationParams params(name, cfg, default_minim_params);

//     // the incremental change in height
//     double dH = get_as<double>("dH", cfg);

//     // the target height (if > 0)
//     auto H_target = get_as<double>("H_target", cfg, 0.);
//     bool relative_target = get_as<bool>("relative_target", cfg, false);
//     std::shared_ptr<double> H_final = std::make_shared<double>(0.);

//     std::shared_ptr<double> update_potential(nullptr);
//     if (cfg["update_potential_constant"]) {
//         update_potential = std::make_shared<double>(
//             get_as<double>("update_potential_constant", cfg));
//     }

//     MinimizationParams minimization(
//         get_as<Config>("minimization", cfg, Config()),
//                        params.minimization_params);

//     Operation operation =
//     [dH, H_target, relative_target, H_final, update_potential,
//      minimization, logger, monitor]
//     (PCPVertex& vertex_model)
//     {
//         if (update_potential) {
//             vertex_model.update_stripe_boundary_potential(*update_potential);
//         }

//         double H = vertex_model.get_stripe_boundary_width();
//         if (H > 1.e12) {
//             throw std::runtime_error("Stripe boundary needs to be initialised "
//                 "before it can be changed in `increment_stipe_width`!");
//         }
//         if (*H_final < 1.e-12) {
//             *H_final = H_target;
//             if (relative_target) {
//                 *H_final *= H;
//             }
//         }

//         // Reduce H by dH until H_final
//         if (*H_final > 1.e-12) {
//             while (H > *H_final + 1.e-12) {
//                 H = vertex_model.increment_stripe_boundary_width(
//                     std::max(dH, *H_final - H));
//                 if (H > *H_final + 1.e-12) {
//                     vertex_model.minimize_energy(minimization, monitor);
//                 }
//             }
//         }

//         // reduce by dH
//         else {
//             vertex_model.increment_stripe_boundary_width(dH);
//         }
//     };

//     return std::make_pair(operation, params);
// }

// /// Increment the width of the stripe
// /** Stripe needs to be initialised!
//  */
// OperationBundle build_increment_stripe_curvature (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double dk = get_as<double>("increment_curvature", cfg);

//     std::shared_ptr<double> update_potential(nullptr);
//     if (cfg["update_potential_constant"]) {
//         update_potential = std::make_shared<double>(
//             get_as<double>("update_potential_constant", cfg));
//     }

//     Operation operation = [dk, update_potential](PCPVertex& vertex_model)
//     {
//         double H = vertex_model.get_stripe_boundary_width();
//         if (H > 1.e12) {
//             throw std::runtime_error("Stripe boundary needs to be initialised "
//                 "before it can be changed in `increment_stipe_curvature`!");
//         }

//         if (update_potential) {
//             vertex_model.update_stripe_boundary_potential(*update_potential);
//         }

//         vertex_model.increment_stripe_boundary_curvature(dk);
//     };

//     return std::make_pair(operation, params);
// }

// template <class PlanarCellPolarity, typename Logger>
// OperationBundle build_iterate_pcp (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params,
//         std::shared_ptr<PlanarCellPolarity> pcp,
//         std::shared_ptr<bool> prolog,
//         std::shared_ptr<Logger> logger = nullptr,
//         std::function<void()> monitor = [](){ })
// {
//     if (not pcp) {
//         throw std::runtime_error("Received nullptr in build_iterate_pcp!");
//     }
    
//     OperationParams params(name, cfg, default_minim_params);

//     auto num_steps(get_as<std::size_t>("num_steps", cfg, 1));
//     auto pcp_cfg(get_as<Config>("PlanarCellPolarity", cfg, Config()));

//     Operation operation = [pcp, prolog, pcp_cfg, num_steps, monitor, logger]
//             ([[maybe_unused]] PCPVertex& vertex_model)
//     {
//         pcp->update_parameters(pcp_cfg);

//         if (not *prolog) {
//             pcp->prolog();
//             *prolog = true;
//         }

//         for (std::size_t step = 0; step < num_steps; step++) {
//             pcp->iterate();
//             monitor();

//             if (stop_now.load()) {
//                 pcp->get_logger()->warn("Was told to stop. Not iterating "
//                     "further ...");
//                 throw GotSignal(received_signum.load());
//             }
//         }
//     };

//     return std::make_pair(operation, params);
// }




// /// A pure shear experiment -- convergence and extension
// OperationBundle build_pure_shear (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double v0 = get_as<double>("shear_velocity", cfg);
//     bool deform_plastic = get_as<bool>("deform_plastic", cfg);

//     auto initialised = std::make_shared<bool>(false);
//     auto rho = std::make_shared<double>(0.);
//     auto Lx0 = std::make_shared<double>(0.);
//     auto Ly0 = std::make_shared<double>(0.);

//     Operation operation =
//         [rho, v0, Lx0, Ly0,initialised, deform_plastic]
//         (PCPVertex& vertex_model)
//     {
//         if (not vertex_model.get_space()->periodic) {
//             throw std::runtime_error("Space needs to be periodic to apply "
//                 "`pure shear` operation");
//         }

//         using SpaceVec = typename PCPVertex::SpaceVec;
//         auto space = vertex_model.get_space();

//         if (not *initialised) {
//             auto domain = space->get_domain_size();
//             *Lx0 = domain[0];
//             *Ly0 = domain[1];
//             vertex_model.get_logger()->info("Initializing the domain shape "
//                 "for pure shear operation. Calling the pure ;shear operation "
//                 "again will overwrite changes.");
            
//             *initialised = true;
//         }

//         auto domain = space->get_domain_size();
//         double Lx = *Lx0 * exp(*rho);
//         double Ly = *Ly0 * exp(-*rho);
//         if (fabs(domain[0] - Lx) > 1.e-10 or fabs(domain[1] - Ly) > 1.e-10) {
//             vertex_model.get_logger()->error("Vertex model domain: {}, {}",
//                                              domain[0], domain[1]);
//             vertex_model.get_logger()->error("Expected domain: {}, {}",
//                                              Lx, Ly);
//             throw std::runtime_error("Failed to apply `pure shear` operation, "
//                 "because the domain size was changed between two operations. "
//                 "Calling this operation would overwrite changes!");
//         }

//         *rho += v0;
//         Lx = *Lx0 * exp(*rho);
//         Ly = *Ly0 * exp(-*rho);

//         const auto& cells = vertex_model.get_am().cells();
//         double max_A0 = 0.;
//         for (const auto& cell : cells) {
//             max_A0 = std::max(max_A0, cell->state.area_preferential());
//         }
//         double L = 2 * sqrt(max_A0 / M_PI); // diameter of a cell

//         if (Lx < 3 * L or Ly < 3 * L) {
//             throw std::runtime_error("Cannot apply `pure shear` operation as "
//                 "the shorter tissue axis does not fit 3 diameters of a cell!");
//         }

//         vertex_model.stretch_domain(SpaceVec({Lx, Ly}) - domain,
//                                     false, false, false,
//                                     deform_plastic);
//         // NOTE area preserving
//     };

//     return std::make_pair(operation, params);
// }
// /// The operation to relax preferential area
// /** \details The following parameter are extracted from cfg 
//  *           (besides those passed to `OperationParams`):
//  *               - `factor` (double): incremental factor c. Should be in [0, 1].
//  *               - `factor_type` (double, default: 0.): incremental factor c'.
//  *                      Should be in [0, 1]. A factor that relaxes all cells
//  *                      of one type to the average cell area
//  *               - `factor_tissue` (double, default: 0.): incremental factor c''.
//  *                      Should be in [0, 1]. A factor that relaxes all cells
//  *                      to aim for a target average cell area.
//  *               - `average_cell_area` (double, default: 1): The target
//  *                      average cell area.
//  *               - `minimum` (double, optional): Don't set preferential area 
//  *                                               below this value.
//  *               - `minimum` (double, optional): Don't set preferential area 
//  *                                               above this value.
//  *               - `minimum_hair` (double, optional): Limit to the preferential 
//  *                          area of hair cells.
//  *               - `maximum_hair` (double, optional): Limit to the preferential 
//  *                          area of hair cells.
//  *               - `minimum_support` (double, optional): Limit to the 
//  *                          preferential area of cells.
//  *               - `maximum_support` (double, optional): Limit to the 
//  *                          preferential area of cells.
//  * 
//  *  Preferential area is relaxed towards the actual cell area:
//  *  \f$ A_0(t+1) - A_0(t) =   c (A - A_0) + c\prime (<A>_{type} - A0)
//  *                          + c\prime\prime (<A> - A_{target}) \f$
//  * 
//  *  Operation is applied to all cells with respective \f$A\f$ and \f$A_0\f$.
//  */
// OperationBundle build_relax_area (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double factor(get_as<double>("factor", cfg));
//     double factor_type(  get_as<double>("factor_type", cfg, 0.));
//     double factor_tissue(get_as<double>("factor_tissue", cfg, 0.));
//     double target_area(  get_as<double>("average_cell_area", cfg, 1.));

//     double minimum(get_as<double>("minimum", cfg, 0.));
//     double maximum(get_as<double>("maximum", cfg, 0.));

//     double minimum_HC(get_as<double>("minimum_hair", cfg, 0.));
//     double maximum_HC(get_as<double>("maximum_hair", cfg, 0.));

//     double minimum_SC(get_as<double>("minimum_support", cfg, 0.));
//     double maximum_SC(get_as<double>("maximum_support", cfg, 0.));

//     Operation operation = [factor, factor_type, factor_tissue, target_area,
//                            minimum, maximum,
//                            minimum_HC, maximum_HC,
//                            minimum_SC, maximum_SC]
//             (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;

//         const auto& am = vertex_model.get_am();
//         const auto& cells = am.cells();

//         double average_HC = 0.;
//         std::size_t num_HC = 0;
//         double average_SC = 0.;
//         std::size_t num_SC = 0;
//         double average = 0.;
//         for (const auto& cell : cells) {
//             double area = am.area_of(cell);
//             if (cell->state.type == CellType::hair) {
//                 average_HC += area;
//                 num_HC++;
//             }
//             else if (cell->state.type == CellType::support) {
//                 average_SC += area;
//                 num_SC++;
//             }
//             average += area;
//         }
//         average_HC /= static_cast<double>(num_HC); 
//         average_SC /= static_cast<double>(num_SC); 
//         average /= static_cast<double>(cells.size());

//         PCPVertex::RuleFuncCell update = 
//         [factor, factor_type, factor_tissue, target_area,
//          average_HC, average_SC, average,
//          minimum, maximum, minimum_HC, maximum_HC, minimum_SC, maximum_SC,
//          am]
//         (const auto& cell)
//         {
//             auto state = cell->state;
//             double A0 = state.area_preferential();
//             state._area_preferential += factor * (am.area_of(cell) - A0);
//             if (state.type == CellType::hair) {
//                 state._area_preferential += factor_type * (average_HC - A0);
//             }
//             else if (state.type == CellType::support) {
//                 state._area_preferential += factor_type * (average_SC - A0);
//             }
//             state._area_preferential += factor_tissue * (target_area - average);

//             // general maximum and minimum
//             state._area_preferential = std::max(state._area_preferential,
//                                                 minimum);
//             if (maximum > 1.e-12) {
//                 state._area_preferential = std::min(state._area_preferential,
//                                                     maximum);
//             }

//             // hair cell maximum and minimum 
//             if (state.type == CellType::hair) {
//                 state._area_preferential = std::max(state._area_preferential,
//                                                     minimum_HC);
//                 if (maximum_HC > 1.e-12) {
//                     state._area_preferential = 
//                         std::min(state._area_preferential, maximum_HC);
//                 }
//             }
//             // hair cell maximum and minimum 
//             else if (state.type == CellType::support) {
//                 state._area_preferential = std::max(state._area_preferential,
//                                                     minimum_SC);
//                 if (maximum_SC > 1.e-12) {
//                     state._area_preferential = 
//                         std::min(state._area_preferential, maximum_SC);
//                 }
//             }

//             return state;
//         };
        
//         apply_rule<Update::sync>(update, cells);

//         return;
//     };

//     return std::make_pair(operation, params);
// }

// /// The operation to set preferential area
// /** \details The following parameter are extracted from cfg 
//  *           (besides those passed to `OperationParams`):
//  *               - `progenitor` (double, default: 0.): new value for
//  *                      cells of type progenitor. Value 0. is ignored and
//  *                      previous value of preferential area kept.
//  *               - `stddev_progenitor` (double, default: 0.): stddev for cells
//  *                      of type progenitor; using lognormal distribution
//  *               - `hair` (double, default: 0.): new value for cells
//  *                      of type hair. Value 0. is ignored and
//  *                      previous value of preferential area kept.
//  *               - `stddev_hair` (double, default: 0.): stddev for cells
//  *                      of type hair; using lognormal distribution
//  *               - `support` (double, default: 0.): new value for cells
//  *                      of type support. Value 0. is ignored and
//  *                      previous value of preferential area kept.
//  *               - `stddev_support` (double, default: 0.): stddev for cells
//  *                      of type support; using lognormal distribution
//  *               - `relax_domain` (bool, default: false): If true, the domain
//  *                      size is adapted to fit the cells as per preferential 
//  *                      area.
//  *               - `relax_domain_PD_axis` (bool, default: false): If true,
//  *                      the domain size is adapted to fit the cells as per
//  *                      preferential area. Domain size changes only in PD (x)
//  *                      axis.
//  */
// OperationBundle build_set_area (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double prog(get_as<double>("progenitor", cfg, 0.));
//     double stddev_prog(get_as<double>("stddev_progenitor", cfg, 0.));
//     double hair(get_as<double>("hair", cfg, 0.));
//     double stddev_hair(get_as<double>("stddev_hair", cfg, 0.));
//     double support(get_as<double>("support", cfg, 0.));
//     double stddev_support(get_as<double>("stddev_support", cfg, 0.));
//     bool relax_domain(get_as<bool>("relax_domain", cfg, false));
//     bool relax_domain_PD_axis(get_as<bool>("relax_domain_PD_axis",
//                                            cfg, false));

//     if (relax_domain and relax_domain_PD_axis) {
//         throw std::invalid_argument("In operation `set area` only one of cfgs "
//             "`relax_domain` and `relax_domain_PD_axis` can be true!");
//     }

//     Operation operation = [prog, stddev_prog,
//                            hair, stddev_hair,
//                            support, stddev_support,
//                            relax_domain, relax_domain_PD_axis]
//             (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;

//         const auto& am = vertex_model.get_am();
//         const auto& cells = am.cells();

//         std::lognormal_distribution<> dist_prog, dist_hair, dist_support;
//         if (prog > 1.e-12 and stddev_prog > 1.e-12) {
//             dist_prog = get_lognormal_distribution(prog, stddev_prog);
//         }
//         if (hair > 1.e-12 and stddev_hair > 1.e-12) {
//             dist_hair = get_lognormal_distribution(hair, stddev_hair);
//         }
//         if (support > 1.e-12 and stddev_support > 1.e-12) {
//             dist_support = get_lognormal_distribution(support, stddev_support);
//         }

//         PCPVertex::RuleFuncCell update = [vertex_model,
//                                           dist_prog{std::move(dist_prog)},
//                                           dist_hair{std::move(dist_hair)},
//                                           dist_support{std::move(dist_support)},
//                                           prog, stddev_prog,
//                                           hair, stddev_hair,
//                                           support, stddev_support]
//                 (const auto& cell) mutable
//         {
//             std::function<double(double, double,
//                                  std::lognormal_distribution<>&)> new_area = 
//                 [vertex_model](double mean, double stddev,
//                                std::lognormal_distribution<>& distr)
//             {
//                 if (stddev > 1.e-12) {
//                     return distr(*vertex_model.get_rng());
//                 }
//                 else {
//                     return mean;
//                 }
//             };

//             auto state = cell->state;
//             double new_A = 0.;
//             double mean = 0.;
//             if (state.type == CellType::progenitor and prog > 1.e-12) {
//                 mean = prog;
//                 new_A = new_area(prog, stddev_prog, dist_prog);
//             }
//             else if (state.type == CellType::hair and hair > 1.e-12) {
//                 mean = hair;
//                 new_A = new_area(hair, stddev_hair, dist_hair);
//             }
//             else if (state.type == CellType::support and support > 1.e-12) {
//                 mean = support;
//                 new_A = new_area(support, stddev_support, dist_support);
//             }

//             if (mean > 1.e-12) {
//                 state._area_preferential = mean;
//                 state._area_preferential_fluctuations = new_A - mean;
//             }

//             return state;
//         };
        
//         apply_rule<Update::sync>(update, cells);

//         if (relax_domain) {
//             double area = std::accumulate(cells.begin(), cells.end(), 0.,
//                             [](const double& val, const auto& cell) {
//                                 return val + cell->state.area_preferential();
//                             });

//             PCPVertex::SpaceVec domain_size = 
//                     vertex_model.get_space()->get_domain_size();
//             vertex_model.increase_domain_size(
//                     area - domain_size[0]*domain_size[1]
//             );
//         }
//         else if (relax_domain_PD_axis) {
//             double area = std::accumulate(cells.begin(), cells.end(), 0.,
//                             [](const double& val, const auto& cell) {
//                                 return val + cell->state.area_preferential();
//                             });

//             PCPVertex::SpaceVec domain_size = 
//                     vertex_model.get_space()->get_domain_size();
//             double dA = area - domain_size[0]*domain_size[1];

//             vertex_model.stretch_domain({dA/domain_size[1], 0.},
//                                         false, false, false);
//         }

//         return;
//     };

//     return std::make_pair(operation, params);
// }

// OperationBundle build_update_boundary_parameter (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     auto boundary_cfg = get_as<Config>("boundary_parameter", cfg);

//     Operation operation = [boundary_cfg] (PCPVertex& vertex_model)
//     {
//         vertex_model.update_boundary_parameter(boundary_cfg);
//     };

//     return std::make_pair(operation, params);
// }

// /// Set a torque to cells
// /** The configuration is passed on to the constructor of RotationCellState
//  *  The config following the keys to the respective cell type:
//  *      - `hair`
//  *      - `support`
//  *      - `progenitor`
//  * 
//  *  If a configuration is not given, torque zero is set.
//  */
// OperationBundle build_set_torque (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     Config cfg_default;
//     cfg_default["torque"] = 0.;

//     Config cfg_hair = get_as<Config>("hair", cfg, cfg_default);
//     Config cfg_support = get_as<Config>("support", cfg, cfg_default);
//     Config cfg_progenitor = get_as<Config>("progenitor", cfg, cfg_default);

//     Operation operation = [cfg_hair, cfg_support, cfg_progenitor]
//                           (PCPVertex& vertex_model)
//     {
//         using CellType = PCPVertex::CellType;

//         const auto& cells = vertex_model.get_am().cells();
//         for (const auto& cell : cells) {
//             if (cell->state.type == CellType::hair) {
//                 auto rot_state = std::make_shared<RotationCellState>(cfg_hair);
//                 cell->custom_links().rotation_state = rot_state;
//             }
//             else if (cell->state.type == CellType::support) {
//                 auto rot_state = std::make_shared<RotationCellState>(
//                                                 cfg_support);
//                 cell->custom_links().rotation_state = rot_state;
//             }
//             else if (cell->state.type == CellType::progenitor) {
//                 auto rot_state = std::make_shared<RotationCellState>(
//                                                 cfg_progenitor);
//                 cell->custom_links().rotation_state = rot_state;
//             }
//         }
//     };

//     return std::make_pair(operation, params);
// }


// /// A bending model in periodic boundary conditions
// OperationBundle build_bending_box_bc (
//         std::string name, const Config& cfg,
//         const MinimizationParams& default_minim_params)
// {
//     OperationParams params(name, cfg, default_minim_params);

//     double v0 = get_as<double>("increment_curvature", cfg);
//     bool deform_plastic = get_as<bool>("deform_plastic", cfg);

//     Operation operation = [v0, deform_plastic](PCPVertex& vertex_model)
//     {
//         vertex_model.curve_boundary(v0, deform_plastic);
//     };

//     return std::make_pair(operation, params);
// }
