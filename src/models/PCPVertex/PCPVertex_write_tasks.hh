#ifndef UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH
#define UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"
#include <utopia/data_io/data_manager/defaults.hh>
#include <complex>
#include <cmath>

using namespace Utopia::DataIO;

/** Available datatree:
 *      - Energy
 *          - Energy_total
 *          - Energy_time (the linkded time for all energy adaptors)
 *          - Energy_linetension
 *          - Energy_areaelasticity
 *          - Energy_cell_contractility
 *          - Energy_edge_contractility
 *          - Energy_cell_cell_polarity
 *          - Energy_polarity_exclusion
 *          - Energy_lagrange_net_polarisation
 *          - Energy_lagrange_const_concentration
 *      - Vertices (time series groups)
 *      - Cells (time series groups)
 *      - Edges (time series groups)
 *      - Statistics
 *          - num T1s
 *          - num_T1s_attempted
 *          - num T2s
 *          - Cell_area
 *          - Statistics_time
 */
namespace Utopia::Models::PCPVertex::DataIO{

/// Datamanager adaptor for total energy
/** PCPVertex::get_energy() normalized to number of cells
 * 
 *  \note excludes boundary energy
 */
auto energy_adaptor = std::make_tuple(

    // name of the task
    "Energy_total",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(
            (model.get_energy() - model.get_boundary_energy())
            / model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Total");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end energy_adaptor

/// Datamanager adaptor for timepoints
/** The dataset to which the other energy adaptors link their coordinate time
 */
auto time_energy_adaptor = std::make_tuple(

    // name of the task
    "Energy_time",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Time");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
    }
); // end time_energy_adaptor

/// Datamanager adaptor for linetension energy
/** PCPVertex::get_energy_linetension() normalized to number of cells
 */
auto linetension_adaptor = std::make_tuple(

    // name of the task
    "Energy_linetension",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_energy_linetension() /
                       model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Linetension");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end linetension_adaptor

/// Datamanager adaptor for area elasticity energy
/** PCPVertex::get_energy_areaelasticity() normalized to number of cells
 */
auto areaelasticity_adaptor = std::make_tuple(

    // name of the task
    "Energy_areaelasticity",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_energy_areaelasticity() /
                       model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Areaelasticity");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end areaelasticity_adaptor

/// Datamanager adaptor for cell contractility energy
/** PCPVertex::get_energy_cell_contractility() normalized to number of cells
 */
auto cell_contractility_adaptor = std::make_tuple(

    // name of the task
    "Energy_cell_contractility",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(  model.get_energy_cell_contractility()
                       / model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_contractility");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end cell_contractility_adaptor

/// Datamanager adaptor for edge contractility energy
/** PCPVertex::get_energy_edge_contractility() normalized to number of cells
 */
auto edge_contractility_adaptor = std::make_tuple(

    // name of the task
    "Energy_edge_contractility",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(  model.get_energy_edge_contractility()
                       / model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Edge_contractility");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end edge_contractility_adaptor

/// Datamanager adaptor for the boundary area elasticity energy
/** PCPVertex::get_boundary_area_energy() normalized to number of cells
 */
auto boundary_area_elasticity_adaptor = std::make_tuple(

    // name of the task
    "Energy_boundary_area_elasticity",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_boundary_area_energy() /
                       model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Boundary_area_elasticity");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end areaelasticity_adaptor

/// Datamanager adaptor for boundary shape elasticity energy
/** PCPVertex::get_boundary_shape_energy() normalized to number of cells
 */
auto boundary_shape_elasticity_adaptor = std::make_tuple(

    // name of the task
    "Energy_boundary_shape_elasticity",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(  model.get_boundary_shape_energy()
                       / model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Boundary_shape_elasticity");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end cell_contractility_adaptor

/// Datamanager adaptor for cell-cell polarity energy
auto cell_cell_polarity_adaptor = std::make_tuple(

    // name of the task
    "Energy_cell_cell_polarity",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    []([[maybe_unused]] auto& dataset, [[maybe_unused]] auto& model) {
        // dataset->write(model.get_energy_cell_cell_polarity() / 
        //                model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_cell_polarity");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end cell_cell_polarity_adaptor

/// Datamanager adaptor for polarity exclusion energy
auto polarity_exclusion_adaptor = std::make_tuple(

    // name of the task
    "Energy_polarity_exclusion",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    []([[maybe_unused]] auto& dataset, [[maybe_unused]] auto& model) {
        // dataset->write(model.get_energy_polarity_exclusion() / 
        //                model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Polarity_exclusion");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end polarity_exclusion_adaptor

/// Datamanager adaptor for lagrange net polarization energy
auto lagrange_net_polarisation_adaptor = std::make_tuple(

    // name of the task
    "Energy_lagrange_net_polarisation",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    []([[maybe_unused]] auto& dataset, [[maybe_unused]] auto& model) {
        // dataset->write(model.get_energy_lagrange_net_polarisation() / 
        //                model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Lagrange_net_polarisation");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end lagrange_net_polarisation_adaptor

/// Datamanager adaptor for lagrange const concentration energy
auto lagrange_const_concentration_adaptor = std::make_tuple(

    // name of the task
    "Energy_lagrange_const_concentration",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    []([[maybe_unused]] auto& dataset, [[maybe_unused]] auto& model) {
        // dataset->write(model.get_energy_lagrange_const_concentration() / 
        //                model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Lagrange_const_concentration");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end lagrange_const_concentration_adaptor


/// Datamanager adaptor for vertex properties
/** \details Properties are
 *      -# absolute x coordinate
 *      -# absolute y coordinate
 * 
 *  Attributes are:
 *      -# Lx: Domain size in x coordinate
 *      -# Ly: Domain size in y coordinate
 */
template <typename SpaceVec>
auto vertices_adaptor = std::make_tuple(
    // name of the task
    "Vertices",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Vertices");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& vertices = am.vertices();
        dataset->write(vertices.begin(), vertices.end(),
                        [am](auto&& vertex) { 
                            return am.position_of(vertex)[0]; });
        dataset->write(vertices.begin(), vertices.end(),
                        [am](auto&& vertex) { 
                            return am.position_of(vertex)[1]; });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_am().vertices().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series"); },

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute("coords__property", 
                                  std::vector<std::string>({"x", "y"}));

        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& vertices = model.get_am().vertices();
        std::vector<std::size_t> ids{};
        ids.reserve(vertices.size());
        std::transform(vertices.begin(), vertices.end(),
                       std::back_inserter(ids),
                       [](const auto& v) { return v->id(); });
        hdfdataset->add_attribute("coords__id", ids);

        const SpaceVec domain = model.get_space()->get_domain_size();
        hdfdataset->add_attribute("Lx", domain[0]);
        hdfdataset->add_attribute("Ly", domain[1]);
    }
    
); // end vertex position adaptor

/// Datamanager adaptor for cell properties
/** \details Properties are
 *      -# absolute position cell_type
 *      -# absolute position x coordinate of cell center
 *      -# absolute position y coordinate of cell center
 *      -# absolute position area
 *      -# absolute position perimeter
 *      -# absolute position shape_index: p = perimeter / sqrt(area)
 *      -# absolute position number of neighbors
 * 
 *  Attributes are:
 *      -# Lx: Domain size in x coordinate
 *      -# Ly: Domain size in y coordinate
 */
template <typename SpaceVec, typename CellType>
auto cells_adaptor = std::make_tuple(
    // name of the task
    "Cells",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Cells");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        dataset->write(cells.begin(), cells.end(),
                       [](const auto& cell) {
                            return static_cast<double>(cell->state.type);
                       });
        
        std::vector<SpaceVec> centers;
        centers.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(centers),
                       [am](const auto& c) { return am.barycenter_of(c); });
        dataset->write(centers.begin(), centers.end(),
                       [](auto&& pos) { return pos[0]; });
        dataset->write(centers.begin(), centers.end(),
                       [](auto&& pos) { return pos[1]; });

        std::vector<double> areas;
        areas.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(areas),
                       [am](const auto& c) { return am.area_of(c); });
        dataset->write(areas);

        std::vector<double> area_preferentials;
        area_preferentials.reserve(cells.size());
        std::transform(cells.begin(), cells.end(),
                       std::back_inserter(area_preferentials),
                       [](const auto& c) { 
                            return c->state.area_preferential; });
        dataset->write(area_preferentials);
        
        std::vector<double> perimeters;
        perimeters.reserve(cells.size());
        std::transform(cells.begin(), cells.end(),
                       std::back_inserter(perimeters),
                       [am](const auto& c) { return am.perimeter_of(c); });
        dataset->write(perimeters);
        
        std::vector<double> shape_indices;
        shape_indices.reserve(cells.size());
        for (unsigned int i = 0; i < cells.size(); i++) {
            shape_indices.push_back(perimeters[i]/sqrt(areas[i]));
        }
        dataset->write(shape_indices);

        dataset->write(cells.begin(), cells.end(),
                       [am](const auto& cell) {
                            return static_cast<double>(
                                   am.neighbors_of(cell).size());
                       });
        dataset->write(cells.begin(), cells.end(),
                       [am](const auto& cell) {
                            unsigned int num_hair_neighbors = 0;
                            for (const auto n : am.neighbors_of(cell)) {
                                num_hair_neighbors += 
                                    (n->state.type == CellType::hair);
                            }
                            return static_cast<double>(num_hair_neighbors);
                       });

        dataset->write(cells.begin(), cells.end(),
                       [am](const auto& cell) {
                            if (cell->state.type != CellType::hair) {
                                return 0.;
                            }

                            auto hair_neighbors = am.hair_neighbors_of(cell);
                            if (hair_neighbors.size() < 4) {
                                return 0.;
                            }

                            using namespace std::complex_literals;
                            std::complex<double> hex_order = std::accumulate(
                                hair_neighbors.begin(), hair_neighbors.end(),
                                std::complex<double>(0., 0.),
                                [am, cell](std::complex<double> val,
                                           const auto& nb)
                                {
                                    using namespace std::complex_literals;

                                    double distance = am.distance(cell, nb);
                                    double dx = am.displacement(cell, nb)[0];
                                    double theta = acos(dx / distance);
                                    return val + std::exp(1i * 6. * theta);
                                }
                            );

                            std::complex<double> N(hair_neighbors.size());                            
                            return std::norm(hex_order / N);
                       });

        dataset->write(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                if (not cell->custom_links().rotation_state) {
                    return 0.;
                }
                return cell->custom_links().rotation_state->tracked_rotation;
            });
        
        dataset->write(
            cells.begin(), cells.end(),
            [am](const auto& cell) {
                return static_cast<double>(am.is_boundary(cell));
            });

        // the polarity
        // std::vector<SpaceVec> polarities;
        // for (const auto& c : cells) {
        //     SpaceVec polarity(arma::fill::zeros);
        //     for (auto [e, flip] : c->custom_links().edges) {
        //         auto a = e->a->position(), b = e->b->position();
        //         if (flip) {
        //             std::swap(a, b);
        //         }
        //         auto disp = model.get_space()->displacement(a, b);

        //         double sigma = e->state.get_sigma(c);
                
        //         polarity += sigma * disp;
        //     }
        //     polarities.push_back(std::make_pair(polarity_x, polarity_y));
        // }
        // dataset->write(polarities.begin(), polarities.end(),
        //                [](auto&& pos) { return pos[0]; });
        // dataset->write(polarities.begin(), polarities.end(),
        //                [](auto&& pos) { return pos[1]; });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {12, m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute("coords__property", 
                std::vector<std::string>({
                    "cell_type",
                    "x",
                    "y",
                    "area",
                    "area_preferential",
                    "perimeter",
                    "shape_index",
                    "num_neighbors",
                    "num_hair_neighbors",
                    "hexatic_order",
                    "rotation",
                    "is_boundary"
                }));
                    // "polarity_x",
                    // "polarity_y"}));
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& cells = model.get_am().cells();
        std::vector<std::size_t> ids{};
        ids.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(ids),
                       [](const auto& c) { return c->id(); });
        hdfdataset->add_attribute("coords__id", ids);

        const SpaceVec domain = model.get_space()->get_domain_size();
        hdfdataset->add_attribute("Lx", domain[0]);
        hdfdataset->add_attribute("Ly", domain[1]);
    }    
); // end cell position adaptor

/// Datamanager adaptor for edges properties
auto edges_adaptor = std::make_tuple(

    // name of the task
    "Edges",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Edges");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& edges = model.get_am().edges();
        dataset->write(edges.begin(), edges.end(),
                       [](const auto& edge) {
                           return static_cast<int>(
                                    edge->custom_links().a->id());
                       });
        dataset->write(edges.begin(), edges.end(),
                       [](const auto& edge) {
                           return static_cast<int>(
                                    edge->custom_links().b->id());
                       });
    },
                
    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_am().edges().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute("coords__property", 
                                  std::vector<std::string>({"vertex_a",
                                                            "vertex_b"}));
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& edges = model.get_am().edges();
        std::vector<std::size_t> ids{};
        ids.reserve(edges.size());
        std::transform(edges.begin(), edges.end(), std::back_inserter(ids),
                       [](const auto& e) { return e->id(); });
        hdfdataset->add_attribute("coords__id", ids);
    }     
); // end edge link adaptor

/// Datamanager adaptor for cell energies
/** Energies are
 *      -# PCPVertex::get_energy_areaelasticity
 *      -# PCPVertex::get_energy_cell_contractility
 */
auto cell_energies_adaptor = std::make_tuple(
    // name of the task
    "Cell_energies",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Cell_energies");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();

        dataset->write(cells.begin(), cells.end(),
                       [model](const auto& c) {
                            return model.get_energy_areaelasticity({c});
                       });
        dataset->write(cells.begin(), cells.end(),
                       [model](const auto& c) {
                            return model.get_energy_cell_contractility({c});
                       });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "energy_term");
        hdfdataset->add_attribute("coords__energy_term", 
                std::vector<std::string>({
                    "area_elasticity",
                    "cell_contractility"
                }));
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& cells = model.get_am().cells();
        std::vector<std::size_t> ids{};
        ids.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(ids),
                       [](const auto& c) { return c->id(); });
        hdfdataset->add_attribute("coords__id", ids);
    }    
); // end cell energy adaptor

/// Datamanager adaptor for edge energies
/** Energies are
 *      -# PCPVertex::get_energy_linetension
 *      -# PCPVertex::get_energy_edge_contractility
 */
auto edge_energies_adaptor = std::make_tuple(

    // name of the task
    "Edge_energies",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Edge_energies");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& edges = am.edges();

        dataset->write(edges.begin(), edges.end(),
                       [model](const auto& e) {
                            return model.get_energy_linetension({e});
                       });
        dataset->write(edges.begin(), edges.end(),
                       [model](const auto& e) {
                            return model.get_energy_edge_contractility({e});
                       });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_am().edges().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "energy_term");
        hdfdataset->add_attribute("coords__energy_term", 
                std::vector<std::string>({
                    "linetension",
                    "edge_contractility"
                }));
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& edges = model.get_am().edges();
        std::vector<std::size_t> ids{};
        ids.reserve(edges.size());
        std::transform(edges.begin(), edges.end(), std::back_inserter(ids),
                       [](const auto& e) { return e->id(); });
        hdfdataset->add_attribute("coords__id", ids);
    }    
); // end edge energy adaptor

/// Datamanager adaptor for the cluster of hair cells
/** Attributes are:
 *      -# Lx: Domain size in x coordinate
 *      -# Ly: Domain size in y coordinate
 */
auto hair_cluster_adaptor = std::make_tuple(
    // name of the task
    "Hair_cluster",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Hair_cluster");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& cells = model.get_am().cells();
        dataset->write(cells.begin(), cells.end(),
            [](const auto& cell) {
                std::size_t cluster_id;
                if (cell->custom_links().nd_cell) {
                    cluster_id = cell->custom_links().nd_cell->state.cluster_id;
                }
                else {
                    cluster_id = 0;
                }
                return cluster_id;
            });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "id");
    }    
); // end hair cluster adaptor

/// Datamanager adaptor for T1 cell intercalation counter
auto T1_adaptor = std::make_tuple(

    // name of the task
    "num_T1s",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_num_T1s());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("num_T1s");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end T1_adaptor

/// Datamanager adaptor for attempted T1 cell intercalation counter
auto T1_attempted_adaptor = std::make_tuple(

    // name of the task
    "num_T1s_attempted",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_num_T1s_attempted());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("num_T1s_attempted");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end T1_attempted_adaptor

/// Datamanager adaptor for T2 cell extrusion counter
auto T2_adaptor = std::make_tuple(

    // name of the task
    "num_T2s",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_num_T2s());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("num_T2s");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end T2_adaptor

/// Datamanager adaptor for cell area statistics
template <typename CellType>
auto cell_area_adaptor = std::make_tuple(

    // name of the task
    "Cell_area",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        std::vector<double> area_cells(CellType::num_cell_types, 0.);
        std::vector<int> num_cells(CellType::num_cell_types, 0.);
        for (auto c : cells) {
            num_cells[c->state.type]++;
            area_cells[c->state.type] += am.area_of(c);
        }
        double area_total = std::accumulate(area_cells.begin(),
                                            area_cells.end(), 0.);
        double area_average = area_total / cells.size();

        for (int i = 0; i < CellType::num_cell_types; i++) {
            area_cells[i] /= num_cells[i];
        }
        double average_hair = area_cells[CellType::hair];
        double average_support = area_cells[CellType::support];

        std::vector<double> data = {area_average,
                                    average_hair, average_support};
        dataset->write(data);
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_area", {H5S_UNLIMITED, 3});
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        hdfdataset->add_attribute("dim_name__1", "property");
        hdfdataset->add_attribute("coords__property", 
            std::vector<std::string>({"area_average", "area_hair_average",
                                      "area_support_average"}));

    }
); // end cell_area_adaptor

/// Datamanager adaptor for timepoints
auto statistics_time_adaptor = std::make_tuple(

    // name of the task
    "Statistics_time",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_time());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Time");
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
    }
); // end statistics_time_adaptor



/// The Datamanager decider on equilibrium condition
template<typename VertexModel>
struct DeciderEquilibriumCondition
:
    Utopia::DataIO::Default::Decider<VertexModel>
{
    bool operator()(VertexModel& m) override {
        return m.equilibrium_condition();
    }

    void set_from_cfg(const Config&) override { }
};

/// Create a Decidermap of custom deciders, including defaults
template <typename VertexModel>
Utopia::DataIO::Default::DefaultDecidermap< VertexModel >
build_custom_deciders ()
{
    auto deciders = Utopia::DataIO::Default::default_deciders< VertexModel >;

    // add custom deciders
    deciders["equilibrium_decider"] =
    []() -> decltype(auto) {
        return std::make_shared<
                    DataIO::DeciderEquilibriumCondition< VertexModel >>();
    };

    return deciders;
}

} // namespace Utopia::Models::PCPVertex::DataIO


#endif // UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH