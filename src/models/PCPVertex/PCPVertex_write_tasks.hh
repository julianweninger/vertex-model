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
 *          - transitions
 *      - Vertices (time series groups)
 *      - Cells (time series groups)
 *      - Edges (time series groups)
 *      - Statistics
 *          - Statistics_time
 *          - Cell_stats
 *          - Hair_cell_stats
 *          - Support_cell_stats
 *          - Bulk_cell_stats
 *          - Bulk_hair_cell_stats
 *          - Bulkd_support_cell_stats
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
                            return c->state.area_preferential(); });
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
                           return am.hexatic_order_of(cell);
                       });

        dataset->write(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                if (not cell->custom_links().rotation_state) {
                    return 0.;
                }
                return cell->custom_links().rotation_state->tracked_rotation;
            });

        dataset->write(cells.begin(), cells.end(),
                       [](const auto& cell) {
                           return cell->state.polarity;
                       });
        
        dataset->write(
            cells.begin(), cells.end(),
            [am](const auto& cell) {
                return static_cast<double>(am.is_boundary(cell));
            });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {13, m.get_am().cells().size()});
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
                    "polarity",
                    "is_boundary"
                }));
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

auto transition_adaptor = std::make_tuple(

    // name of the task
    "transitions",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        std::vector<double> stats{};
        stats.reserve(3);
        
        stats.push_back(model.get_num_T1s());
        stats.push_back(model.get_num_T1s_attempted());
        stats.push_back(model.get_num_T2s());
        
        dataset->write(stats);
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("transitions", { H5S_UNLIMITED, 3 });
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
                std::vector<std::string>({
                    "num_T1s",
                    "num_T1s_attempted",
                    "num_T2s"
                }));
    }

); // transitions

template <typename CellContainer, typename AgentManager>
std::vector<double> generate_statistics (const CellContainer& cells,
                                         const AgentManager& am)
{
    std::function<double(const std::vector<double>&)> average =
    [](const std::vector<double>& values) {
        if (values.size() == 0) {
            return std::nan("1");
        }
        return (  std::accumulate(values.begin(), values.end(), 0.)
                / static_cast<double>(values.size()));
    };
    std::function<double(const std::vector<double>&, double)> stddev =
    [](const std::vector<double>& values, double mean) {
        if (values.size() == 0) {
            return std::nan("1");
        }
        return sqrt(  std::accumulate(values.begin(), values.end(), 0.,
                            [mean](const double& sum, const double& val) {
                                return (  sum
                                        + std::pow(val - mean, 2));
                            })
                    / static_cast<double>(values.size()));
    };
    std::function<double(const std::vector<double>&)> max =
    [](const std::vector<double>& values) {
        if (values.size() == 0) {
            return std::nan("1");
        }
        return *std::max_element(values.begin(), values.end());
    };
    std::function<double(const std::vector<double>&)> min =
    [](const std::vector<double>& values) {
        if (values.size() == 0) {
            return std::nan("1");
        }
        return *std::min_element(values.begin(), values.end());
    };

    std::vector<double> values;
    values.reserve(cells.size());

    std::vector<double> stats;
    stats.reserve(16);

    // area
    std::transform(cells.begin(), cells.end(), std::back_inserter(values),
                    [am](const auto& cell) {
                        return am.area_of(cell);
                    });
    stats.push_back(average(values));
    stats.push_back(stddev(values, stats.back()));
    stats.push_back(max(values));
    stats.push_back(min(values));
    values.clear();
    values.reserve(cells.size());

    // area preferential
    std::transform(cells.begin(), cells.end(), std::back_inserter(values),
                    [](const auto& cell) {
                        return cell->state.area_preferential();
                    });
    stats.push_back(average(values));
    stats.push_back(stddev(values, stats.back()));
    stats.push_back(min(values));
    stats.push_back(max(values));
    values.clear();
    values.reserve(cells.size());

    // perimeter
    std::transform(cells.begin(), cells.end(), std::back_inserter(values),
                    [am](const auto& cell) {
                        return am.perimeter_of(cell);
                    });
    stats.push_back(average(values));
    stats.push_back(stddev(values, stats.back()));
    stats.push_back(max(values));
    stats.push_back(min(values));
    values.clear();
    values.reserve(cells.size());

    // shape index
    std::transform(cells.begin(), cells.end(), std::back_inserter(values),
                    [am](const auto& cell) {
                        return am.shape_index_of(cell);
                    });
    stats.push_back(average(values));
    stats.push_back(stddev(values, stats.back()));
    stats.push_back(max(values));
    stats.push_back(min(values));
    values.clear();
    values.reserve(cells.size());

    // num neighbors
    std::transform(cells.begin(), cells.end(), std::back_inserter(values),
                    [am](const auto& cell) {
                        return am.neighbors_of(cell).size();
                    });
    stats.push_back(average(values));
    stats.push_back(stddev(values, stats.back()));
    stats.push_back(max(values));
    stats.push_back(min(values));
    values.clear();
    values.reserve(cells.size());

    // num_hair_neighbors
    std::transform(cells.begin(), cells.end(), std::back_inserter(values),
                    [am](const auto& cell) {
                        return am.hair_neighbors_of(cell).size();
                    });
    stats.push_back(average(values));
    stats.push_back(stddev(values, stats.back()));
    stats.push_back(max(values));
    stats.push_back(min(values));
    values.clear();
    values.reserve(cells.size());

    // hexatic order
    std::transform(cells.begin(), cells.end(), std::back_inserter(values),
                    [am](const auto& cell) {
                        return am.hexatic_order_of(cell);
                    });
    stats.push_back(average(values));
    stats.push_back(stddev(values, stats.back()));
    stats.push_back(max(values));
    stats.push_back(min(values));
    values.clear();
    values.reserve(cells.size());

    // rotation
    std::transform(cells.begin(), cells.end(), std::back_inserter(values),
                    [am](const auto& cell) {
                        if (not cell->custom_links().rotation_state) {
                            return 0.;
                        }
                        return cell->custom_links().rotation_state->\
                                        tracked_rotation;
                    });
    stats.push_back(average(values));
    stats.push_back(stddev(values, stats.back()));
    stats.push_back(max(values));
    stats.push_back(min(values));
    
    values.clear();
    values.reserve(cells.size());

    stats.push_back(cells.size());

    return stats;
}

/// Datamanager adaptor for cell area statistics
template <typename Cell, typename CellType = typename Cell::State::CellType>
auto cell_stats_adaptor = std::make_tuple(

    // name of the task
    "Cell_stats",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        
        dataset->write(generate_statistics(cells, am));
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_stats", {H5S_UNLIMITED, 33});
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
            std::vector<std::string>({"area",
                                      "area__stddev",
                                      "area__max",
                                      "area__min",
                                      "area_preferential",
                                      "area_preferential__stddev",
                                      "area_preferential__max",
                                      "area_preferential__min",
                                      "perimeter",
                                      "perimeter__stddev",
                                      "perimeter__max",
                                      "perimeter__min",
                                      "shape_index",
                                      "shape_index__stddev",
                                      "shape_index__max",
                                      "shape_index__min",
                                      "num_neighbors",
                                      "num_neighbors__stddev",
                                      "num_neighbors__max",
                                      "num_neighbors__min",
                                      "num_hair_neighbors",
                                      "num_hair_neighbors__stddev",
                                      "num_hair_neighbors__max",
                                      "num_hair_neighbors__min",
                                      "hexatic_order",
                                      "hexatic_order__stddev",
                                      "hexatic_order__max",
                                      "hexatic_order__min",
                                      "rotation",
                                      "rotation__stddev",
                                      "rotation__max",
                                      "rotation__min",
                                      "count"
                                      }));


    }
); // end cell_stats_adaptor

/// Datamanager adaptor for cell area statistics
template <typename Cell, typename CellType = typename Cell::State::CellType>
auto hair_cell_stats_adaptor = std::make_tuple(

    // name of the task
    "Hair_cell_stats",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        
        AgentContainer<Cell> hair_cells;
        hair_cells.reserve(cells.size());
        std::copy_if(cells.begin(), cells.end(), std::back_inserter(hair_cells),
                     [](const auto& cell) {
                         return cell->state.type == CellType::hair;
                     });
        hair_cells.shrink_to_fit();
        
        dataset->write(generate_statistics(hair_cells, am));
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Hair_cell_stats", {H5S_UNLIMITED, 33});
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
            std::vector<std::string>({"area",
                                      "area__stddev",
                                      "area__max",
                                      "area__min",
                                      "area_preferential",
                                      "area_preferential__stddev",
                                      "area_preferential__max",
                                      "area_preferential__min",
                                      "perimeter",
                                      "perimeter__stddev",
                                      "perimeter__max",
                                      "perimeter__min",
                                      "shape_index",
                                      "shape_index__stddev",
                                      "shape_index__max",
                                      "shape_index__min",
                                      "num_neighbors",
                                      "num_neighbors__stddev",
                                      "num_neighbors__max",
                                      "num_neighbors__min",
                                      "num_hair_neighbors",
                                      "num_hair_neighbors__stddev",
                                      "num_hair_neighbors__max",
                                      "num_hair_neighbors__min",
                                      "hexatic_order",
                                      "hexatic_order__stddev",
                                      "hexatic_order__max",
                                      "hexatic_order__min",
                                      "rotation",
                                      "rotation__stddev",
                                      "rotation__max",
                                      "rotation__min",
                                      "count"
                                      }));


    }
); // end hair_cell_stats_adaptor

/// Datamanager adaptor for cell area statistics
template <typename Cell, typename CellType = typename Cell::State::CellType>
auto support_cell_stats_adaptor = std::make_tuple(

    // name of the task
    "Support_cell_stats",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        
        AgentContainer<Cell> support_cells;
        support_cells.reserve(cells.size());
        std::copy_if(cells.begin(), cells.end(),
                     std::back_inserter(support_cells),
                     [](const auto& cell) {
                         return cell->state.type == CellType::support;
                     });
        support_cells.shrink_to_fit();
        
        dataset->write(generate_statistics(support_cells, am));
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Support_cell_stats", {H5S_UNLIMITED, 33});
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
            std::vector<std::string>({"area",
                                      "area__stddev",
                                      "area__max",
                                      "area__min",
                                      "area_preferential",
                                      "area_preferential__stddev",
                                      "area_preferential__max",
                                      "area_preferential__min",
                                      "perimeter",
                                      "perimeter__stddev",
                                      "perimeter__max",
                                      "perimeter__min",
                                      "shape_index",
                                      "shape_index__stddev",
                                      "shape_index__max",
                                      "shape_index__min",
                                      "num_neighbors",
                                      "num_neighbors__stddev",
                                      "num_neighbors__max",
                                      "num_neighbors__min",
                                      "num_hair_neighbors",
                                      "num_hair_neighbors__stddev",
                                      "num_hair_neighbors__max",
                                      "num_hair_neighbors__min",
                                      "hexatic_order",
                                      "hexatic_order__stddev",
                                      "hexatic_order__max",
                                      "hexatic_order__min",
                                      "rotation",
                                      "rotation__stddev",
                                      "rotation__max",
                                      "rotation__min",
                                      "count"
                                      }));


    }
); // end support_cell_stats_adaptor

/// Datamanager adaptor for cell area statistics
template <typename Cell, typename CellType = typename Cell::State::CellType>
auto bulk_cell_stats_adaptor = std::make_tuple(

    // name of the task
    "Bulk_cell_stats",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        
        AgentContainer<Cell> bulk_cells;
        bulk_cells.reserve(cells.size());
        std::copy_if(cells.begin(), cells.end(), std::back_inserter(bulk_cells),
                     [am](const auto& cell) {
                         return not am.is_boundary(cell);
                     });
        bulk_cells.shrink_to_fit();
        
        dataset->write(generate_statistics(bulk_cells, am));
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Bulk_cell_stats", {H5S_UNLIMITED, 33});
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
            std::vector<std::string>({"area",
                                      "area__stddev",
                                      "area__max",
                                      "area__min",
                                      "area_preferential",
                                      "area_preferential__stddev",
                                      "area_preferential__max",
                                      "area_preferential__min",
                                      "perimeter",
                                      "perimeter__stddev",
                                      "perimeter__max",
                                      "perimeter__min",
                                      "shape_index",
                                      "shape_index__stddev",
                                      "shape_index__max",
                                      "shape_index__min",
                                      "num_neighbors",
                                      "num_neighbors__stddev",
                                      "num_neighbors__max",
                                      "num_neighbors__min",
                                      "num_hair_neighbors",
                                      "num_hair_neighbors__stddev",
                                      "num_hair_neighbors__max",
                                      "num_hair_neighbors__min",
                                      "hexatic_order",
                                      "hexatic_order__stddev",
                                      "hexatic_order__max",
                                      "hexatic_order__min",
                                      "rotation",
                                      "rotation__stddev",
                                      "rotation__max",
                                      "rotation__min",
                                      "count"
                                      }));


    }
); // end bulk_cell_stats_adaptor

/// Datamanager adaptor for cell area statistics
template <typename Cell, typename CellType = typename Cell::State::CellType>
auto bulk_hair_cell_stats_adaptor = std::make_tuple(

    // name of the task
    "Bulk_hair_cell_stats",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        
        AgentContainer<Cell> bulk_cells;
        bulk_cells.reserve(cells.size());
        std::copy_if(cells.begin(), cells.end(), std::back_inserter(bulk_cells),
                     [am](const auto& cell) {
                         return not am.is_boundary(cell);
                     });
        bulk_cells.shrink_to_fit();
        
        AgentContainer<Cell> hair_bulk_cells;
        hair_bulk_cells.reserve(bulk_cells.size());
        std::copy_if(bulk_cells.begin(), bulk_cells.end(),
                     std::back_inserter(hair_bulk_cells),
                     [](const auto& cell) {
                         return cell->state.type == CellType::hair;
                     });
        hair_bulk_cells.shrink_to_fit();
        
        dataset->write(generate_statistics(hair_bulk_cells, am));
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Bulk_hair_cell_stats", {H5S_UNLIMITED, 33});
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
            std::vector<std::string>({"area",
                                      "area__stddev",
                                      "area__max",
                                      "area__min",
                                      "area_preferential",
                                      "area_preferential__stddev",
                                      "area_preferential__max",
                                      "area_preferential__min",
                                      "perimeter",
                                      "perimeter__stddev",
                                      "perimeter__max",
                                      "perimeter__min",
                                      "shape_index",
                                      "shape_index__stddev",
                                      "shape_index__max",
                                      "shape_index__min",
                                      "num_neighbors",
                                      "num_neighbors__stddev",
                                      "num_neighbors__max",
                                      "num_neighbors__min",
                                      "num_hair_neighbors",
                                      "num_hair_neighbors__stddev",
                                      "num_hair_neighbors__max",
                                      "num_hair_neighbors__min",
                                      "hexatic_order",
                                      "hexatic_order__stddev",
                                      "hexatic_order__max",
                                      "hexatic_order__min",
                                      "rotation",
                                      "rotation__stddev",
                                      "rotation__max",
                                      "rotation__min",
                                      "count"
                                      }));


    }
); // end bulk_hair_cell_stats_adaptor

/// Datamanager adaptor for cell area statistics
template <typename Cell, typename CellType = typename Cell::State::CellType>
auto bulk_support_cell_stats_adaptor = std::make_tuple(

    // name of the task
    "Bulk_support_cell_stats",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Statistics");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& cells = am.cells();
        
        AgentContainer<Cell> bulk_cells;
        bulk_cells.reserve(cells.size());
        std::copy_if(cells.begin(), cells.end(), std::back_inserter(bulk_cells),
                     [am](const auto& cell) {
                         return not am.is_boundary(cell);
                     });
        bulk_cells.shrink_to_fit();
        
        AgentContainer<Cell> support_bulk_cells;
        support_bulk_cells.reserve(bulk_cells.size());
        std::copy_if(bulk_cells.begin(), bulk_cells.end(),
                     std::back_inserter(support_bulk_cells),
                     [](const auto& cell) {
                         return cell->state.type == CellType::support;
                     });
        support_bulk_cells.shrink_to_fit();
        
        dataset->write(generate_statistics(support_bulk_cells, am));
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Bulk_support_cell_stats", {H5S_UNLIMITED, 33});
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
            std::vector<std::string>({"area",
                                      "area__stddev",
                                      "area__max",
                                      "area__min",
                                      "area_preferential",
                                      "area_preferential__stddev",
                                      "area_preferential__max",
                                      "area_preferential__min",
                                      "perimeter",
                                      "perimeter__stddev",
                                      "perimeter__max",
                                      "perimeter__min",
                                      "shape_index",
                                      "shape_index__stddev",
                                      "shape_index__max",
                                      "shape_index__min",
                                      "num_neighbors",
                                      "num_neighbors__stddev",
                                      "num_neighbors__max",
                                      "num_neighbors__min",
                                      "num_hair_neighbors",
                                      "num_hair_neighbors__stddev",
                                      "num_hair_neighbors__max",
                                      "num_hair_neighbors__min",
                                      "hexatic_order",
                                      "hexatic_order__stddev",
                                      "hexatic_order__max",
                                      "hexatic_order__min",
                                      "rotation",
                                      "rotation__stddev",
                                      "rotation__max",
                                      "rotation__min",
                                      "count"
                                      }));


    }
); // end bulk_support_cell_stats_adaptor

/// Datamanager adaptor for timepoints
/** The dataset to which the other energy adaptors link their coordinate time
 */
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