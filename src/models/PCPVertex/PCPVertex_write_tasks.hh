#ifndef UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH
#define UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::PCPVertex::DataIO{

/** Available datatree:
 *      - Energy
 *          - Energy_total
 *          - Energy_time
 *          - Energy_linetension
 *          - Energy_areaelasticity
 *          - Energy_contractility
 *          - Energy_cell_cell_polarity
 *          - Energy_polarity_exclusion
 *          - Energy_lagrange_net_polarisation
 *          - Energy_lagrange_const_concentration
 *      - Vertices (time series groups)
 *      - Cells (time series groups)
 *      - Edges (time series groups)
 *      - Statistics
 *          - Cell_area
 *          - Statistics_time
 */

/// Datamanager adaptor for total energy
auto energy_adaptor = std::make_tuple(

    // name of the task
    "Energy_total",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_energy() / model.get_am().cells().size());
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

/// Datamanager adaptor for total energy
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

/// Datamanager adaptor for total energy
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

/// Datamanager adaptor for total energy
auto contractility_adaptor = std::make_tuple(

    // name of the task
    "Energy_contractility",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        dataset->write(model.get_energy_edge_contractility() / 
                       model.get_am().cells().size() + 
                       model.get_energy_cell_contractility() / 
                       model.get_am().cells().size());
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Contractility");
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
); // end contractility_adaptor

/// Datamanager adaptor for total energy
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

/// Datamanager adaptor for total energy
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

/// Datamanager adaptor for total energy
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

/// Datamanager adaptor for total energy
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
        const auto domain = model.get_space()->get_domain_size();
        hdfdataset->add_attribute("Lx", domain[0]);
        hdfdataset->add_attribute("Ly", domain[1]);
        // For ids, the dimensions are trivial
        // hdfdataset->add_attribute("coords__coordinate", std::vector<std::size_t>{1, 1});
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
        for (const auto& c : cells) {
            centers.push_back(am.barycenter_of(c));
        }
        dataset->write(centers.begin(), centers.end(),
                       [](auto&& pos) { return pos[0]; });
        dataset->write(centers.begin(), centers.end(),
                       [](auto&& pos) { return pos[1]; });

        std::vector<double> areas;
        areas.reserve(cells.size());
        for (const auto& c : cells) { 
            areas.push_back(am.area_of(c)); 
        }
        dataset->write(areas);
        
        std::vector<double> perimeters;
        perimeters.reserve(cells.size());
        for (const auto& c : cells) {
            perimeters.push_back(am.perimeter_of(c));
        }
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
            {8, m.get_am().cells().size()});
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
                    "perimeter",
                    "shape_index",
                    "num_neighbors",
                    "num_hair_neighbors"
                }));
                    // "polarity_x",
                    // "polarity_y"}));
        hdfdataset->add_attribute("dim_name__1", "id");
        const auto domain = model.get_space()->get_domain_size();
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
        const auto& vertices = model.get_am().vertices();
        int running_id = 0;
        for (auto& v : vertices) {
            v->state.current_id = running_id++; // set identity within container
        }

        const auto& edges = model.get_am().edges();
        dataset->write(edges.begin(), edges.end(),
                       [](const auto& edge) {
                           return static_cast<int>(
                                    edge->custom_links().a->state.current_id);
                       });
        dataset->write(edges.begin(), edges.end(),
                       [](const auto& edge) {
                           return static_cast<int>(
                                    edge->custom_links().b->state.current_id);
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
        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute("coords__property", 
                                  std::vector<std::string>({"vertex_a",
                                                            "vertex_b"}));
        hdfdataset->add_attribute("dim_name__1", "id");

        
    }    
); // end edge link adaptor

/// Datamanager adaptor for total energy
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
); // end energy_adaptor

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

} // namespace Utopia::Models::PCPVertex::DataIO


#endif // UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH