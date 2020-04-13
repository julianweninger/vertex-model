#ifndef UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH
#define UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"
#include "geometry.hh"

using namespace Utopia::DataIO;

namespace Utopia::Models::PCPVertex::DataIO{

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
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Total");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
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
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Time");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
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
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Linetension");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
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
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Areaelasticity");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
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
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Contractility");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
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
    [](auto& dataset, auto& model) {
        dataset->write(model.get_energy_cell_cell_polarity() / 
                       model.get_am().cells().size());
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_cell_polarity");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
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
    [](auto& dataset, auto& model) {
        dataset->write(model.get_energy_polarity_exclusion() / 
                       model.get_am().cells().size());
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Polarity_exclusion");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
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
    [](auto& dataset, auto& model) {
        dataset->write(model.get_energy_lagrange_net_polarisation() / 
                       model.get_am().cells().size());
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Lagrange_net_polarisation");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
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
    [](auto& dataset, auto& model) {
        dataset->write(model.get_energy_lagrange_const_concentration() / 
                       model.get_am().cells().size());
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Lagrange_const_concentration");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {}
    ,

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
    }
); // end lagrange_const_concentration_adaptor


/// Datamanager adaptor for vertex-position
auto vertex_position_adaptor = std::make_tuple(
    // name of the task
    "Vertex_position",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Vertex_position");
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
    [](auto& grp, auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "coordinate");
        hdfdataset->add_attribute("coords__coordinate", 
                                  std::vector<std::string>({"x", "y"}));
        hdfdataset->add_attribute("dim_name__1", "id");
        const auto domain = model.get_space()->get_domain_size();
        hdfdataset->add_attribute("Lx", domain[0]);
        hdfdataset->add_attribute("Ly", domain[1]);
        // For ids, the dimensions are trivial
        // hdfdataset->add_attribute("coords__coordinate", std::vector<std::size_t>{1, 1});
    }
    
); // end vertex position adaptor

/// Datamanager adaptor for cell position
template <typename SpaceVec>
auto cell_position_adaptor = std::make_tuple(
    // name of the task
    "Cell_position",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Cell_position");
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
        for (const auto& c : cells) {
            centers.push_back(am.barycenter_of(c));
        }
        dataset->write(centers.begin(), centers.end(),
                       [](auto&& pos) { return pos[0]; });
        dataset->write(centers.begin(), centers.end(),
                       [](auto&& pos) { return pos[1]; });

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
            {3, m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "coordinate");
        hdfdataset->add_attribute("coords__coordinate", 
                                  std::vector<std::string>({"cell_type",
                                                            "x", "y"}));
                                                            // "polarity_x",
                                                            // "polarity_y"}));
        hdfdataset->add_attribute("dim_name__1", "id");
        const auto domain = model.get_space()->get_domain_size();
        hdfdataset->add_attribute("Lx", domain[0]);
        hdfdataset->add_attribute("Ly", domain[1]);
    }    
); // end cell position adaptor

/// Datamanager adaptor for position
auto edge_link_adaptor = std::make_tuple(

    // name of the task
    "Edge_link",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Edge_link");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& vertices = model.get_am().vertices();
        const auto& edges = model.get_am().edges();
        dataset->write(edges.begin(), edges.end(),
                       [](const auto& edge) {
                           return static_cast<int>(edge->custom_links().a->id());
                       });
        dataset->write(edges.begin(), edges.end(),
                       [](const auto& edge) {
                           return static_cast<int>(edge->custom_links().b->id());
                       });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_am().edges().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "vertex");
        hdfdataset->add_attribute("coords__vertex", 
                                  std::vector<std::string>({"a", "b"}));
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
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Cell_area", {H5S_UNLIMITED, 3});
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        hdfdataset->add_attribute("dim_name__1", "properties");
        hdfdataset->add_attribute("coords__properties", 
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
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset("Time");
    },
    
    // attribute writer for basegroup
    [](auto& grp, auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
    }
); // end statistics_time_adaptor

} // namespace Utopia::Models::PCPVertex::DataIO


#endif // UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH