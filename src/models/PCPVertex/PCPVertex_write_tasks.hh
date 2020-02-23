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
        dataset->write(model.get_energy_normalised());
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
        dataset->write(model.get_energy_linetension_normalised());
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
        dataset->write(model.get_energy_areaelasticity_normalised());
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
        dataset->write(model.get_energy_contractility_normalised());
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
        dataset->write(model.get_energy_cell_cell_polarity_normalised());
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
        dataset->write(model.get_energy_polarity_exclusion_normalised());
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
        dataset->write(model.get_energy_lagrange_net_polarisation_normalised());
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
        dataset->write(model.get_energy_lagrange_const_concentration_normalised());
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
        const auto& vertices = model.get_vertices();
        dataset->write(vertices.begin(), vertices.end(),
                        [](auto&& vertex) {
                            return static_cast<double>(vertex.lock()->x);
                        });
        dataset->write(vertices.begin(), vertices.end(),
                        [](auto&& vertex) {
                            return static_cast<double>(vertex.lock()->y);
                        });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_vertices().size()});
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
        auto [Lx, Ly] = model.get_domain_size();
        hdfdataset->add_attribute("Lx", Lx);
        hdfdataset->add_attribute("Ly", Ly);
        // For ids, the dimensions are trivial
        // hdfdataset->add_attribute("coords__coordinate", std::vector<std::size_t>{1, 1});
    }
    
); // end vertex position adaptor

/// Datamanager adaptor for cell position
template <bool periodic_bc>
auto cell_position_adaptor = std::make_tuple(

    // name of the task
    "Cell_position",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Cell_position");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& cells = model.get_cells();
        dataset->write(cells.begin(), cells.end(),
                        [](auto&& cell) {
                            return static_cast<double>(cell.lock()->type);
                        });

        dataset->write(cells.begin(), cells.end(),
                        [](auto&& cell) {
                            return static_cast<double>(cell.lock()->s->x);
                        });
        dataset->write(cells.begin(), cells.end(),
                        [](auto&& cell) {
                            return static_cast<double>(cell.lock()->s->y);
                        });
        // the polarity
        auto polarities_y = std::make_shared<std::vector<double>>();
        dataset->write(cells.begin(), cells.end(), [polarities_y](auto & cell)
        {
            double polarity_x = 0.;
            double polarity_y = 0.;
            for (auto [e, flip] : cell.lock()->edges_ordered) {
                auto a = *e->a, b = *e->b;
                if (flip) {
                    std::swap(a, b);
                }
                auto disp = displacement<periodic_bc>(a, b);

                double sigma = e->get_sigma(cell.lock());
                
                polarity_x += sigma * std::get<0>(disp);
                polarity_y += sigma * std::get<1>(disp);
            }

            if (cell.lock()->area_sgn == 1) {
                polarities_y->push_back(polarity_y);
                return polarity_x;
            }
            else {
                polarities_y->push_back(-1. * polarity_y);
                return -1. * polarity_x;
            }
        });
        dataset->write(polarities_y->begin(), polarities_y->end(), [](auto val) {
            return val;
        });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {5, m.get_cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        hdfdataset->add_attribute("dim_name__0", "coordinate");
        hdfdataset->add_attribute("coords__coordinate", 
                                  std::vector<std::string>({"cell_type",
                                                            "x", "y",
                                                            "polarity_x",
                                                            "polarity_y"}));
        hdfdataset->add_attribute("dim_name__1", "id");
        auto [Lx, Ly] = model.get_domain_size();
        hdfdataset->add_attribute("Lx", Lx);
        hdfdataset->add_attribute("Ly", Ly);
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
        auto vertices = model.get_vertices();
        int running_id = 0;
        for (auto &v : vertices) {
            v.lock()->current_id = running_id++;
        }

        const auto& edges = model.get_edges();
        dataset->write(edges.begin(), edges.end(),
                        [](auto&& edge) {
                            return static_cast<int>(edge.lock()->a->current_id);
                        });
        dataset->write(edges.begin(), edges.end(),
                        [](auto&& edge) {
                            return static_cast<int>(edge.lock()->b->current_id);
                        });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {2, m.get_edges().size()});
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

} // namespace Utopia::Models::PCPVertex::DataIO


#endif // UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH