#ifndef UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH
#define UTOPIA_MODELS_PCPVERTEX_WRITETASKS_HH

#include "utopia/data_io/hdfgroup.hh"
#include <utopia/data_io/data_manager/defaults.hh>
#include <complex>
#include <cmath>

using namespace Utopia::DataIO;

/** Available datatree:
 *      - Energy
 *          - Terms
 *          - Time (the linkded time for all energy adaptors)
 *          - Transitions
 *          - Interface_length
 *      - Vertices (time series groups)
 *      - Cells (time series groups)
 *      - Edges (time series groups)
 *      - Cell_energies (time series groups)
 *      - Edge_energies (time series groups)
 */
namespace Utopia::Models::PCPVertex::DataIO{

/// Datamanager adaptor for total energy
auto energy_adaptor = std::make_tuple(

    // name of the task
    "Energy",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        double E_total = 0.;
        std::map<std::string, double> values({});
        const auto& terms = model.get_work_function_terms();
        for (const auto& name : model.get_work_function_term_register()) {
            values[name] = 0.;
        }
        for (const auto& [name, term] : terms) {
            values[name] = term->compute_energy();
            E_total += values[name];
        }
        values["total"] = E_total;

        std::vector<double> write({});
        for (const auto& [n, v] : values) {
            write.push_back(v);
        }

        dataset->write(write);
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& model) -> decltype(auto) {
        std::size_t N = model.fix_number_work_function_terms();
        return group->open_dataset("Terms", {H5S_UNLIMITED, N+1});
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");
        
        std::set<std::string> names({});
        for (const auto& name : model.get_work_function_term_register()) {
            names.insert(name);
        }
        names.insert("total");

        hdfdataset->add_attribute("dim_name__1", "term");
        hdfdataset->add_attribute("coords_mode__term", "values");
        hdfdataset->add_attribute("coords__term", 
                                  std::vector(names.begin(), names.end()));
        
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
template <typename SpaceVec>
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
                            return static_cast<float>(cell->state.type);
                       });
        
        std::vector<SpaceVec> centers;
        centers.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(centers),
                       [am](const auto& c) { 
                           return am.barycenter_of(c);
                        });
        dataset->write(centers.begin(), centers.end(),
                       [](auto&& pos) { 
                           return static_cast<float>(pos[0]); 
                        });
        dataset->write(centers.begin(), centers.end(),
                       [](auto&& pos) { 
                           return static_cast<float>(pos[1]); 
                        });

        std::vector<float> areas;
        areas.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(areas),
                       [am](const auto& c) { 
                           return static_cast<float>(am.area_of(c));
                        });
        dataset->write(areas);

        dataset->write(cells.begin(), cells.end(),
            [](const auto& c) { 
                return static_cast<float>(c->state.area_preferential);
            });
        
        std::vector<float> perimeters;
        perimeters.reserve(cells.size());
        std::transform(cells.begin(), cells.end(),
                       std::back_inserter(perimeters),
                       [am](const auto& c) {
                           return static_cast<float>(am.perimeter_of(c));
                        });
        dataset->write(perimeters);

        dataset->write(cells.begin(), cells.end(),
                       [am](const auto& cell) {
                            return static_cast<float>(
                                   am.neighbors_of(cell).size());
                       });
        dataset->write(cells.begin(), cells.end(),
                       [am](const auto& cell) {
                            auto nbs = am.neighbors_of(cell);
                            std::size_t N = std::count_if(
                                nbs.begin(), nbs.end(),
                                [cell](const auto& nb) {
                                    return nb->state.type == cell->state.type;
                                }
                            );
                            return static_cast<float>(N);
                       });


        std::vector<float> hex_order;
        hex_order.reserve(cells.size());
        std::vector<float> hex_order_corr;
        hex_order_corr.reserve(cells.size());
        std::vector<float> hex_order_N;
        hex_order_N.reserve(cells.size());
        for (const auto& cell : cells) {
            const auto [hex, hex_corr, hex_N] = am.hexatic_order_of(
                cell, cell->state.type, 2, cell->state.type
            );
            hex_order.push_back(hex);
            hex_order_corr.push_back(hex_corr);
            hex_order_N.push_back(hex_N);
        }
        dataset->write(hex_order);
        dataset->write(hex_order_corr);
        dataset->write(hex_order_N);

        // // write cell polarity wrt curved x axis
        // dataset->write(cells.begin(), cells.end(),
        //                [model](const auto& cell) {
        //                     double curvature = model.get_curvature_boundary();
        //                     if (fabs(curvature) < 1.e-8) {
        //                         double pol = cell->state.polarity();
        //                         return static_cast<float>(pol);
        //                     }

        //                     if (model.get_space()->periodic) {
        //                         throw std::runtime_error("Not implemented "
        //                             "cell polarity in curved periodic BC.");
        //                     }

        //                     SpaceVec origin({0., -1./curvature});

        //                     const auto& am = model.get_am();
        //                     SpaceVec pos = am.barycenter_of(cell) - origin;

        //                     double theta = std::atan2(pos[0], pos[1]);

        //                     double pol = cell->state.polarity(theta);
        //                     return static_cast<float>(pol);
        //                });

        std::vector<SpaceVec> qs;
        std::transform(cells.begin(), cells.end(), std::back_inserter(qs),
                       [am](const auto& c) {
                           arma::mat22 q = am.elongation_of(c);

                            double abs_q = sqrt(arma::trace(q * q.t())/2.);
                            double ori = atan2(q[2], q[0]) / 2.;

                            return abs_q * SpaceVec({cos(ori), sin(ori)});
                       });
        dataset->write(qs.begin(), qs.end(), 
                       [](auto&& q) {
                           return static_cast<float>(q[0]);
                        });
        dataset->write(qs.begin(), qs.end(), 
                       [](auto&& q) {
                           return static_cast<float>(q[1]);
                        });

        dataset->write(
            cells.begin(), cells.end(),
            [am](const auto& cell) {
                return static_cast<float>(am.is_boundary(cell));
            });

        
        for (const auto& [name, term] : model.get_work_function_terms()) {
            for (const auto& data : term->write_cell_properties()) {
                dataset->write(std::vector<float>(data.begin(), data.end()));
            }
        }
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        std::size_t cnt = 14;
        for (const auto& [name, term] : m.get_work_function_terms()) {
            cnt += term->write_task_cell_properties_names().size();
        }

        return group->open_dataset(std::to_string(m.get_time()), 
            {cnt, m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        std::vector<std::string> properties({
            "cell_type",
            "x",
            "y",
            "area",
            "area_preferential",
            "perimeter",
            "num_neighbors",
            "num_neighbors__self",
            "hexatic_order",
            "hexatic_order__corr",
            "hexatic_order__N",
            "q_x",
            "q_y",
            "is_boundary"
        });
        for (const auto& [name, term] : model.get_work_function_terms()) {
            const auto task_properties = term->write_task_cell_properties_names();
            for (const auto& task_property : task_properties) {
                properties.insert(properties.end(), name + "__" + task_property);
            }
        }

        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute(
            "coords__property", 
            std::vector<std::string>(properties.begin(), properties.end())
        );

        const auto& cells = model.get_am().cells();
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
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
template <typename SpaceVec>
auto edges_adaptor = std::make_tuple(

    // name of the task
    "Edges",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Edges");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& am = model.get_am();
        const auto& edges = am.edges();
        dataset->write(edges.begin(), edges.end(),
                    [](const auto& edge) {
                        return static_cast<float>(edge->custom_links().a->id());
                    });
        dataset->write(edges.begin(), edges.end(),
                    [](const auto& edge) {
                        return static_cast<float>(edge->custom_links().b->id());
                    });
        
        dataset->write(edges.begin(), edges.end(),
                       [am](const auto& edge) {
                           return static_cast<float>(am.length_of(edge));
                       });        
        dataset->write(edges.begin(), edges.end(),
                       [am](const auto& edge) {
                           SpaceVec displ = am.displacement(edge);
                           displ /= arma::norm(displ);
                           SpaceVec axis({1., 0.});

                            // nematic angle in [0, pi / 2]
                           double angle = std::acos(arma::dot(displ, axis));

                           return static_cast<float>(angle);
                       });
                       
        dataset->write(edges.begin(), edges.end(),
                       [am](const auto& edge) {
                           auto [a, b] = am.template adjoints_of<true>(edge);

                           if (a == nullptr) { return static_cast<float>(-1); }
                           return static_cast<float>(a->id());
                       });
        dataset->write(edges.begin(), edges.end(),
                       [am](const auto& edge) {
                           auto [a, b] = am.template adjoints_of<true>(edge);

                           if (b == nullptr) { return static_cast<float>(-1); }
                           return static_cast<float>(b->id());
                       });
                       
        dataset->write(edges.begin(), edges.end(),
                       [am](const auto& edge) {
                           auto [a, b] = am.template adjoints_of<true>(edge);

                           if (a == nullptr) { return static_cast<float>(-1); }
                           return static_cast<float>(a->state.type);
                       });
        dataset->write(edges.begin(), edges.end(),
                       [am](const auto& edge) {
                           auto [a, b] = am.template adjoints_of<true>(edge);

                           if (b == nullptr) { return static_cast<float>(-1); }
                           return static_cast<float>(b->state.type);
                       });

        for (const auto& [name, term] : model.get_work_function_terms()) {
            for (const auto& data : term->write_edge_properties()) {
                dataset->write(std::vector<float>(data.begin(), data.end()));
            }
        }
    },
                
    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        std::size_t cnt = 8;
        for (const auto& [name, term] : m.get_work_function_terms()) {
            cnt += term->write_task_edge_properties_names().size();
        }

        return group->open_dataset(std::to_string(m.get_time()), 
            {cnt, m.get_am().edges().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        std::vector<std::string> properties({
            "vertex_a",
            "vertex_b",
            "length",
            "orientation",
            "id_alpha",
            "id_beta",
            "type_alpha",
            "type_beta"
        });
        for (const auto& [name, term] : model.get_work_function_terms()) {
            const auto task_properties = term->write_task_edge_properties_names();
            for (const auto& task_propery : task_properties) {
                properties.insert(properties.end(), name + "__" + task_propery);
            }
        }

        const auto& edges = model.get_am().edges();
        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute(
            "coords__property", 
            std::vector<std::string>(properties.begin(), properties.end())
        );
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        std::vector<std::size_t> ids{};
        ids.reserve(edges.size());
        std::transform(edges.begin(), edges.end(), std::back_inserter(ids),
                       [](const auto& e) { return e->id(); });
        hdfdataset->add_attribute("coords__id", ids);


        const SpaceVec skew = model.get_space()->get_skew();
        const double curvature = model.get_space()->get_curvature();

        hdfdataset->add_attribute("skew_x", skew[0]);
        hdfdataset->add_attribute("skew_y", skew[1]);
        hdfdataset->add_attribute("curvature", curvature);
    }     
); // end edge link adaptor

/// Datamanager adaptor for edge energies and tensions
auto edge_energies_adaptor = std::make_tuple(

    // name of the task
    "Edge_energies",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Edge_energies");
    },

    // writer function
    [](auto& dataset, auto& model) {
        for (const auto& [name, term] : model.get_work_function_terms()) {
            for (const auto& data : term->write_edge_energies()) {
                dataset->write(std::vector<float>(data.begin(), data.end()));
            }
        }
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        std::size_t cnt = 0;
        for (const auto& [name, term] : m.get_work_function_terms()) {
            cnt += term->write_task_edge_energies_names().size();
        }
        return group->open_dataset(std::to_string(m.get_time()), 
            {cnt, m.get_am().edges().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        std::vector<std::string> terms({});
        for (const auto& [name, term] : model.get_work_function_terms()) {
            for (const auto& n : term->write_task_edge_energies_names()) {
                terms.push_back(name + "__" + n);
            }
        }
        hdfdataset->add_attribute("dim_name__0", "term");
        hdfdataset->add_attribute("coords__term", terms);

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

/// Datamanager adaptor for cell energies and pressures
auto cell_energies_adaptor = std::make_tuple(
    // name of the task
    "Cell_energies",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Cell_energies");
    },

    // writer function
    [](auto& dataset, auto& model) {
        for (const auto& [name, term] : model.get_work_function_terms()) {
            for (const auto& data : term->write_cell_energies()) {
                dataset->write(std::vector<float>(data.begin(), data.end()));
            }
        }
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        std::size_t cnt = 0;
        for (const auto& [name, term] : m.get_work_function_terms()) {
            cnt += term->write_task_cell_energies_names().size();
        }
        return group->open_dataset(std::to_string(m.get_time()), 
            {cnt, m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, auto& model) {
        std::vector<std::string> terms({});
        for (const auto& [name, term] : model.get_work_function_terms()) {
            for (const auto& n : term->write_task_cell_energies_names()) {
                terms.push_back(name + "__" + n);
            }
        }
        hdfdataset->add_attribute("dim_name__0", "term");
        hdfdataset->add_attribute("coords__term", terms);

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


/// Datamanager adaptor for the cluster of cells
/** Cluster id labels of cells with type 1
 *  Attributes are:
 *      -# Lx: Domain size in x coordinate
 *      -# Ly: Domain size in y coordinate
 */
auto cell_cluster_adaptor = std::make_tuple(
    // name of the task
    "Cell_cluster",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Cell_cluster");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& cells = model.get_am().cells();
        auto clusters = model.get_cluster_ids(1);
        dataset->write(cells.begin(), cells.end(),
            [clusters](const auto& cell) {
                return clusters.at(cell);
            });
    },

    // builder function
    [](auto& group, auto& m) -> decltype(auto) {
        return group->open_dataset(std::to_string(m.get_time()), 
            {1, m.get_am().cells().size()});
    },

    // attribute writer for basegroup
    [](auto& grp, [[maybe_unused]] auto& m) {
        grp->add_attribute("content", "time_series");},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "property");
        hdfdataset->add_attribute("coords__property", 
                                  std::vector<std::string>({"cluster_id"}));
        hdfdataset->add_attribute("dim_name__1", "id");
        hdfdataset->add_attribute("coords_mode__id", "values");
        const auto& cells = model.get_am().cells();
        std::vector<std::size_t> ids{};
        ids.reserve(cells.size());
        std::transform(cells.begin(), cells.end(), std::back_inserter(ids),
                       [](const auto& c) { return c->id(); });
        hdfdataset->add_attribute("coords__id", ids);

    }    
); // end hair cluster adaptor

auto transition_adaptor = std::make_tuple(

    // name of the task
    "Transitions",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        std::vector<double> stats{};
        stats.reserve(9);
        
        stats.push_back(model.get_num_T1s());
        stats.push_back(model.get_T1_frequency());
        stats.push_back(model.get_T1_frequency_accumulated());
        stats.push_back(model.get_num_T1s_attempted());
        stats.push_back(model.get_T1_attempt_frequency());
        stats.push_back(model.get_T1_attempt_frequency_accumulated());
        stats.push_back(model.get_num_T2s());
        stats.push_back(model.get_T2_frequency());
        stats.push_back(model.get_T2_frequency_accumulated());
        
        dataset->write(stats);
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Transitions", { H5S_UNLIMITED, 9 });
    },
    
    // attribute writer for basegroup
    []([[maybe_unused]] auto& grp, [[maybe_unused]] auto& m) {},

    // attribute writer for dataset
    [](auto& hdfdataset, [[maybe_unused]] auto& model) {
        hdfdataset->add_attribute("dim_name__0", "time");
        hdfdataset->add_attribute("coords_mode__time", "linked");
        hdfdataset->add_attribute("coords__time", "Time");

        hdfdataset->add_attribute("dim_name__1", "type");
        hdfdataset->add_attribute("coords__type", 
                std::vector<std::string>({
                    "num_T1s",
                    "T1_frequency",
                    "T1_frequency_accumulated",
                    "num_T1s_attempted",
                    "T1_attempt_frequency",
                    "T1_attempt_frequency_accumulated",
                    "num_T2s",
                    "T2_frequency",
                    "T2_frequency_accumulated"
                }));
    }
); // Transitions


/// Datamanager adaptor for length of heterotypic cell interfaces
auto interface_length_adaptor = std::make_tuple(

    // name of the task
    "Interface_length",

    // basegroup builder
    [](std::shared_ptr<HDFGroup>&& grp) -> std::shared_ptr<HDFGroup> {
        return grp->open_group("Energy");
    },

    // writer function
    [](auto& dataset, auto& model) {
        const auto& cm = model.get_am();
        const auto& edges = cm.edges();

        double L = 0.;
        for (const auto& edge : edges) {
            const auto [a, b] = cm.adjoints_of(edge);
            if (not a or not b) { continue; }
            if (a->state.type != b->state.type) {
                L += cm.length_of(edge);
            }
        }

        dataset->write(L);
    },

    // builder function
    [](auto& group, [[maybe_unused]] auto& m) -> decltype(auto) {
        return group->open_dataset("Interface_length");
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
); // end interface_length_adaptor



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