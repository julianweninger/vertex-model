#define BOOST_TEST_MODULE PCPTopologyTestOperations

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "utils.hh"

#include "../../PCPVertex/PCPVertex.hh"
#include "../../PCPVertex/energy.hh"
#include "../../PCPVertex/algorithm.hh"
#include "../../PCPVertex/operations.hh"
#include "../../PCPVertex/PCPVertex_write_tasks.hh"

#include "../PCPTopology.hh"
#include "../operations.hh"
#include "../PCPTopology_write_tasks.hh"

#include "../../Collier/Collier.hh"
#include "../../NotchDelta/NotchDelta.hh"
#include "../../NotchDelta/Differentiation_write_tasks.hh"
#include "../../Collier/Collier_write_tasks.hh"

using namespace Utopia;
using namespace Utopia::Models::PCPVertex::OperationCollection;
using Utopia::DataIO::Config;

using CellType = Models::PCPVertex::PCPVertex::CellType;
using SpaceVec = Models::PCPVertex::PCPVertex::SpaceVec;

struct Fixture : ModelFixture {
    PseudoParent<> pp;
    
    Models::PCPVertex::PCPVertex vertex_model;
    
    Config cfg;

    Utopia::Models::PCPVertex::MinimizationParams default_minim_params;

    Fixture ()
    :
        ModelFixture(),
        pp("vertex_model_cfg.yml"),
        vertex_model("PCPVertex", pp, {},
                     std::make_tuple(
                        Models::PCPVertex::DataIO::time_energy_adaptor)),
        cfg(YAML::LoadFile("test_operations.yml")),
        default_minim_params(cfg["default_minimization"])
    { }
};

std::pair<double, double> get_statistics (std::vector<double>& values)
{
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    double mean = sum / values.size();

    double sq_sum = std::inner_product(values.begin(), values.end(),
                                       values.begin(), 0.0);
    double stddev = std::sqrt(sq_sum / values.size() - mean * mean);

    return std::make_pair(mean, stddev);
}

BOOST_AUTO_TEST_CASE(test_PCPTopology_apply_operation) {
    Config cfg(YAML::LoadFile("test_operations.yml"));
    Utopia::Models::PCPVertex::MinimizationParams default_minim_params(
        get_as<Config>("default_minimization", cfg));

    PseudoParent<> pseudo_parent("topology_model_cfg.yml");
    Models::PCPVertex::PCPTopology model("PCPTopology", pseudo_parent, {},
        std::make_tuple(Models::PCPVertex::DataIO::time_energy_adaptor));
    std::string name;
    std::size_t prev_time = model.get_continuous_time();
    std::size_t time;
        
    name = "apply_operation_prolog";
    model.register_operation(build_jiggle(
        name, get_as<Config>(name, cfg), default_minim_params));
    model.prolog();
    time = model.get_continuous_time();
    BOOST_TEST(time - prev_time == 8); // 4 iterations, minimize every
    prev_time = time;

    std::function<std::size_t(Models::PCPVertex::PCPTopology&,
                              std::string)> apply = 
        [cfg, default_minim_params]
        (Models::PCPVertex::PCPTopology& model, std::string name)
    {
        std::size_t time = model.get_continuous_time();
        model.register_operation(build_jiggle(
            name, get_as<Config>(name, cfg), default_minim_params));
        model.iterate();

        return (model.get_continuous_time() - time);
    };

    BOOST_TEST(apply(model, "apply_operation") == 10);
    // NOTE 5 iterations, minimize every

    BOOST_TEST(apply(model, "apply_operation_every") == 10);
    // NOTE 5 iterations, minimize every

    BOOST_TEST(apply(model, "apply_operation_once") == 2);
    // NOTE 5 iterations, minimize once

    BOOST_TEST(apply(model, "apply_operation_manual") == 0);
    // NOTE 5 iterations, minimize every

    BOOST_TEST(apply(model, "apply_operation_repeat") == 15);
    // NOTE 5 iterations, minimize every 3 times
        
    name = "apply_operation_epilog";
    model.register_operation(build_jiggle(
        name, get_as<Config>(name, cfg), default_minim_params));
    prev_time = model.get_continuous_time();
    model.epilog();
    time = model.get_continuous_time();
    BOOST_TEST(time - prev_time == 8); // 4 iterations, minimize every
    prev_time = time;
}

BOOST_FIXTURE_TEST_SUITE (test_PCPTopology_operations, Fixture)

    BOOST_AUTO_TEST_CASE(test_PCPTopology_jiggle) {
        const std::string name = "jiggle";
        auto [operation, params] = build_jiggle(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_differentiate_random) {
        const std::string name = "differentiate_random";
        auto [operation, params] = build_differentiate_random(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        auto cnt = std::count_if(cells.begin(), cells.end(),
                                 [](const auto& cell) {
                                        return cell->state.type == 
                                               CellType::hair; });

        BOOST_CHECK_CLOSE(double(cnt)/cells.size(), 0.25, 25);

        cnt = std::count_if(cells.begin(), cells.end(),
                            [](const auto& cell) {
                                return cell->state.type == 
                                        CellType::progenitor; });
        BOOST_TEST(cnt == 0);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_differentiate_random_fraction) {
        const std::string name = "differentiate_random_fraction";
        auto [operation, params] = build_differentiate_random(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        auto cnt = std::count_if(cells.begin(), cells.end(),
                                 [](const auto& cell) {
                                        return cell->state.type == 
                                               CellType::hair; });

        BOOST_TEST(cnt == 14); // Rounded from 0.23 * 64 = 14.72 cells 

        cnt = std::count_if(cells.begin(), cells.end(),
                            [](const auto& cell) {
                                return cell->state.type == 
                                        CellType::progenitor; });
        BOOST_TEST(cnt == 0);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_differentiate_Collier) {
        using Utopia::Models::Collier::Collier;

        std::string name = "differentiate_Collier";

        std::shared_ptr<Collier> collier(new Collier(
            "Collier", pp, 
            get_as<Config>("Collier", 
                           get_as<Config>(name, cfg)),
            std::make_tuple(
                Utopia::Models::Differentiation::DataIO::density_time)));
        auto collier_prolog = std::make_shared<bool>(false);

        auto [operation, params] = build_differentiate_Collier(
            name, get_as<Config>(name, cfg), default_minim_params,
            collier, collier_prolog);

        // This should call the prolog and iterate a first 100 steps
        operation(vertex_model);

        BOOST_TEST(collier_prolog);
        BOOST_TEST(collier->get_time() == 100);

        const auto& cells = vertex_model.get_am().cells();
        auto cnt = std::count_if(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->custom_links().c_cell; });
        BOOST_TEST(cnt == cells.size());

        // Continue iteration
        operation(vertex_model);
        BOOST_TEST(collier->get_time() == 200);

        cnt = std::count_if(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->custom_links().c_cell; });
        BOOST_TEST(cnt == cells.size());
        
        
        // Continue iteration with new task
        name = "differentiate_Collier_second_task";
        auto [new_operation, new_params] = build_differentiate_Collier(
            name, get_as<Config>(name, cfg), default_minim_params,
            collier, collier_prolog);
        
        new_operation(vertex_model);
        
        BOOST_TEST(collier->get_time() == 210);

        cnt = std::count_if(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->custom_links().c_cell; });
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_differentiate_NotchDelta) {
        using Utopia::Models::PCPVertex::PCPVertex;
        using Utopia::Models::NotchDelta::NotchDelta;

        // Function to compare the correct mapping of cells
        std::function<void(const PCPVertex&)> check_mapping = []
                (const PCPVertex& vertex_model)
        {
            const auto& am = vertex_model.get_am();
            const auto& cells = am.cells();

            // test that all cells have been mapped
            auto cnt = std::count_if(cells.begin(), cells.end(),
                [](const auto& cell) {
                    return cell->custom_links().nd_cell; });
            BOOST_TEST(cnt == cells.size());

            // test that links are correct
            for (const auto& cell : cells) {
                const auto& neighbors = am.neighbors_of(cell);
                const auto& nd_cell = cell->custom_links().nd_cell;
                const auto& nd_neighbors =  nd_cell->custom_links().neighbors;
                
                // have same number of neighbors
                BOOST_TEST(neighbors.size() == nd_neighbors.size());

                // the cell's neighbors appear in nd_cell's neighborhood
                for (auto& n : neighbors) {
                    const auto& nd_n = n->custom_links().nd_cell;
                    bool found = (   std::find(nd_neighbors.begin(),
                                            nd_neighbors.end(), nd_n)
                                  != nd_neighbors.end());
                    BOOST_TEST(found);
                }

                // the state's have been synchronized
                auto type = cell->state.type;
                if (type == PCPVertex::CellType::progenitor) {
                    BOOST_TEST(   nd_cell->state.cell_type
                               == NotchDelta::CellType::progenitor);
                }
                else if (type == PCPVertex::CellType::hair) {
                    BOOST_TEST(   nd_cell->state.cell_type
                               == NotchDelta::CellType::hair);
                }
                else if (type == PCPVertex::CellType::support) {
                    BOOST_TEST(   nd_cell->state.cell_type
                               == NotchDelta::CellType::support);
                }
            }
        };

        std::string name = "differentiate_NotchDelta";

        std::shared_ptr<NotchDelta> notch_delta(new NotchDelta(
            "NotchDelta", pp, 
            get_as<Config>("NotchDelta", 
                           get_as<Config>(name, cfg)),
            std::make_tuple(
                Utopia::Models::Differentiation::DataIO::density_time)));
        auto notch_delta_prolog = std::make_shared<bool>(false);

        auto [operation, params] = build_differentiate_NotchDelta(
            name, get_as<Config>(name, cfg), default_minim_params,
            notch_delta, notch_delta_prolog);

        // This should call the prolog and iterate a first 100 steps
        operation(vertex_model);

        BOOST_TEST(notch_delta_prolog);
        BOOST_TEST(notch_delta->get_time() == 100);

        check_mapping(vertex_model);


        // Continue iteration
        operation(vertex_model);
        BOOST_TEST(notch_delta->get_time() == 200);
        
        check_mapping(vertex_model);
        
        
        // Continue iteration with new task
        name = "differentiate_NotchDelta_second_task";
        auto [new_operation, new_params] = build_differentiate_NotchDelta(
            name, get_as<Config>(name, cfg), default_minim_params,
            notch_delta, notch_delta_prolog);
        
        new_operation(vertex_model);
        BOOST_TEST(notch_delta->get_time() == 210);
        
        check_mapping(vertex_model);
    }

    // Test the NotchDelta differentiation when NotchInitial condition is 
    // overwritten
    BOOST_AUTO_TEST_CASE(test_PCPTopology_differentiate_NotchDelta_skip_init) {
        using Utopia::Models::PCPVertex::PCPVertex;
        using Utopia::Models::NotchDelta::NotchDelta;

        // Function to compare the correct mapping of cells
        std::function<void(const PCPVertex&)> check_mapping = []
                (const PCPVertex& vertex_model)
        {
            const auto& am = vertex_model.get_am();
            const auto& cells = am.cells();

            // test that all cells have been mapped
            auto cnt = std::count_if(cells.begin(), cells.end(),
                [](const auto& cell) {
                    return cell->custom_links().nd_cell; });
            BOOST_TEST(cnt == cells.size());

            // test that links are correct
            for (const auto& cell : cells) {
                const auto& neighbors = am.neighbors_of(cell);
                const auto& nd_cell = cell->custom_links().nd_cell;
                const auto& nd_neighbors =  nd_cell->custom_links().neighbors;
                
                // have same number of neighbors
                BOOST_TEST(neighbors.size() == nd_neighbors.size());

                // the cell's neighbors appear in nd_cell's neighborhood
                for (auto& n : neighbors) {
                    const auto& nd_n = n->custom_links().nd_cell;
                    bool found = (   std::find(nd_neighbors.begin(),
                                            nd_neighbors.end(), nd_n)
                                  != nd_neighbors.end());
                    BOOST_TEST(found);
                }

                // the state's have been synchronized
                auto type = cell->state.type;
                if (type == PCPVertex::CellType::progenitor) {
                    BOOST_TEST(   nd_cell->state.cell_type
                               == NotchDelta::CellType::progenitor);
                }
                else if (type == PCPVertex::CellType::hair) {
                    BOOST_TEST(   nd_cell->state.cell_type
                               == NotchDelta::CellType::hair);
                }
                else if (type == PCPVertex::CellType::support) {
                    BOOST_TEST(   nd_cell->state.cell_type
                               == NotchDelta::CellType::support);
                }
            }
        };

        std::string name = "differentiate_NotchDelta_skip_init";
        auto process_cfg(get_as<Config>(name, cfg));
        auto Random_prolif(get_as<Config>("Random", process_cfg));
        auto NotchDelta_prolif(get_as<Config>("NotchDelta", process_cfg));
        auto NotchDelta_2nd_prolif(get_as<Config>("NotchDelta_2nd_task",
                                                     process_cfg));



        // the random differentiation
        auto [rnd_operation, rnd_params] = build_differentiate_random(
            name, Random_prolif, default_minim_params);


        // run the random differentiation
        // Terminally differentiated cells
        rnd_operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        auto cnt_hair = std::count_if(cells.begin(), cells.end(),
                                      [](const auto& cell) {
                                            return (   cell->state.type
                                                    == CellType::hair); });

        BOOST_CHECK_CLOSE(double(cnt_hair)/cells.size(), 0.35, 35);

        auto cnt = std::count_if(cells.begin(), cells.end(),
                            [](const auto& cell) {
                                return cell->state.type == 
                                        CellType::progenitor; });
        BOOST_TEST(cnt == 0);

        // save states for next try
        std::vector<std::size_t> cell_states;
        cell_states.reserve(cells.size());
        for (const auto& cell : cells) {
            cell_states.push_back(cell->state.type);
        }


        // The NotchDelta model
        // the notch differentiation
        std::shared_ptr<NotchDelta> notch_delta(new NotchDelta(
            "NotchDelta", pp, 
            get_as<Config>("NotchDelta", NotchDelta_prolif),
            std::make_tuple(
                Utopia::Models::Differentiation::DataIO::density_time)));
        auto notch_delta_prolog = std::make_shared<bool>(false);
        auto [operation, params] = build_differentiate_NotchDelta(
            name, NotchDelta_prolif, default_minim_params,
            notch_delta, notch_delta_prolog);
        

        operation(vertex_model);

        BOOST_TEST(notch_delta_prolog);
        BOOST_TEST(notch_delta->get_time() == 0);

        std::vector<std::size_t> new_cell_states;
        new_cell_states.reserve(cells.size());
        for (const auto& cell : cells) {
            new_cell_states.push_back(cell->state.type);
        }
        
        BOOST_CHECK_EQUAL_COLLECTIONS(
            cell_states.begin(), cell_states.end(), 
            new_cell_states.begin(), new_cell_states.end());

        check_mapping(vertex_model);

        auto [scnd_operation, scnd_params] = build_differentiate_NotchDelta(
            name, NotchDelta_2nd_prolif, default_minim_params,
            notch_delta, notch_delta_prolog);

        scnd_operation(vertex_model);

        BOOST_TEST(notch_delta_prolog);
        BOOST_TEST(notch_delta->get_time() == 100);

        auto new_cnt_hair = std::count_if(cells.begin(), cells.end(),
                                      [](const auto& cell) {
                                            return (   cell->state.type
                                                    == CellType::hair); });

        BOOST_TEST(new_cnt_hair == cnt_hair);

        cnt = std::count_if(cells.begin(), cells.end(),
                            [](const auto& cell) {
                                return cell->state.type == 
                                        CellType::progenitor; });
        BOOST_TEST(cnt == 0);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_area) {
        const std::string name = "increment_area";
        auto [operation, params] = build_increment_area(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();

        std::vector<double> areas;
        areas.reserve(cells.size());
        for (const auto& cell : cells) {
            areas.push_back(cell->state.area_preferential);
        }
        auto [mean, stddev] = get_statistics(areas);

        BOOST_CHECK_CLOSE(mean, 2, 10);
        BOOST_CHECK_CLOSE(stddev, 0.1, 10);

                   
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);

        op_diff(vertex_model);
        
        areas.clear();
        areas.reserve(cells.size());
        for (const auto& cell : cells) {
            areas.push_back(cell->state.area_preferential);
        }
        std::tie(mean, stddev) = get_statistics(areas);

        BOOST_CHECK_CLOSE(mean, 2, 10);
        BOOST_CHECK_CLOSE(stddev, 0.1, 10);


        operation(vertex_model);

        std::vector<double> areas_HCs;
        std::vector<double> areas_SCs;
        areas_HCs.reserve(cells.size());
        areas_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                areas_HCs.push_back(cell->state.area_preferential);
            }
            else if (cell->state.type == CellType::support) {
                areas_SCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs.shrink_to_fit();
        areas_SCs.shrink_to_fit();
        
        std::tie(mean, stddev) = get_statistics(areas_HCs);

        BOOST_CHECK_CLOSE(mean, 4, 10);
        BOOST_CHECK_CLOSE(stddev, 0.2, 25);
        
        std::tie(mean, stddev) = get_statistics(areas_SCs);

        BOOST_CHECK_CLOSE(mean, 5, 10);
        BOOST_CHECK_CLOSE(stddev, 0.3, 15);
        
    }
    
    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_area_adapt)
    {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "increment_area_adapt";
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);
        auto [operation, params] = build_increment_area(
            name, get_as<Config>(name, cfg), default_minim_params);

        const auto& cells = vertex_model.get_am().cells();
        double area = 0.;
        for (const auto& c : cells) {
            area += c->state.area_preferential;
        }

        op_diff(vertex_model);
        operation(vertex_model);
        
        double new_area = 0.;
        for (const auto& c : cells) {
            new_area += c->state.area_preferential;
        }

        BOOST_CHECK_CLOSE(area, new_area, 5);

        std::vector<double> areas_HCs;
        areas_HCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                areas_HCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs.shrink_to_fit();
        
        auto [mean, stddev] = get_statistics(areas_HCs);

        BOOST_CHECK_CLOSE(mean, 2, 10);
        BOOST_CHECK_CLOSE(stddev, 0.1, 25);
    }
    
    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_area_relax)
    {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "increment_area_relax";
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);
        auto [operation, params] = build_increment_area(
            name, get_as<Config>(name, cfg), default_minim_params);

        const auto& cells = vertex_model.get_am().cells();

        // apply the differentiation and operation
        op_diff(vertex_model);
        operation(vertex_model);
        
        // check that the cells fit in the domain
        double cell_area = 0.;
        for (const auto& c : cells) {
            cell_area += c->state.area_preferential;
        }

        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(domain[0] * domain[1], cell_area, 2.e-1);

        // check the statistics of HCs
        std::vector<double> areas_HCs;
        areas_HCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                areas_HCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs.shrink_to_fit();
        
        auto [mean, stddev] = get_statistics(areas_HCs);

        BOOST_CHECK_CLOSE(mean, 2, 10); // incremented by 1
        BOOST_CHECK_CLOSE(stddev, 0.1, 25);

        // check statistics of SCs
        std::vector<double> areas_SCs;
        areas_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::support) {
                areas_SCs.push_back(cell->state.area_preferential);
            }
        }
        areas_SCs.shrink_to_fit();
        
        std::tie(mean, stddev) = get_statistics(areas_SCs);

        BOOST_CHECK_CLOSE(mean, 3, 10); // incremented by 2
        BOOST_CHECK_CLOSE(stddev, 0.3, 25);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_cell_contractility) {
        const std::string name = "increment_cell_contractility";
        auto [operation, params] = build_increment_cell_contractility(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();



        std::vector<double> contractilities;
        contractilities.reserve(cells.size());
        for (const auto& cell : cells) {
            contractilities.push_back(cell->state.contractility);
        }

        // increment by one
        BOOST_TEST(   contractilities
                   == std::vector<double>(contractilities.size(), 2.));

        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);


        op_diff(vertex_model);
        
        contractilities.clear();
        contractilities.reserve(cells.size());
        for (const auto& cell : cells) {
            contractilities.push_back(cell->state.contractility);
        }

        // differentiation does not change
        BOOST_TEST(   contractilities
                   == std::vector<double>(contractilities.size(), 2.));
        

        operation(vertex_model);

        // build statistics with different cell types
        std::vector<double> contractilities_HCs;
        std::vector<double> contractilities_SCs;
        contractilities_HCs.reserve(cells.size());
        contractilities_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                contractilities_HCs.push_back(cell->state.contractility);
            }
            else if (cell->state.type == CellType::support) {
                contractilities_SCs.push_back(cell->state.contractility);
            }
        }
        contractilities_HCs.shrink_to_fit();
        contractilities_SCs.shrink_to_fit();

        // increment by 2
        BOOST_TEST(   contractilities_HCs
                   == std::vector<double>(contractilities_HCs.size(), 4.));
        // increment by 3
        BOOST_TEST(   contractilities_SCs
                   == std::vector<double>(contractilities_SCs.size(), 5.));   
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_domain)
    {
        const std::string name = "increment_domain";
        auto [operation, params] = build_increment_domain(
            name, get_as<Config>(name, cfg), default_minim_params);

        const auto& cells = vertex_model.get_am().cells();
        double area = 0.;
        for (const auto& c : cells) {
            area += c->state.area_preferential;
        }

        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(domain[0] * domain[1], area, 2.e-1);

        operation(vertex_model);
        
        double new_area = 0.;
        for (const auto& c : cells) {
            new_area += c->state.area_preferential;
        }

        BOOST_CHECK_CLOSE(area, new_area, 1e-2);

        SpaceVec new_domain = vertex_model.get_space()->get_domain_size();
        SpaceVec d_domain = new_domain - domain;

        BOOST_CHECK_CLOSE(d_domain[0], 0.1, 1e-5);
        BOOST_CHECK_CLOSE(d_domain[1], 0.2, 1e-5);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_domain_area)
    {
        const std::string name = "increment_domain_area";
        auto [operation, params] = build_increment_domain(
            name, get_as<Config>(name, cfg), default_minim_params);

        const auto& cells = vertex_model.get_am().cells();
        double area = 0.;
        for (const auto& c : cells) {
            area += c->state.area_preferential;
        }

        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(domain[0] * domain[1], area, 2.e-1);

        operation(vertex_model);
        
        double new_area = 0.;
        for (const auto& c : cells) {
            new_area += c->state.area_preferential;
        }
        BOOST_CHECK_CLOSE(area, new_area, 1e-2);

        SpaceVec new_domain = vertex_model.get_space()->get_domain_size();

        // incremented by area
        BOOST_CHECK_CLOSE(new_domain[0] * new_domain[1], 
                          domain[0] * domain[1] + 0.02, 1e-5);
        // ratio constant
        BOOST_CHECK_CLOSE(new_domain[0] / new_domain[1],
                          domain[0] / domain[1], 1.e-5);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_domain_compensate)
    {
        const std::string name = "increment_domain_compensate";
        auto [operation, params] = build_increment_domain(
            name, get_as<Config>(name, cfg), default_minim_params);

        const auto& cells = vertex_model.get_am().cells();
        double area = 0.;
        for (const auto& c : cells) {
            area += c->state.area_preferential;
        }
        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(domain[0] * domain[1], area, 2.e-1);

        operation(vertex_model);
        
        double new_area = 0.;
        for (const auto& c : cells) {
            new_area += c->state.area_preferential;
        }
        SpaceVec new_domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(new_domain[0] * new_domain[1], new_area, 2.e-1);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_domain_compensate_fix_hc)
    {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "increment_domain_compensate_fix_hc";
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);
        auto [operation, params] = build_increment_domain(
            name, get_as<Config>(name, cfg), default_minim_params);

        op_diff(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        double area = 0.;
        double area_hc = 0;
        for (const auto& c : cells) {
            area += c->state.area_preferential;
            if (c->state.type == CellType::hair) {
                area_hc += c->state.area_preferential;
            }
        }
        SpaceVec domain = vertex_model.get_space()->get_domain_size();

        operation(vertex_model);
        
        double new_area = 0.;
        double new_area_hc = 0;
        for (const auto& c : cells) {
            new_area += c->state.area_preferential;
            if (c->state.type == CellType::hair) {
                new_area_hc += c->state.area_preferential;
            }
        }
        SpaceVec new_domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(new_domain[0] * new_domain[1], new_area, 2.e-1);
        BOOST_CHECK_CLOSE(area_hc, new_area_hc, 1e-7);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_domain_compensate_fix_sc)
    {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "increment_domain_compensate_fix_sc";
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);
        auto [operation, params] = build_increment_domain(
            name, get_as<Config>(name, cfg), default_minim_params);

        op_diff(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        double area = 0.;
        double area_sc = 0;
        for (const auto& c : cells) {
            area += c->state.area_preferential;
            if (c->state.type == CellType::support) {
                area_sc += c->state.area_preferential;
            }
        }
        SpaceVec domain = vertex_model.get_space()->get_domain_size();

        operation(vertex_model);
        
        double new_area = 0.;
        double new_area_sc = 0;
        for (const auto& c : cells) {
            new_area += c->state.area_preferential;
            if (c->state.type == CellType::support) {
                new_area_sc += c->state.area_preferential;
            }
        }
        SpaceVec new_domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(new_domain[0] * new_domain[1], new_area, 2.e-1);
        BOOST_CHECK_CLOSE(area_sc, new_area_sc, 1e-7);
    }
    
    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_edge_contractility) {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "increment_edge_contractility";
        auto [operation, params] = build_increment_edge_contractility(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& am = vertex_model.get_am();
        const auto& edges = am.edges();

        for (const auto& edge : edges) {
            BOOST_TEST(edge->state.contractility == 0.11);
        }

        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);

        op_diff(vertex_model);
        // NOTE terminal differentiation, hence progenitor - non-progenitor
        //      not tested
        // NOTE The properties of the entities (edges) do not change during
        //      differentiation.

        for (const auto& edge : edges) {
            const auto [a, b] = am.adjoints_of(edge);
            if (a->state.type == CellType::hair and
                b->state.type == CellType::hair)
            {
                BOOST_TEST(edge->state.contractility == 0.11);
            }
            else if (a->state.type == CellType::support and
                     b->state.type == CellType::support)
            {
                BOOST_TEST(edge->state.contractility == 0.11);
            }
            else if (a->state.type != b->state.type)
            {
                BOOST_TEST(edge->state.contractility == 0.11);
            }
        }

        // apply operation once more
        operation(vertex_model);
        // NOTE this sets all edges contractility to the same matrix entries

        for (const auto& edge : edges) {
            const auto [a, b] = am.adjoints_of(edge);
            if (a->state.type == CellType::hair and
                b->state.type == CellType::hair)
            {
                BOOST_TEST(edge->state.contractility / 2 == 0.22);
            }
            else if (a->state.type == CellType::support and
                     b->state.type == CellType::support)
            {
                BOOST_TEST(edge->state.contractility / 2 == 0.33);
            }
            else if (a->state.type != b->state.type)
            {
                BOOST_TEST(edge->state.contractility / 2 == 0.23);
            }
        }
    }
    
    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_linetension) {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "increment_linetension";

        auto [operation, params] = build_increment_linetension(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& am = vertex_model.get_am();
        const auto& edges = am.edges();

        for (const auto& edge : edges) {
            BOOST_TEST(edge->state.linetension == 0.11);
        }

        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);

        op_diff(vertex_model);
        // NOTE terminal differentiation, hence progenitor - non-progenitor
        //      not tested
        // NOTE The properties of the entities (edges) do not change during
        //      differentiation.

        for (const auto& edge : edges) {
            const auto [a, b] = am.adjoints_of(edge);
            if (a->state.type == CellType::hair and
                b->state.type == CellType::hair)
            {
                BOOST_TEST(edge->state.linetension == 0.11);
            }
            else if (a->state.type == CellType::support and
                     b->state.type == CellType::support)
            {
                BOOST_TEST(edge->state.linetension == 0.11);
            }
            else if (a->state.type != b->state.type)
            {
                BOOST_TEST(edge->state.linetension == 0.11);
            }
        }

        // apply operation once more
        operation(vertex_model);
        // NOTE this sets all edges linetension to the same matrix entries

        for (const auto& edge : edges) {
            const auto [a, b] = am.adjoints_of(edge);
            if (a->state.type == CellType::hair and
                b->state.type == CellType::hair)
            {
                BOOST_TEST(edge->state.linetension / 2 == 0.22);
            }
            else if (a->state.type == CellType::support and
                     b->state.type == CellType::support)
            {
                BOOST_TEST(edge->state.linetension / 2 == 0.33);
            }
            else if (a->state.type != b->state.type)
            {
                BOOST_TEST(edge->state.linetension / 2 == 0.23);
            }
        }
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_increment_shape_index) {
        const std::string name = "increment_shape_index";
        auto [operation, params] = build_increment_shape_index(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();

        std::vector<double> shape_indices;
        shape_indices.reserve(cells.size());
        for (const auto& cell : cells) {
            shape_indices.push_back(cell->state.shape_index_preferential);
        }
        
        auto [mean, stddev] = get_statistics(shape_indices);
        BOOST_CHECK_CLOSE(mean, 4.6, 10); // initially is 3.6
        BOOST_CHECK_SMALL(stddev, 1e-6);

        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);


        op_diff(vertex_model);
        
        shape_indices.clear();
        shape_indices.reserve(cells.size());
        for (const auto& cell : cells) {
            shape_indices.push_back(cell->state.shape_index_preferential);
        }

        std::tie(mean, stddev) = get_statistics(shape_indices);
        BOOST_CHECK_CLOSE(mean, 4.6, 10); // differentiation does not change
        BOOST_CHECK_SMALL(stddev, 1e-6);
        

        operation(vertex_model);

        std::vector<double> shape_indices_HCs;
        std::vector<double> shape_indices_SCs;
        shape_indices_HCs.reserve(cells.size());
        shape_indices_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                shape_indices_HCs.push_back(cell->state.shape_index_preferential);
            }
            else if (cell->state.type == CellType::support) {
                shape_indices_SCs.push_back(cell->state.shape_index_preferential);
            }
        }
        shape_indices_HCs.shrink_to_fit();
        shape_indices_SCs.shrink_to_fit();
        
        std::tie(mean, stddev) = get_statistics(shape_indices_HCs);
        BOOST_CHECK_CLOSE(mean, 6.6, 10); // increment from 4.6
        BOOST_CHECK_SMALL(stddev, 1e-6);
        
        std::tie(mean, stddev) = get_statistics(shape_indices_SCs);
        BOOST_CHECK_CLOSE(mean, 7.6, 10); // increment from 4.6
        BOOST_CHECK_SMALL(stddev, 1e-6);        
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_proliferate)
    {
        const std::string name = "proliferate";
        auto [operation, params] = build_proliferate(
            name, get_as<Config>(name, cfg), default_minim_params);

        const auto time = vertex_model.get_time();
        const auto& cells = vertex_model.get_am().cells();
        const auto num_cells = cells.size();

        operation(vertex_model);

        BOOST_TEST(vertex_model.get_time() > time);
        BOOST_TEST(cells.size() == num_cells + 1);
        // NOTE the procedure of division itself is tested in PCPVertex
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_relax_area)
    {        
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "relax_area";
        auto [operation, params] = build_relax_area(
            name, get_as<Config>(name, cfg), default_minim_params);

        const auto& am = vertex_model.get_am();
        const auto& cells = vertex_model.get_am().cells();

        // random cell for that A_0 is different from A
        // A_0 / A = 2
        auto cell = cells[0];
        auto cell_hair = cells[1]; 
        cell_hair->state.type = CellType::hair;
        auto cell_support = cells[2];
        cell_support->state.type = CellType::support;

        cell->state.area_preferential = 2. * am.area_of(cell);
        cell_hair->state.area_preferential = 2. * am.area_of(cell_hair);
        cell_support->state.area_preferential = 2. * am.area_of(cell_support);

        operation(vertex_model);

        // A_0' / A = 2 - 0.1 * (A_0 / A - 1) = 1.9
        BOOST_CHECK_CLOSE(cell->state.area_preferential / am.area_of(cell), 1.9,
                          1e-8);
        BOOST_CHECK_CLOSE((cell_hair->state.area_preferential
                           / am.area_of(cell_hair)),
                          1.9, 1e-8);
        BOOST_CHECK_CLOSE((  cell_support->state.area_preferential
                           / am.area_of(cell_support)),
                          1.9, 1e-8);

        
        // use minimum area preferential
        cell->state.area_preferential = 100.;
        cell_hair->state.area_preferential = 100.;
        cell_support->state.area_preferential = 100.;

        const std::string name_min = "relax_area_minimum";
        auto [operation_min, params_min] = build_relax_area(
            name_min, get_as<Config>(name_min, cfg), default_minim_params);

        operation_min(vertex_model);

        // relax instantly, but only to minimum
        BOOST_CHECK_CLOSE(cell->state.area_preferential, 12., 1e-8);
        BOOST_CHECK_CLOSE(cell_hair->state.area_preferential, 18., 1e-8);
        BOOST_CHECK_CLOSE(cell_support->state.area_preferential, 16., 1e-8);
        
        // use maximum area preferential
        cell->state.area_preferential = 0.001;
        cell_hair->state.area_preferential = 0.001;
        cell_support->state.area_preferential = 0.001;

        const std::string name_max = "relax_area_maximum";
        auto [operation_max, params_max] = build_relax_area(
            name_max, get_as<Config>(name_max, cfg), default_minim_params);

        operation_max(vertex_model);

        // relax instantly, but only to maximum
        BOOST_CHECK_CLOSE(cell->state.area_preferential, 0.4, 1e-8);
        BOOST_CHECK_CLOSE(cell_hair->state.area_preferential, 0.3, 1e-8);
        BOOST_CHECK_CLOSE(cell_support->state.area_preferential, 0.2, 1e-8);
    }

    BOOST_AUTO_TEST_CASE(test_PCPTopology_set_area) {
        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        double init_area = domain[0] * domain[1];

        std::string name = "set_area";
        auto [operation, params] = build_set_area(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();

        std::vector<double> areas;
        areas.reserve(cells.size());
        for (const auto& cell : cells) {
            areas.push_back(cell->state.area_preferential);
        }
        auto [mean, stddev] = get_statistics(areas);

        BOOST_CHECK_CLOSE(mean, 2, 10);
        BOOST_CHECK_CLOSE(stddev, 0.1, 15);


        // differentiate progenitor cells, test again                   
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);

        op_diff(vertex_model);
        operation(vertex_model);

        std::vector<double> areas_HCs;
        std::vector<double> areas_SCs;
        areas_HCs.reserve(cells.size());
        areas_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
        auto cell_hair = cells[0];
            if (cell->state.type == CellType::hair) {
                areas_HCs.push_back(cell->state.area_preferential);
            }
            else if (cell->state.type == CellType::support) {
                areas_SCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs.shrink_to_fit();
        areas_SCs.shrink_to_fit();
        
        std::tie(mean, stddev) = get_statistics(areas_HCs);
        BOOST_CHECK_CLOSE(mean, 4, 10);
        BOOST_CHECK_CLOSE(stddev, 0.2, 10);
        
        std::tie(mean, stddev) = get_statistics(areas_SCs);
        BOOST_CHECK_CLOSE(mean, 5, 10);
        BOOST_CHECK_CLOSE(stddev, 0.3, 15);

        domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(domain[0] * domain[1], init_area, 1.e-5);


        // change only the area of support cells
        name = "set_area_partial";
        auto [operation_partial, params_partial] = build_set_area(
            name, get_as<Config>(name, cfg), default_minim_params);
       
        operation_partial(vertex_model);

        std::vector<double> areas_HCs_partial;
        areas_SCs.clear();
        areas_HCs_partial.reserve(cells.size());
        areas_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                areas_HCs_partial.push_back(cell->state.area_preferential);
            }
            else if (cell->state.type == CellType::support) {
                areas_SCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs_partial.shrink_to_fit();
        areas_SCs.shrink_to_fit();
        
        std::tie(mean, stddev) = get_statistics(areas_HCs_partial);
        BOOST_CHECK_CLOSE(mean, 4, 10);
        BOOST_CHECK_CLOSE(stddev, 0.2, 10);
        BOOST_CHECK_EQUAL_COLLECTIONS(
            areas_HCs_partial.begin(), areas_HCs_partial.end(), 
            areas_HCs.begin(), areas_HCs.end());
        
        std::tie(mean, stddev) = get_statistics(areas_SCs);
        BOOST_CHECK_CLOSE(mean, 6, 10);
        BOOST_CHECK_CLOSE(stddev, 0.2, 15);

        domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(domain[0] * domain[1], init_area, 1.e-5);


        // set area of support 
        name = "set_area_fixed";
        auto [operation_fixed, params_fixed] = build_set_area(
            name, get_as<Config>(name, cfg), default_minim_params);
       
        operation_fixed(vertex_model);

        areas_HCs.clear();
        areas_SCs.clear();
        areas_HCs.reserve(cells.size());
        areas_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                areas_HCs.push_back(cell->state.area_preferential);
            }
            else if (cell->state.type == CellType::support) {
                areas_SCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs.shrink_to_fit();
        areas_SCs.shrink_to_fit();
        
        std::tie(mean, stddev) = get_statistics(areas_HCs_partial);
        BOOST_TEST(areas_HCs == std::vector<double>(areas_HCs.size(), 2.),
                   boost::test_tools::per_element());
        BOOST_TEST(areas_SCs == std::vector<double>(areas_SCs.size(), 7.),
                   boost::test_tools::per_element());

        domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(domain[0] * domain[1], init_area, 1.e-5);
        
    }
    
    BOOST_AUTO_TEST_CASE(test_PCPTopology_set_area_relax)
    {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "set_area_relax";
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);
        auto [operation, params] = build_set_area(
            name, get_as<Config>(name, cfg), default_minim_params);

        SpaceVec domain_0 = vertex_model.get_space()->get_domain_size();

        op_diff(vertex_model);
        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        std::vector<double> areas_HCs;
        std::vector<double> areas_SCs;
        areas_HCs.reserve(cells.size());
        areas_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                areas_HCs.push_back(cell->state.area_preferential);
            }
            else if (cell->state.type == CellType::support) {
                areas_SCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs.shrink_to_fit();
        areas_SCs.shrink_to_fit();
        
        auto [mean, stddev] = get_statistics(areas_HCs);

        BOOST_CHECK_CLOSE(mean, 4, 10);
        BOOST_CHECK_CLOSE(stddev, 0.2, 25);
        
        std::tie(mean, stddev) = get_statistics(areas_SCs);

        BOOST_CHECK_CLOSE(mean, 5, 10);
        BOOST_CHECK_CLOSE(stddev, 0.3, 15);

        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        double cell_area = std::accumulate(areas_HCs.begin(),
                                           areas_HCs.end(), 0.);
        cell_area = std::accumulate(areas_SCs.begin(), areas_SCs.end(),
                                    cell_area);
        BOOST_CHECK_CLOSE(domain[0] * domain[1], cell_area, 1.e-5);

        // Test no shear implied
        BOOST_TEST(domain[0] / domain[1] == domain_0[0] / domain_0[1]);
    }
    
    BOOST_AUTO_TEST_CASE(test_PCPTopology_set_area_relax_PD_axis)
    {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "set_area_relax_PD_axis";
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);
        auto [operation, params] = build_set_area(
            name, get_as<Config>(name, cfg), default_minim_params);

        SpaceVec domain_0 = vertex_model.get_space()->get_domain_size();

        op_diff(vertex_model);
        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        std::vector<double> areas_HCs;
        std::vector<double> areas_SCs;
        areas_HCs.reserve(cells.size());
        areas_SCs.reserve(cells.size());
        for (const auto& cell : cells) {
            if (cell->state.type == CellType::hair) {
                areas_HCs.push_back(cell->state.area_preferential);
            }
            else if (cell->state.type == CellType::support) {
                areas_SCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs.shrink_to_fit();
        areas_SCs.shrink_to_fit();
        
        auto [mean, stddev] = get_statistics(areas_HCs);

        BOOST_CHECK_CLOSE(mean, 4, 10);
        BOOST_CHECK_CLOSE(stddev, 0.2, 25);
        
        std::tie(mean, stddev) = get_statistics(areas_SCs);

        BOOST_CHECK_CLOSE(mean, 5, 10);
        BOOST_CHECK_CLOSE(stddev, 0.3, 15);

        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        double cell_area = std::accumulate(areas_HCs.begin(),
                                           areas_HCs.end(), 0.);
        cell_area = std::accumulate(areas_SCs.begin(), areas_SCs.end(),
                                    cell_area);
        BOOST_CHECK_CLOSE(domain[0] * domain[1], cell_area, 1.e-5);

        // PD axis unchanged
        BOOST_TEST(domain[1] == domain_0[1]);
    }
    
    BOOST_AUTO_TEST_CASE(test_PCPTopology_set_area_relax_PD_axis_FAIL)
    {
        const std::string name = "set_area_relax_PD_axis_FAIL";
        BOOST_CHECK_THROW(
            build_set_area(name, get_as<Config>(name, cfg),
                           default_minim_params),
            std::invalid_argument
        );
    }
    
BOOST_AUTO_TEST_SUITE_END()
