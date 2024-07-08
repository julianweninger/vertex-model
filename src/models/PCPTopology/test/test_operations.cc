#define BOOST_TEST_MODULE PCPTopologyTestOperations

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "utils.hh"

#include "../../PCPVertex/PCPVertex.hh"
#include "../../PCPVertex/algorithm.hh"
#include "../../PCPVertex/operations.hh"
#include "../../PCPVertex/PCPVertex_write_tasks.hh"
#include "../../PCPVertex/work_function.hh"

#include "../PCPTopology.hh"
#include "../operations.hh"
#include "../PCPTopology_write_tasks.hh"

#include "../../Collier/Collier.hh"
#include "../../NotchDelta/NotchDelta.hh"
#include "../../NotchDelta/Differentiation_write_tasks.hh"
#include "../../Collier/Collier_write_tasks.hh"

using namespace Utopia;
using namespace Utopia::Models::PCPTopology::OperationCollection;
using Utopia::DataIO::Config;

using VertexModel = Models::PCPVertex::PCPVertex;
using SpaceVec = VertexModel::SpaceVec;

using TopologyModel = Models::PCPTopology::PCPTopology;

using NotchModel = Models::NotchDelta::NotchDelta;


struct Fixture : ModelFixture {
    PseudoParent<> pp;
    
    VertexModel vertex_model;
    
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

BOOST_AUTO_TEST_CASE(test_PCPTopology_apply_operation) {
    Config cfg(YAML::LoadFile("test_operations.yml"));
    Utopia::Models::PCPVertex::MinimizationParams default_minim_params(
        get_as<Config>("default_minimization", cfg));

    PseudoParent<> pseudo_parent("topology_model_cfg.yml");
    TopologyModel model("PCPTopology", pseudo_parent, {},
        std::make_tuple(Models::PCPVertex::DataIO::time_energy_adaptor));
    std::string name;
    std::size_t prev_time = model.get_continuous_time();
    std::size_t time;
        
    name = "apply_operation_prolog";
    model.register_operation(build_void(
        name, get_as<Config>(name, cfg), default_minim_params));
    model.prolog();
    time = model.get_continuous_time();
    BOOST_TEST(time - prev_time == 20);

    // 5 steps per minimization
    // 2 repeats
    // jiggle + step each repeat,
    // 1 for initial condition per minimization
    // 2 x 2 + 1 = 5 iterations
    std::function<std::size_t(TopologyModel&, std::string)> apply = 
        [cfg, default_minim_params]
        (TopologyModel& model, std::string name)
    {
        std::size_t time = model.get_continuous_time();
        model.register_operation(build_void(
            name, get_as<Config>(name, cfg), default_minim_params));
        model.iterate();

        return (model.get_continuous_time() - time);
    };

    BOOST_TEST(apply(model, "apply_operation") == 25);
    // NOTE 5 iterations, minimize every

    BOOST_TEST(apply(model, "apply_operation_every") == 25);
    // NOTE 5 iterations, minimize every

    BOOST_TEST(apply(model, "apply_operation_once") == 5);
    // NOTE 5 iterations, minimize once

    BOOST_TEST(apply(model, "apply_operation_manual") == 0);
    // NOTE 5 iterations, no minimization

    BOOST_TEST(apply(model, "apply_operation_repeat") == 35);
    // NOTE 5 iterations, jiggle 3 times (2x3 + 1)
        
    name = "apply_operation_epilog";
    model.register_operation(build_void(
        name, get_as<Config>(name, cfg), default_minim_params));
    
    prev_time = model.get_continuous_time();
    model.epilog();
    
    time = model.get_continuous_time();
    BOOST_TEST(time - prev_time == 20);
}

BOOST_FIXTURE_TEST_SUITE (test_PCPTopology_operations, Fixture)

    BOOST_AUTO_TEST_CASE(test_PCPTopology_jiggle) {
        const std::string name = "void";
        auto [operation, params] = build_void(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);
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


    BOOST_AUTO_TEST_CASE(test_PCPTopology_differentiate_domain) {
        const std::string name = "differentiate_domain";
        auto [operation, params] = build_differentiate_domain(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();
        auto cnt_0 = std::count_if(cells.begin(), cells.end(),
                                   [](const auto& cell) {
                                        return cell->state.type == 0; });
        auto cnt_1 = std::count_if(cells.begin(), cells.end(),
                                   [](const auto& cell) {
                                        return cell->state.type == 1; });
        auto cnt_2 = std::count_if(cells.begin(), cells.end(),
                                   [](const auto& cell) {
                                        return cell->state.type == 2; });
        auto cnt_3 = std::count_if(cells.begin(), cells.end(),
                                   [](const auto& cell) {
                                        return cell->state.type == 3; });
        auto cnt_4 = std::count_if(cells.begin(), cells.end(),
                                   [](const auto& cell) {
                                        return cell->state.type == 4; });

        double f_accu = 0; // size of largest square containing type-1
        double f_area = 0.2 * 0.3 * 64.; // size of largest square containing value
        BOOST_CHECK_SMALL(cnt_0 - f_area + f_accu, 2.);
        f_accu = 0.2 * 0.3 * 64.;
        f_area = 0.6 * 0.8 * 64.;
        BOOST_CHECK_SMALL(cnt_1 - f_area + f_accu, 2.);
        f_accu = 0.6 * 0.8 * 64.;
        f_area = 0.8 * 1.0 * 64.;
        BOOST_CHECK_SMALL(cnt_2 - f_area + f_accu, 2.);
        f_accu = 0.8 * 1.0 * 64.;
        f_area = 1.0 * 1.0 * 64.;
        BOOST_CHECK_SMALL(cnt_3 - f_area + f_accu, 2.);
        BOOST_CHECK_EQUAL(cnt_4, 0);
    }


    BOOST_AUTO_TEST_CASE(test_PCPTopology_differentiate_NotchDelta) {
        // Function to compare the correct mapping of cells
        std::function<void(const VertexModel&)> check_mapping = []
                (const VertexModel& vertex_model)
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
                BOOST_TEST(   nd_cell->state.cell_type
                           == NotchModel::CellType(type));
            }
        };

        std::string name = "differentiate_NotchDelta";

        std::shared_ptr<NotchModel> notch_delta(new NotchModel(
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

    // Test the NotchModel differentiation when NotchInitial condition is 
    // overwritten
    BOOST_AUTO_TEST_CASE(test_PCPTopology_differentiate_NotchDelta_skip_init) {
        // Function to compare the correct mapping of cells
        std::function<void(const VertexModel&)> check_mapping = []
                (const VertexModel& vertex_model)
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
                BOOST_TEST(   nd_cell->state.cell_type
                           == NotchModel::CellType(type));
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
        auto cnt_prog = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return (   cell->state.type
                        == int(NotchModel::CellType::progenitor)); }
        );
        BOOST_CHECK_SMALL(double(cnt_prog), 1.e-6);

        auto cnt_support = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return (   cell->state.type
                        == int(NotchModel::CellType::support)); }
        );
        BOOST_CHECK_SMALL(cnt_support - 0.65 * cells.size(), 1.);

        auto cnt_hair = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return (   cell->state.type
                        == int(NotchModel::CellType::hair)); }
        );
        BOOST_CHECK_SMALL(cnt_hair - 0.35 * cells.size(), 1.);

        // save states for next try
        std::vector<std::size_t> cell_states;
        cell_states.reserve(cells.size());
        for (const auto& cell : cells) {
            cell_states.push_back(cell->state.type);
        }


        // The NotchModel
        // the notch differentiation
        std::shared_ptr<NotchModel> notch_delta(new NotchModel(
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

        auto new_cnt_prog = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return (   cell->state.type
                        == int(NotchModel::CellType::progenitor)); }
        );
        BOOST_CHECK_SMALL(double(new_cnt_prog), 1.e-6);

        auto new_cnt_support = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return (   cell->state.type
                        == int(NotchModel::CellType::support)); }
        );
        BOOST_CHECK_EQUAL(cnt_support, new_cnt_support);

        auto new_cnt_hair = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return (   cell->state.type
                        == int(NotchModel::CellType::hair)); }
        );
        BOOST_CHECK_EQUAL(cnt_hair, new_cnt_hair);
    }


    BOOST_AUTO_TEST_CASE(test_PCPTopology__differentiate_random) {
        const std::string name = "differentiate_random";
        auto [operation, params] = build_differentiate_random(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        auto cnt_0 = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.type == 0;
            }
        );
        auto cnt_1 = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.type == 1;
            }
        );
        auto cnt_2 = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.type == 2;
            }
        );
        auto cnt_3 = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.type > 2;
            }
        );

        BOOST_CHECK_SMALL(cnt_0 - 0.2 * 64, 2.);
        BOOST_CHECK_SMALL(cnt_1 - 0.5 * 64, 2.);
        BOOST_CHECK_SMALL(cnt_2 - 0.3 * 64, 2.);
        BOOST_CHECK_EQUAL(cnt_3, 0);
    }


    BOOST_AUTO_TEST_CASE(test_PCPTopology__differentiate_random_fraction) {
        const std::string name = "differentiate_random_fraction";
        auto [operation, params] = build_differentiate_random(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();
        auto cnt_0 = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.type == 0;
            }
        );
        auto cnt_1 = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.type == 1;
            }
        );
        auto cnt_2 = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.type == 2;
            }
        );
        auto cnt_3 = std::count_if(
            cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->state.type > 2;
            }
        );

        BOOST_CHECK_EQUAL(cnt_0, std::ceil(0.2 * 64));
        BOOST_CHECK_EQUAL(cnt_1, std::ceil(0.5 * 64));
        BOOST_CHECK_EQUAL(cnt_2, 64 - cnt_1 - cnt_0);
        BOOST_CHECK_EQUAL(cnt_3, 0);
    }


    BOOST_AUTO_TEST_CASE(test_PCPTopology_minimize_cell_contacts) {
        const std::string name = "minimize_cell_contacts";
        auto config = get_as<Config>(name, cfg);
        auto [operation, params] = build_minimize_cell_contacts(
            name, config, default_minim_params);

        auto [differentiate, _params] = build_differentiate_random(
            name, get_as<Config>("differentiate_vanish", config),
            default_minim_params
        );
    	
        differentiate(vertex_model);

        operation(vertex_model);

        const auto& am = vertex_model.get_am();
        const auto& cells = am.cells();
        std::size_t contacts = 0;
        for (const auto& cell : cells) {
            if (cell->state.type != 1) {
                continue;
            }

            for (const auto& n : am.neighbors_of(cell)) {
                if (n->state.type == cell->state.type) {
                    contacts++;
                }
            }
        }
        contacts /= 2;

        BOOST_CHECK_SMALL(static_cast<double>(contacts), 1.1);


        for (const auto& cell : cells) {
            cell->state.type = 0;
        }

        auto [_differentiate, __params] = build_differentiate_random(
            name, get_as<Config>("differentiate_minimize", config),
            default_minim_params
        );

        _differentiate(vertex_model);

        contacts = 0;
        for (const auto& cell : cells) {
            if (cell->state.type != 1) {
                continue;
            }

            for (const auto& n : am.neighbors_of(cell)) {
                if (n->state.type == cell->state.type) {
                    contacts++;
                }
            }
        }

        // in 50 / 50 expect 3 / 6 neighbours to be type 1 in ranodm
        BOOST_CHECK_SMALL(static_cast<double>(contacts), 2.75 * 32);
        // NOTE have 52 / 48 ratio

        operation(vertex_model);

        std::size_t new_contacts = 0;
        for (const auto& cell : cells) {
            if (cell->state.type != 1) {
                continue;
            }

            for (const auto& n : am.neighbors_of(cell)) {
                if (n->state.type == cell->state.type) {
                    new_contacts++;
                }
            }
        }
        // in 50 / 50 expect 2 / 6 neighbours to be type 1 at best
        BOOST_CHECK_SMALL(static_cast<double>(new_contacts), 1.8 * 32);
    }


    BOOST_AUTO_TEST_CASE(test_PCPTopology_proliferate)
    {
        const std::string name = "proliferate";
        auto [operation, params] = build_proliferate(
            name, get_as<Config>(name, cfg), default_minim_params,
            vertex_model.get_logger());

        auto domain_0 = vertex_model.get_space()->get_domain_size();
        double area_0 = domain_0[0] * domain_0[1];

        const auto time = vertex_model.get_time();
        const auto& cells = vertex_model.get_am().cells();
        const auto num_cells = cells.size();

        operation(vertex_model);

        BOOST_TEST(vertex_model.get_time() > time);
        BOOST_TEST(cells.size() == num_cells + 1);
        // NOTE the procedure of division itself is tested in VertexModel


        auto domain = vertex_model.get_space()->get_domain_size();
        double area = domain[0] * domain[1];

        if (vertex_model.get_space()->periodic) {
            BOOST_CHECK_CLOSE(area - area_0, 1., 1.e-5);
        }
    }
    
BOOST_AUTO_TEST_SUITE_END()
