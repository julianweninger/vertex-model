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
    double stdev = std::sqrt(sq_sum / values.size() - mean * mean);

    return std::make_pair(mean, stdev);
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
        using Utopia::Models::NotchDelta::NotchDelta;

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

        const auto& cells = vertex_model.get_am().cells();
        auto cnt = std::count_if(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->custom_links().nd_cell; });
        BOOST_TEST(cnt == cells.size());

        // Continue iteration
        operation(vertex_model);
        BOOST_TEST(notch_delta->get_time() == 200);

        cnt = std::count_if(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->custom_links().nd_cell; });
        BOOST_TEST(cnt == cells.size());
        
        
        // Continue iteration with new task
        name = "differentiate_NotchDelta_second_task";
        auto [new_operation, new_params] = build_differentiate_NotchDelta(
            name, get_as<Config>(name, cfg), default_minim_params,
            notch_delta, notch_delta_prolog);
        
        new_operation(vertex_model);
        
        BOOST_TEST(notch_delta->get_time() == 210);

        cnt = std::count_if(cells.begin(), cells.end(),
            [](const auto& cell) {
                return cell->custom_links().nd_cell; });
        BOOST_TEST(cnt == cells.size());
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
        auto [mean, stdev] = get_statistics(areas);

        BOOST_CHECK_CLOSE(mean, 2, 10);
        BOOST_CHECK_CLOSE(stdev, 0.1, 10);

                   
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);

        op_diff(vertex_model);
        
        areas.clear();
        areas.reserve(cells.size());
        for (const auto& cell : cells) {
            areas.push_back(cell->state.area_preferential);
        }
        std::tie(mean, stdev) = get_statistics(areas);

        BOOST_CHECK_CLOSE(mean, 2, 10);
        BOOST_CHECK_CLOSE(stdev, 0.1, 10);


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
        
        std::tie(mean, stdev) = get_statistics(areas_HCs);

        BOOST_CHECK_CLOSE(mean, 4, 10);
        BOOST_CHECK_CLOSE(stdev, 0.2, 25);
        
        std::tie(mean, stdev) = get_statistics(areas_SCs);

        BOOST_CHECK_CLOSE(mean, 5, 10);
        BOOST_CHECK_CLOSE(stdev, 0.3, 15);
        
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
        
        auto [mean, stdev] = get_statistics(areas_HCs);

        BOOST_CHECK_CLOSE(mean, 2, 10);
        BOOST_CHECK_CLOSE(stdev, 0.1, 25);
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

        SpaceVec d_domain = new_domain - domain;
        BOOST_CHECK_CLOSE(d_domain[0], 0.1, 1e-5);
        BOOST_CHECK_CLOSE(d_domain[1], 0.2, 1e-5);
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

        SpaceVec d_domain = new_domain - domain;
        BOOST_CHECK_CLOSE(d_domain[0], 0.1, 1e-5);
        BOOST_CHECK_CLOSE(d_domain[1], 0.2, 1e-5);
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

        SpaceVec d_domain = new_domain - domain;
        BOOST_CHECK_CLOSE(d_domain[0], 0.1, 1e-5);
        BOOST_CHECK_CLOSE(d_domain[1], 0.2, 1e-5);
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
        
        auto [mean, stdev] = get_statistics(shape_indices);
        BOOST_CHECK_CLOSE(mean, 4.6, 10); // initially is 3.6
        BOOST_CHECK_SMALL(stdev, 1e-6);

        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);


        op_diff(vertex_model);
        
        shape_indices.clear();
        shape_indices.reserve(cells.size());
        for (const auto& cell : cells) {
            shape_indices.push_back(cell->state.shape_index_preferential);
        }

        std::tie(mean, stdev) = get_statistics(shape_indices);
        BOOST_CHECK_CLOSE(mean, 4.6, 10); // differentiation does not change
        BOOST_CHECK_SMALL(stdev, 1e-6);
        

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
        
        std::tie(mean, stdev) = get_statistics(shape_indices_HCs);
        BOOST_CHECK_CLOSE(mean, 6.6, 10); // increment from 4.6
        BOOST_CHECK_SMALL(stdev, 1e-6);
        
        std::tie(mean, stdev) = get_statistics(shape_indices_SCs);
        BOOST_CHECK_CLOSE(mean, 7.6, 10); // increment from 4.6
        BOOST_CHECK_SMALL(stdev, 1e-6);        
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

    BOOST_AUTO_TEST_CASE(test_PCPTopology_set_area) {
        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        double init_area = domain[0] * domain[1];

        const std::string name = "set_area";
        auto [operation, params] = build_set_area(
            name, get_as<Config>(name, cfg), default_minim_params);

        operation(vertex_model);

        const auto& cells = vertex_model.get_am().cells();

        std::vector<double> areas;
        areas.reserve(cells.size());
        for (const auto& cell : cells) {
            areas.push_back(cell->state.area_preferential);
        }
        auto [mean, stdev] = get_statistics(areas);

        BOOST_CHECK_CLOSE(mean, 2, 10);
        BOOST_CHECK_CLOSE(stdev, 0.1, 10);


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
            if (cell->state.type == CellType::hair) {
                areas_HCs.push_back(cell->state.area_preferential);
            }
            else if (cell->state.type == CellType::support) {
                areas_SCs.push_back(cell->state.area_preferential);
            }
        }
        areas_HCs.shrink_to_fit();
        areas_SCs.shrink_to_fit();
        
        std::tie(mean, stdev) = get_statistics(areas_HCs);
        BOOST_CHECK_CLOSE(mean, 4, 10);
        BOOST_CHECK_CLOSE(stdev, 0.2, 10);
        
        std::tie(mean, stdev) = get_statistics(areas_SCs);
        BOOST_CHECK_CLOSE(mean, 5, 10);
        BOOST_CHECK_CLOSE(stdev, 0.3, 15);

        domain = vertex_model.get_space()->get_domain_size();
        BOOST_CHECK_CLOSE(domain[0] * domain[1], init_area, 1.e-5);
        
    }
    
    BOOST_AUTO_TEST_CASE(test_PCPTopology_set_area_adapt)
    {
        using CellType = Models::PCPVertex::PCPVertex::CellType;

        const std::string name = "set_area_adapt";
        auto [op_diff, params_diff] = build_differentiate_random(
            "differentiate_random", get_as<Config>("differentiate_random", cfg),
            default_minim_params);
        auto [operation, params] = build_set_area(
            name, get_as<Config>(name, cfg), default_minim_params);

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
        
        auto [mean, stdev] = get_statistics(areas_HCs);

        BOOST_CHECK_CLOSE(mean, 4, 10);
        BOOST_CHECK_CLOSE(stdev, 0.2, 25);
        
        std::tie(mean, stdev) = get_statistics(areas_SCs);

        BOOST_CHECK_CLOSE(mean, 5, 10);
        BOOST_CHECK_CLOSE(stdev, 0.3, 15);

        SpaceVec domain = vertex_model.get_space()->get_domain_size();
        double cell_area = std::accumulate(areas_HCs.begin(),
                                             areas_HCs.end(), 0.);
        cell_area += std::accumulate(areas_SCs.begin(), areas_SCs.end(), 0.);
        BOOST_CHECK_CLOSE(domain[0] * domain[1], cell_area, 1.e-5);
    }
    
BOOST_AUTO_TEST_SUITE_END()
