#define BOOST_TEST_MODULE PCPVertexTest

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include "utils.hh"

bool test_equal (double a, double b, double precision = 1e-14) {
    return std::abs(a - b) < precision;
}

template<bool periodic_bc>
void test_initialisation_hexagon (std::string cfg)
{
    auto model = model_factory<periodic_bc>(cfg);
    auto model_cfg = model.get_cfg();

    auto vertices = model.get_vertices();
    auto edges = model.get_edges();
    auto cells = model.get_cells();

    int num_rows = get_as<int>("lattice_rows", model_cfg);
    int num_columns = get_as<int>("lattice_columns", model_cfg);
    BOOST_TEST (cells.size() == num_rows * num_columns);

    if constexpr (periodic_bc) {
        for (auto v : vertices) {
            BOOST_TEST (v.lock()->adj_edges.size() == 3);
            BOOST_TEST (v.lock()->adj_cells.size() == 3);
        }
        for (auto e : edges) {
            BOOST_TEST (not e.lock()->adj_cell_a.expired());
            BOOST_TEST (not e.lock()->adj_cell_b.expired());
        }
        for (auto c : cells) {
            BOOST_TEST (c.lock()->vertices.size() == 6);
            BOOST_TEST (c.lock()->edges_ordered.size() == 6);
        }
    }
    else {
        for (auto v : vertices) {
            BOOST_TEST (v.lock()->adj_edges.size() <= 3);
            BOOST_TEST (v.lock()->adj_cells.size() <= 3);
        }
        for (auto e : edges) {
            BOOST_TEST ((not e.lock()->adj_cell_a.expired() or 
                         not e.lock()->adj_cell_b.expired()) == true);
        }
        for (auto c : cells) {
            BOOST_TEST (c.lock()->vertices.size() == 6);
            BOOST_TEST (c.lock()->edges_ordered.size() == 6);
        }
    }

    destruct_model_factory(model);
}

BOOST_AUTO_TEST_CASE(Initialisation_hexagon_periodic)
{
    test_initialisation_hexagon<true>("test_periodic.yml");
}

BOOST_AUTO_TEST_CASE(Initialisation_hexagon_non_periodic)
{
    test_initialisation_hexagon<false>("test_non_periodic.yml");
}


template<bool periodic_bc>
void test_topological_increase_domain (std::string cfg) {
    auto model = model_factory<periodic_bc>(cfg);
    auto [Lx, Ly] = model.get_domain_size();

    double dA = 0.1;

    model.increase_domain_size(dA);
    auto [Lx_prime, Ly_prime] = model.get_domain_size();
    BOOST_TEST (test_equal(Lx_prime * Ly_prime, Lx * Ly + dA));
    BOOST_TEST (test_equal(Lx_prime / Ly_prime, Lx / Ly));

    model.increase_domain_size(-dA);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx_prime * Ly_prime, Lx * Ly));
    BOOST_TEST (test_equal(Lx_prime / Ly_prime, Lx / Ly));

    destruct_model_factory(model);
}

BOOST_AUTO_TEST_CASE(Topological_increase_domain_periodic)
{
    test_topological_increase_domain<true>("test_periodic.yml");
}

BOOST_AUTO_TEST_CASE(Topological_increase_domain_non_periodic)
{
    test_topological_increase_domain<false>("test_non_periodic.yml");
}

template<bool periodic_bc>
void test_topological_stretch (std::string cfg) {
    auto model = model_factory<periodic_bc>(cfg);
    auto model_cfg = model.get_cfg();
    auto cells = model.get_cells();

    auto [Lx, Ly] = model.get_domain_size();
    double Lx_prime, Ly_prime;

    double dx = 0.1;
    double dy = 0.2;

    double area_preferential = cells[0].lock()->area_preferential;

    // without compensation
    model.stretch_domain(dx, 0., false);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx + dx, Lx_prime));
    BOOST_TEST (test_equal(Ly, Ly_prime));
    model.stretch_domain(-dx, 0., false);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx, Lx_prime));
    BOOST_TEST (test_equal(Ly, Ly_prime));
    BOOST_TEST (test_equal(area_preferential, cells[0].lock()->area_preferential));

    model.stretch_domain(0., dy, false);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx, Lx_prime));
    BOOST_TEST (test_equal(Ly + dy, Ly_prime));
    model.stretch_domain(0., -dy, false);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx, Lx_prime));
    BOOST_TEST (test_equal(Ly, Ly_prime));
    BOOST_TEST (test_equal(area_preferential, cells[0].lock()->area_preferential));

    // with compensation
    model.stretch_domain(dx, 0., true);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx + dx, Lx_prime));
    BOOST_TEST (test_equal(Ly, Ly_prime));
    BOOST_TEST (area_preferential < cells[0].lock()->area_preferential);
    model.stretch_domain(-dx, 0., true);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx, Lx_prime));
    BOOST_TEST (test_equal(Ly, Ly_prime));
    BOOST_TEST (test_equal(area_preferential, cells[0].lock()->area_preferential));

    model.stretch_domain(0., dy, true);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx, Lx_prime));
    BOOST_TEST (test_equal(Ly + dy, Ly_prime));
    BOOST_TEST (area_preferential < cells[0].lock()->area_preferential);
    model.stretch_domain(0., -dy, true);
    std::tie(Lx_prime, Ly_prime) = model.get_domain_size();
    BOOST_TEST (test_equal(Lx, Lx_prime));
    BOOST_TEST (test_equal(Ly, Ly_prime));
    BOOST_TEST (test_equal(area_preferential, cells[0].lock()->area_preferential));

    destruct_model_factory(model);
}

BOOST_AUTO_TEST_CASE(Topological_stretch_periodic)
{
    test_topological_stretch<true>("test_periodic.yml");
}

BOOST_AUTO_TEST_CASE(Topological_stretch_non_periodic)
{
    test_topological_stretch<false>("test_non_periodic.yml");
}