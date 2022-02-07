#define BOOST_TEST_MODULE PCPVertexTestSpace

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "../space.hh"


using namespace Utopia::Models::PCPVertex::Space;

const double precision = 1e-12;

struct Fixture {
    Utopia::DataIO::Config cfg;

    CustomSpace<2> space;
    CustomSpace<2> space_periodic;

    using SpaceVec = typename CustomSpace<2>::SpaceVec;

    Fixture ()
    :
        cfg(YAML::LoadFile("space_test.yml")),
        space(cfg["2D"]["simple"]),
        space_periodic(cfg["2D"]["simple_periodic"])
    {}
};

// ++ Tests +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex_Space, Fixture)
    BOOST_FIXTURE_TEST_CASE(test_periodicity, Fixture)
    {
        // absolute coordinates
        space.set_domain_size({2., 1.});

        if (space.bound) {
            BOOST_CHECK_CLOSE(space.map_into_space({1.2, 0.1}).at(0), 1.2,
                              precision);
            BOOST_CHECK_CLOSE(space.map_into_space({2.2, 0.1}).at(0), 0.2,
                              precision);
        }
        else {
            BOOST_CHECK_CLOSE(space.map_into_space({1.2, 0.1}).at(0), 1.2,
                              precision);
            BOOST_CHECK_CLOSE(space.map_into_space({2.2, 0.1}).at(0), 2.2,
                              precision);
        }

        BOOST_CHECK_CLOSE(space.distance({0.1, 0.1}, {1.1, 0.1}), 1.,
                          precision);
        BOOST_CHECK_CLOSE(space.distance({0.1, 0.1}, {1.9, 0.1}), 1.8,
                          precision);
        BOOST_CHECK_CLOSE(space_periodic.distance({0.1, 0.1}, {1.9, 0.1}), 0.2,
                          precision);
    }

    BOOST_FIXTURE_TEST_CASE(test_intersection, Fixture)
    {
        SpaceVec pos_0 = {0.1, 0.1};
        SpaceVec vec_0 = {0.1, 0.1};

        SpaceVec pos_1 = {0.1, 0.5};
        SpaceVec vec_1 = {0.1, 0.};
        SpaceVec intersection;
        bool success;
        std::tie(intersection, success) = space.intersection(
            pos_0, vec_0, pos_1, vec_1, false, false);
        BOOST_CHECK_CLOSE(intersection.at(0), 0.5, precision);
        BOOST_CHECK_CLOSE(intersection.at(1), 0.5, precision);
        BOOST_TEST(success);
        
        // is outside finite length 
        std::tie(intersection, success) = space.intersection(
            pos_0, vec_0, pos_1, vec_1, true, true);
        BOOST_TEST(not success);
        
        // parallel
        vec_1 = vec_0;
        std::tie(intersection, success) = space.intersection(
            pos_0, vec_0, pos_1, vec_1, true, true);
        BOOST_TEST(not success);

        // vertical, test precision
        pos_1 = {0.5, 0.1};
        vec_1 = {0., 0.1};
        std::tie(intersection, success) = space.intersection(
            pos_0, vec_0, pos_1, vec_1, false, false);
        BOOST_CHECK_CLOSE(intersection.at(0), 0.5, 1e-10);
        BOOST_CHECK_CLOSE(intersection.at(1), 0.5, 1e-10);
        BOOST_TEST(success);

        // periodic case
        pos_0 = {0., 0.0};
        vec_0 = {0.1, 0.1};
        pos_1 = {0.9, 0.1};
        vec_1 = {0.3, 0.05};
        std::tie(intersection, success) = space_periodic.intersection(
            pos_0, vec_0, pos_1, vec_1, false, false);
        BOOST_CHECK_CLOSE(intersection.at(0), 0.14, precision);
        BOOST_CHECK_CLOSE(intersection.at(1), 0.14, precision);
        BOOST_TEST(success);
    }

    BOOST_AUTO_TEST_CASE(test_curved_periodic_bc)
    {
        // test the transformation
        space_periodic.set_curvature(0.02);
        SpaceVec pos({0.1, 0.1});
        auto [rho, theta] = space_periodic.transform_radial(pos);
        SpaceVec _pos = space_periodic.transform_cartesian(rho, theta);
        BOOST_CHECK_CLOSE(_pos[0], pos[0], 1.e-4);
        BOOST_CHECK_CLOSE(_pos[1], pos[1], 1.e-4);
    }

    BOOST_AUTO_TEST_CASE(test_skewed_periodic_bc)
    {
        space_periodic.set_domain_size(SpaceVec({10., 10.}));
        
        CustomSpace<2> space_sx(cfg["2D"]["simple_periodic"]);
        space_sx.set_domain_size(SpaceVec({10., 10.}));
        space_sx.set_skew(SpaceVec({1., 0.}));

        CustomSpace<2> space_sy(cfg["2D"]["simple_periodic"]);
        space_sy.set_domain_size(SpaceVec({10., 10.}));
        space_sy.set_skew(SpaceVec({0., 1.}));
        

        // right boundary
        SpaceVec pos({2.2, 12.2});
        SpaceVec _pos = space_periodic.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        _pos = space_sx.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 1.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        _pos = space_sy.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        // left boundary
        pos = SpaceVec({2.2, -7.8});
        _pos = space_periodic.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        _pos = space_sx.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 3.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        _pos = space_sy.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);
        
        // upper boundary
        pos = SpaceVec({12.2, 2.2});
        _pos = space_periodic.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        _pos = space_sx.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        _pos = space_sy.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 1.2, 1.e-6);

        // lower boundary
        pos = SpaceVec({-7.8, 2.2});
        _pos = space_periodic.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        _pos = space_sx.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);

        _pos = space_sy.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 3.2, 1.e-6);

        // double boundary
        pos = SpaceVec({12.2, 12.2});
        _pos = space_periodic.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);
        _pos = space_sx.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 1.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 2.2, 1.e-6);
        _pos = space_sy.map_into_space(pos);
        BOOST_CHECK_CLOSE(_pos[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(_pos[1], 1.2, 1.e-6);


        // displacements
        // left boundary
        SpaceVec pos_0({1.2, 2.2});
        SpaceVec pos_1({8.8, 3.4});
        SpaceVec displ = space_periodic.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], -2.4, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 1.2, 1.e-6);
        displ = space_sx.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], -2.4, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 1.2, 1.e-6);
        displ = space_sy.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], -2.4, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 0.2, 1.e-6);

        // right boundary
        pos_0 = SpaceVec({8.8, 2.2});
        pos_1 = SpaceVec({1.2, 3.4});
        displ = space_periodic.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 2.4, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 1.2, 1.e-6);
        displ = space_sx.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 2.4, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 1.2, 1.e-6);
        displ = space_sy.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 2.4, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 2.2, 1.e-6);

        // upper boundary
        pos_0 = SpaceVec({1.2, 2.2});
        pos_1 = SpaceVec({2.4, 8.8});
        displ = space_periodic.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 1.2, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], -3.4, 1.e-6);
        displ = space_sx.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 0.2, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], -3.4, 1.e-6);
        displ = space_sy.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 1.2, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], -3.4, 1.e-6);

        // lower boundary
        pos_0 = SpaceVec({1.2, 8.8});
        pos_1 = SpaceVec({2.4, 2.2});
        displ = space_periodic.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 1.2, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 3.4, 1.e-6);
        displ = space_sx.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 2.2, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 3.4, 1.e-6);
        displ = space_sy.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 1.2, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 3.4, 1.e-6);

        // upper right boundary
        pos_0 = SpaceVec({9.2, 9.2});
        pos_1 = SpaceVec({0.8, 0.8});
        displ = space_periodic.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 1.6, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 1.6, 1.e-6);
        displ = space_sx.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 2.6, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 1.6, 1.e-6);
        displ = space_sy.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], 1.6, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], 2.6, 1.e-6);

        std::swap(pos_0, pos_1);
        displ = space_periodic.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], -1.6, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], -1.6, 1.e-6);
        displ = space_sx.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], -2.6, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], -1.6, 1.e-6);
        displ = space_sy.displacement(pos_0, pos_1);
        BOOST_CHECK_CLOSE(displ[0], -1.6, 1.e-6);
        BOOST_CHECK_CLOSE(displ[1], -2.6, 1.e-6);      
    }

BOOST_AUTO_TEST_SUITE_END()
