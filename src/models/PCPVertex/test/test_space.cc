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

BOOST_FIXTURE_TEST_CASE(test_CustomSpace, Fixture)
{
    SpaceVec pos_0 = {0.1, 0.1};
    SpaceVec vec_0 = {0.1, 0.1};

    SpaceVec pos_1 = {0.1, 0.5};
    SpaceVec vec_1 = {0.1, 0.};
    auto [intersection, success] = space.intersection(pos_0, vec_0, 
                                                      pos_1, vec_1,
                                                      false, false);
    BOOST_CHECK_CLOSE(intersection.at(0), 0.5, precision);
    BOOST_CHECK_CLOSE(intersection.at(1), 0.5, precision);
    BOOST_TEST(success);
    
    // is outside finite length 
    std::tie(intersection, success) = space.intersection(pos_0, vec_0,
                                                         pos_1, vec_1,
                                                         true, true);
    BOOST_TEST(not success);
    
    // parallel
    vec_1 = vec_0;
    std::tie(intersection, success) = space.intersection(pos_0, vec_0,
                                                         pos_1, vec_1,
                                                         true, true);
    BOOST_TEST(not success);

    // vertical, test precision
    pos_1 = {0.5, 0.1};
    vec_1 = {0., 0.1};
    std::tie(intersection, success) = space.intersection(pos_0, vec_0, 
                                                         pos_1, vec_1,
                                                         false, false);
    BOOST_CHECK_CLOSE(intersection.at(0), 0.5, 1e-10);
    BOOST_CHECK_CLOSE(intersection.at(1), 0.5, 1e-10);
    BOOST_TEST(success);

    // periodic case
    pos_0 = {0., 0.0};
    vec_0 = {0.1, 0.1};
    pos_1 = {0.9, 0.1};
    vec_1 = {0.3, 0.05};
    std::tie(intersection, success) = space_periodic.intersection(pos_0, vec_0, 
                                                                  pos_1, vec_1,
                                                                  false, false);
    BOOST_CHECK_CLOSE(intersection.at(0), 0.14, precision);
    BOOST_CHECK_CLOSE(intersection.at(1), 0.14, precision);
    BOOST_TEST(success);

    // absolute coordinates
    space.set_domain_size({2., 1.});

    BOOST_CHECK_CLOSE(space.map_to_absolute_space({0.1, 0.1}).at(0), 0.2, 
                      precision);

    BOOST_CHECK_CLOSE(space.distance({0.1, 0.1}, {0.3, 0.1}), 0.4, precision);
}