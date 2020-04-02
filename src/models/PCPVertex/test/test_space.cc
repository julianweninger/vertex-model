#define BOOST_TEST_MODULE PCPVertexTestSpace

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "utils.hh"
#include "../space.hh"


using namespace Utopia;
using namespace Utopia::Models::PCPVertex::Space;

using CSpace = CustomSpace<2>;
using SpaceVec = CSpace::SpaceVec;
using ACSpace = AbsoluteCustomSpace<2>;
using ASpaceVec = ACSpace::SpaceVec;

const double precision = 1e-12;

// ++ Tests +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
BOOST_FIXTURE_TEST_SUITE (test_Space, ModelFixture)

BOOST_AUTO_TEST_CASE(CustomSpace)
{
    CSpace space;

    BOOST_TEST(space.contains<false>({0.1, 0.2}));
    BOOST_CHECK_CLOSE(space.map_into_space({2.1, 0.2}).at(0), 0.1, precision);
    BOOST_CHECK_CLOSE(space.displacement({0.1, 0.1}, {0.3, 0.1}).at(0), 0.2,
                      precision);
    BOOST_CHECK_CLOSE(space.displacement({0.1, 0.1}, {0.9, 0.1}).at(0), -0.2,
                      precision);
    
    BOOST_CHECK_CLOSE(space.distance({0.1, 0.1}, {0.3, 0.1}), 0.2, precision);
    
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
    vec_1 = {0.3, 0.};
    std::tie(intersection, success) = space.intersection(pos_0, vec_0, 
                                                         pos_1, vec_1,
                                                         false, false);
    BOOST_CHECK_CLOSE(intersection.at(0), 0.1, precision);
    BOOST_CHECK_CLOSE(intersection.at(1), 0.1, precision);
    BOOST_TEST(success);
}


BOOST_AUTO_TEST_CASE(AbsoluteCustomSpace)
{
    ACSpace space;
    space.domain_size = {2., 1.};

    BOOST_CHECK_CLOSE(space.convert_absolute({0.1, 0.1}).at(0), 0.2, precision);

    BOOST_CHECK_CLOSE(space.distance({0.1, 0.1}, {0.3, 0.1}), 0.4, precision);
}

BOOST_AUTO_TEST_SUITE_END()