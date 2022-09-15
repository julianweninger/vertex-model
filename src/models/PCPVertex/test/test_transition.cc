#define BOOST_TEST_MODULE PCPVertexTestTransitions

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <numeric>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "utils.hh"
#include "../PCPVertex.hh"
#include "../algorithm.hh"
#include "../operations.hh"
#include "../PCPVertex_write_tasks.hh"

using namespace Utopia;
using namespace Utopia::Models::PCPVertex;
using Utopia::DataIO::Config;

using Vertex = Utopia::Models::PCPVertex::PCPVertex::Vertex;
using Edge = Utopia::Models::PCPVertex::PCPVertex::Edge;
using Cell = Utopia::Models::PCPVertex::PCPVertex::Cell;

using SpaceVec = Utopia::Models::PCPVertex::PCPVertex::SpaceVec;

template<bool periodic>
struct Fixture {    
    Models::PCPVertex::PCPVertex vertex_model;

    const std::shared_ptr<spdlog::logger> log;

    Fixture ()
    :
        vertex_model(model_factory()),
        log(spdlog::stdout_color_mt("TransitionTests"))
    {
        log->debug("Setting log level to '{}' ...", "trace");
        log->set_level(spdlog::level::from_str("trace"));
    }

    ~Fixture()
    {
        vertex_model.get_logger()->info("Tearing down ...");
        std::remove("test_data.h5");
        spdlog::drop_all();
    }
    
    PCPVertex model_factory() {
        using Utopia::Models::PCPVertex::DataIO::time_energy_adaptor;

        if constexpr (periodic) {
            Utopia::PseudoParent pp("test_transitions_periodic.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else {
            Utopia::PseudoParent pp("test_transitions.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
    }
};

typedef boost::mpl::vector<Fixture<false>, Fixture<true>> Fixtures;

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex_transitions, ModelFixture)

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_weak_links, F, Fixtures) {
        F fixture;
        fixture.vertex_model.prolog();
        test_custom_links(fixture.vertex_model);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_T1, F, Fixtures) {
        F fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        const auto& am = model.get_am();

        const auto& cells = am.cells();
        const auto& edges = am.edges();
        const auto& vertices = am.vertices();

        auto num_cells = cells.size();
        auto num_edges = edges.size();
        auto num_vertices = vertices.size();

        // use the first non boundary edge
        auto edge_it = std::find_if(edges.begin(), edges.end(),
                                  [am](const auto& edge) {
                                      return (not am.is_boundary(edge));
                                  });
        BOOST_TEST((edge_it != edges.end()));
        auto edge = *edge_it;
        BOOST_TEST(not am.is_boundary(edge));

        // remember the adjoint cells
        const auto [adj_a, adj_b] = am.adjoints_of(edge);
        const std::size_t adj_a_num_es = adj_a->custom_links().edges.size();
        const std::size_t adj_b_num_es = adj_b->custom_links().edges.size();

        auto adj_cs_a = am.adjoint_cells_of(edge->custom_links().a);
        BOOST_TEST(adj_cs_a.size() == 3);
        const std::size_t adj_c_num_es = std::accumulate(
            adj_cs_a.begin(), adj_cs_a.end(), 0.,
            [](double val, const auto& c) {
                return val + c->custom_links().edges.size();
            }) - adj_a_num_es - adj_b_num_es;
        auto adj_cs_b = am.adjoint_cells_of(edge->custom_links().b);
        BOOST_TEST(adj_cs_b.size() == 3);
        const std::size_t adj_d_num_es = std::accumulate(
            adj_cs_b.begin(), adj_cs_b.end(), 0.,
            [](double val, const auto& c) {
                return val + c->custom_links().edges.size();
            }) - adj_a_num_es - adj_b_num_es;

        // move both vertices to the center of the edge
        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;
        PCPVertex::SpaceVec center = .5*(am.position_of(a) + am.position_of(b));
        am.move_to(a, 0.999*center);
        am.move_to(b, center);

        // this should remove the edges since L_edge < threshold
        model.iterate();

        // check that last instance of edge in this function
        BOOST_TEST(edge.use_count() == 1);

        // check that the number of objects did not change
        BOOST_TEST(cells.size() == num_cells);
        BOOST_TEST(edges.size() == num_edges);
        BOOST_TEST(vertices.size() == num_vertices);

        // check the general links in the model
        test_custom_links(model);

        // check that the number of edges and vertices in the adjoint cells a
        // and b decreased by 1
        const std::size_t new_adj_a_num_es = adj_a->custom_links().edges.size();
        const std::size_t new_adj_b_num_es = adj_b->custom_links().edges.size();

        const std::size_t new_adj_c_num_es = std::accumulate(
            adj_cs_a.begin(), adj_cs_a.end(), 0.,
            [](double val, const auto& c) {
                if (not c) { throw; }
                return val + c->custom_links().edges.size();
            }) - new_adj_a_num_es - new_adj_b_num_es;

        const std::size_t new_adj_d_num_es = std::accumulate(
            adj_cs_b.begin(), adj_cs_b.end(), 0.,
            [](double val, const auto& c) {
                return val + c->custom_links().edges.size();
            }) - new_adj_a_num_es - new_adj_b_num_es;
        
        // check that the number of edges and vertices in the adjoint cells a
        // and b decreased by 1
        BOOST_TEST(adj_a_num_es - 1 == adj_a->custom_links().edges.size());
        BOOST_TEST(adj_b_num_es - 1 == adj_b->custom_links().edges.size());
        // NOTE num_edges == num_vertices (check in test_custom_links())

        // check that the number of edges and vertices in the adjoint cells c
        // and d increased by 1
        BOOST_TEST(adj_c_num_es + 1 == new_adj_c_num_es);
        BOOST_TEST(adj_d_num_es + 1 == new_adj_d_num_es);
        // NOTE num_edges == num_vertices (check in test_custom_links())
    }

    BOOST_AUTO_TEST_CASE(test_T1_2_cell_boundary_edge) {
        Fixture<false> fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        const auto& am = model.get_am();

        const auto& cells = am.cells();
        const auto& edges = am.edges();
        const auto& vertices = am.vertices();

        auto num_cells = cells.size();
        auto num_edges = edges.size();
        auto num_vertices = vertices.size();

        // use the first non boundary edge
        auto edge_it = std::find_if(edges.begin(), edges.end(),
                                  [am](const auto& edge) {
                                      return am.is_2_cell_boundary_edge(edge);
                                  });
        BOOST_TEST((edge_it != edges.end()));
        auto edge = *edge_it;
        BOOST_TEST(am.is_2_cell_boundary_edge(edge));
        BOOST_TEST(am.is_boundary(edge));

        // remember the adjoint cells
        const auto [adj_a, adj_b] = am.adjoints_of(edge);
        const std::size_t adj_a_num_es = adj_a->custom_links().edges.size();
        const std::size_t adj_b_num_es = adj_b->custom_links().edges.size();

        // restrict test
        auto vertex_a = edge->custom_links().a;
        auto vertex_b = edge->custom_links().b;
        if (not am.is_boundary(vertex_b)) {
            std::swap(vertex_a, vertex_b);
        }
        BOOST_TEST(not am.is_boundary(vertex_a));
        BOOST_TEST(am.is_boundary(vertex_b));
        BOOST_TEST(am.is_3_fold_boundary_vertex(vertex_b));

        auto adj_cs_a = am.adjoint_cells_of(vertex_a);
        BOOST_TEST(adj_cs_a.size() == 3);
        const std::size_t adj_c_num_es = std::accumulate(
            adj_cs_a.begin(), adj_cs_a.end(), 0.,
            [](double val, const auto& c) {
                return val + c->custom_links().edges.size();
            }) - adj_a_num_es - adj_b_num_es;
        auto adj_cs_b = am.adjoint_cells_of(vertex_b);
        BOOST_TEST(adj_cs_b.size() == 2);

        // move both vertices to the center of the edge
        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;
        PCPVertex::SpaceVec center = .5*(am.position_of(a) + am.position_of(b));
        am.move_to(a, 0.999*center);
        am.move_to(b, center);

        // information to identify new objects
        const std::size_t max_edge_id = (*std::max_element(edges.begin(),
                                            edges.end(),
                                            [](const std::shared_ptr<Edge>& a,
                                               const std::shared_ptr<Edge>& b)
                                            {
                                                return (a->id() < b->id());
                                            }))->id();
        const std::size_t max_vertex_id = (*std::max_element(vertices.begin(),
                                            vertices.end(),
                                            [](const std::shared_ptr<Vertex>& a,
                                               const std::shared_ptr<Vertex>& b) {
                                                return (a->id() < b->id());
                                            }))->id();

        // this should remove the edges since L_edge < threshold
        model.iterate();

        // check that last instance of edge in this function
        BOOST_TEST(edge.use_count() == 1);

        // check that the number of objects did not change
        BOOST_TEST(cells.size() == num_cells);
        BOOST_TEST(edges.size() == num_edges);
        BOOST_TEST(vertices.size() == num_vertices);

        // check the general links in the model
        test_custom_links(model);

        // check that the number of edges and vertices in the adjoint cells a
        // and b decreased by 1
        const std::size_t new_adj_a_num_es = adj_a->custom_links().edges.size();
        const std::size_t new_adj_b_num_es = adj_b->custom_links().edges.size();

        const std::size_t new_adj_c_num_es = std::accumulate(
            adj_cs_a.begin(), adj_cs_a.end(), 0.,
            [](double val, const auto& c) {
                if (not c) { throw; }
                return val + c->custom_links().edges.size();
            }) - new_adj_a_num_es - new_adj_b_num_es;
        // NOTE adj_cs_b is only cells a and b
        
        // check that the number of edges and vertices in the adjoint cells a
        // and d decreased by 1
        BOOST_TEST(adj_a_num_es - 1 == adj_a->custom_links().edges.size());
        BOOST_TEST(adj_b_num_es - 1 == adj_b->custom_links().edges.size());
        // NOTE num_edges == num_vertices (check in test_custom_links())

        // check that the number of edges and vertices in the adjoint cells c
        // and c increased by 1
        BOOST_TEST(adj_c_num_es + 1 == new_adj_c_num_es);
        // NOTE num_edges == num_vertices (check in test_custom_links())

        // test the new objects
        for (const auto& c : adj_cs_a) {
            BOOST_TEST(am.is_boundary(c));
        }
        AgentContainer<Edge> new_edges;
        std::copy_if(edges.begin(), edges.end(), std::back_inserter(new_edges),
                     [am, max_edge_id](const auto& e) {
                         return e->id() > max_edge_id;
                     });
        for (const auto& e : new_edges) {
            BOOST_TEST(am.is_boundary(e));
        }

        AgentContainer<Vertex> new_vertices;
        std::copy_if(vertices.begin(), vertices.end(),
                     std::back_inserter(new_vertices),
                     [am, max_vertex_id](const auto& v) {
                         return v->id() > max_vertex_id;
                     });
        for (const auto& v : new_vertices) {
            BOOST_TEST(am.is_boundary(v));
        }
    }

    BOOST_AUTO_TEST_CASE(test_T1_remove_1_cell_boundary_edge_two_fold_vertices)
    {
        Fixture<false> fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        const auto& am = model.get_am();

        const auto& cells = am.cells();
        const auto& edges = am.edges();
        const auto& vertices = am.vertices();

        auto num_cells = cells.size();
        auto num_edges = edges.size();
        auto num_vertices = vertices.size();

        // use the first non boundary edge
        auto edge_it = std::find_if(edges.begin(), edges.end(),
            [am](const auto& edge) {
                const auto& a = edge->custom_links().a;
                const auto& b = edge->custom_links().b;
                return (    am.is_1_cell_boundary_edge(edge)
                        and (    am.is_2_fold_boundary_vertex(a)
                             and am.is_2_fold_boundary_vertex(b)));
            });
        BOOST_TEST((edge_it != edges.end()));
        auto edge = *edge_it;
        BOOST_TEST(am.is_1_cell_boundary_edge(edge));
        BOOST_TEST(am.is_boundary(edge));

        // remember the adjoint cells
        auto [adj_a, adj_b] = am.adjoints_of(edge);
        if (not adj_a) {
            std::swap(adj_a, adj_b);
        }
        BOOST_TEST(adj_a);
        BOOST_TEST(not adj_b);
        const std::size_t adj_a_num_es = adj_a->custom_links().edges.size();

        // restrict test
        auto vertex_a = edge->custom_links().a;
        auto vertex_b = edge->custom_links().b;
        BOOST_TEST(am.is_boundary(vertex_a));
        BOOST_TEST(am.is_2_fold_boundary_vertex(vertex_a));
        BOOST_TEST(am.is_boundary(vertex_b));
        BOOST_TEST(am.is_2_fold_boundary_vertex(vertex_b));

        auto adj_cs_a = am.adjoint_cells_of(vertex_a);
        BOOST_TEST(adj_cs_a.size() == 1);
        auto adj_cs_b = am.adjoint_cells_of(vertex_b);
        BOOST_TEST(adj_cs_b.size() == 1);

        // move both vertices to the center of the edge
        auto a = edge->custom_links().a;
        auto b = edge->custom_links().b;
        PCPVertex::SpaceVec center = .5*(am.position_of(a) + am.position_of(b));
        am.move_to(a, 0.999*center);
        am.move_to(b, center);

        // information to identify new objects
        const std::size_t max_vertex_id = (*std::max_element(vertices.begin(),
                                            vertices.end(),
                                            [](const std::shared_ptr<Vertex>& a,
                                               const std::shared_ptr<Vertex>& b) {
                                                return (a->id() < b->id());
                                            }))->id();

        // this should remove the edges since L_edge < threshold
        model.iterate();

        // check that last instance of edge in this function
        BOOST_TEST(edge.use_count() == 1);

        // check that the number of objects did not change
        BOOST_TEST(cells.size() == num_cells);
        BOOST_TEST(edges.size() == num_edges - 1);
        BOOST_TEST(vertices.size() == num_vertices - 1);

        // check the general links in the model
        test_custom_links(model);
        
        // check that the number of edges and vertices in the adjoint cells a
        // decreased by 1
        BOOST_TEST(adj_a_num_es - 1 == adj_a->custom_links().edges.size());
        // NOTE num_edges == num_vertices (check in test_custom_links())

        // test the new objects
        for (const auto& c : adj_cs_a) {
            BOOST_TEST(am.is_boundary(c));
        }

        AgentContainer<Vertex> new_vertices;
        std::copy_if(vertices.begin(), vertices.end(),
                     std::back_inserter(new_vertices),
                     [am, max_vertex_id](const auto& v) {
                         return v->id() > max_vertex_id;
                     });

        BOOST_TEST(new_vertices.size() == 1);
        for (const auto& v : new_vertices) {
            BOOST_TEST(am.is_boundary(v));
            BOOST_TEST(am.is_2_fold_boundary_vertex(v));
        }
    } // end test T1 1 cell boundary edge 2 x 2 fold vertices

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_T2, F, Fixtures) {
        F fixture;
        auto& model = fixture.vertex_model;
        model.prolog();
        
        auto& am = model.get_am();

        auto& cells = am.cells();
        auto& edges = am.edges();
        auto& vertices = am.vertices();

        auto num_cells = cells.size();
        auto num_edges = edges.size();
        auto num_vertices = vertices.size();

        // find the first cell that is not a boundary cell
        auto cell = *std::find_if(cells.begin(), cells.end(),
                                  [am](const auto& cell) {
                                      return (not am.is_boundary(cell));
                                  });
        BOOST_TEST(not am.is_boundary(cell));

        cell->state.type = 1;

        // remove edges until cell is triangular
        std::size_t cnt = 0;
        while (cell->custom_links().edges.size() > 3) {
            const auto& edges = cell->custom_links().edges;
            const auto& [edge, flip] = edges[123456789 % edges.size()];

            SpaceVec displ = am.displacement(edge);
            am.move_by(edge->custom_links().a,  0.49 * displ);
            am.move_by(edge->custom_links().b, -0.49 * displ);
    
            model.iterate();

            if (cnt > 10) {
                throw std::runtime_error("Unable to remove junctions to reach "
                    "triangular cell.");
            }
        }
        
        // move all vertices close to center
        SpaceVec center = am.barycenter_of(cell);
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            std::shared_ptr<Vertex> vertex;
            if (not flip) {
                vertex = edge->custom_links().a;
            }
            else {
                vertex = edge->custom_links().b;
            }

            SpaceVec pos = am.position_of(vertex);
            SpaceVec displ = model.get_space()->displacement(pos, center);

            am.move_by(vertex, 0.99 * displ);
        }

        model.iterate();

        BOOST_TEST((   std::find(cells.begin(), cells.end(), cell)
                    == cells.end()));

        BOOST_TEST(cell.use_count() == 1);

        BOOST_TEST(cells.size() == num_cells - 1);
        BOOST_TEST(edges.size() == num_edges - 3);
        BOOST_TEST(vertices.size() == num_vertices - 2);

        test_custom_links(model);
    } // test T2 cell removal

    BOOST_AUTO_TEST_CASE(test_T2_boundary) {
        Fixture<false> fixture;
        auto& model = fixture.vertex_model;
        model.prolog();
        
        auto& am = model.get_am();

        auto& cells = am.cells();
        auto& edges = am.edges();
        auto& vertices = am.vertices();

        auto num_cells = cells.size();
        auto num_edges = edges.size();
        auto num_vertices = vertices.size();

        // find the first cell that is a boundary cell
        auto cell = *std::find_if(cells.begin(), cells.end(),
                                [am](const auto& cell) {
                                    return am.is_boundary(cell);
                                });
        BOOST_TEST(am.is_boundary(cell));
        
        // move all vertices close to center
        SpaceVec center = am.barycenter_of(cell);
        for (const auto& [edge, flip] : cell->custom_links().edges) {
            std::shared_ptr<Vertex> vertex;
            if (not flip) {
                vertex = edge->custom_links().a;
            }
            else {
                vertex = edge->custom_links().b;
            }

            SpaceVec pos = am.position_of(vertex);
            SpaceVec displ = model.get_space()->displacement(pos, center);

            am.move_by(vertex, 0.99 * displ);
        }

        model.iterate();

        BOOST_TEST((   std::find(cells.begin(), cells.end(), cell)
                    == cells.end()));

        BOOST_TEST(cell.use_count() == 1);

        BOOST_TEST(cells.size() == num_cells - 1);
        BOOST_TEST(edges.size() <= num_edges - 3);
        BOOST_TEST(vertices.size() <= num_vertices - 2);
        // NOTE T1 boundary edge removal possible

        test_custom_links(model);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_divide_cell, F, Fixtures) {
        F fixture;
        auto& model = fixture.vertex_model;
        bool periodic_BC = model.get_space()->periodic;

        model.prolog();

        const auto& am = model.get_am();

        const auto cells = am.cells();
        const auto edges = am.edges();
        const auto vertices = am.vertices();
        

        // find the first cell that is not a boundary cell
        auto cell = *std::find_if(cells.begin(), cells.end(),
                                  [am, periodic_BC](const auto& cell) {
                                      if (periodic_BC) {
                                          return (not am.is_boundary(cell));
                                      }
                                      else {
                                          return am.is_boundary(cell);
                                      }
                                  });
        if (periodic_BC) {
            BOOST_TEST(not am.is_boundary(cell));
            // NOTE this is a bulk cell
            // NOTE implicit test of removal of bulk cell in non-periodic BC
        }
        else {
            BOOST_TEST(am.is_boundary(cell));
        }

        // set arbitrary values and check inheritance
        cell->state.type = 1;
        cell->state.register_parameter("some_cell_value", {1.});
        model.divide_cell(cell, 0.);


        const auto cells_new = am.cells();
        const auto edges_new = am.edges();
        const auto vertices_new = am.vertices();

        BOOST_TEST(cells.size() + 1 == cells_new.size());
        BOOST_TEST(edges.size() + 3 == edges_new.size());
        BOOST_TEST(vertices.size() + 2 == vertices_new.size());

        for (auto new_cell : {cells_new[cells_new.size() - 2],
                              cells_new[cells_new.size() - 1]})
        {
            BOOST_TEST(new_cell->state.list_parameters().size()
                       ==  cell->state.list_parameters().size());
            BOOST_TEST(new_cell->state.type
                       ==  cell->state.type);
            BOOST_TEST(new_cell->state.get_parameter("some_cell_value")[0]
                       ==  cell->state.get_parameter("some_cell_value")[0]);

            new_cell->state.type = 0;
            new_cell->state.unregister_parameter("some_cell_value");
        }

        test_custom_links(model);

        BOOST_TEST_MESSAGE("Iterating model after cell division.");
        
        // check that it does not fail somewhere ...
        model.iterate();
        test_custom_links(model);
    } // test divide cell

    
    BOOST_AUTO_TEST_CASE(test_triangle_deformation) {
        using SpaceVec = PCPVertex::SpaceVec;

        Fixture<true> fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        const auto& am = model.get_am();
        

        // the test
        std::function<bool(SpaceVec, SpaceVec, SpaceVec, double, double, double,
                           std::string, SpaceVec, SpaceVec, SpaceVec)>
        _test_triangle = [am](SpaceVec a, SpaceVec b, SpaceVec c,
                              double area, double ratio, double theta,
                              std::string message,
                              SpaceVec _a, SpaceVec _b, SpaceVec _c)
        {
            am.get_logger()->debug(message);     

            /// Transform vertices to shape matrix
            auto T = std::make_tuple(a, b, c);
            arma::mat22 shape = am.shape<false, false>(T);

            // test that shape is the transformation matrix
            arma::mat22 base = arma::mat22({{(_b-_a)[0], (_c-_b)[0]},
                                            {(_b-_a)[1], (_c-_b)[1]}});
            bool test_transformation = arma::approx_equal(
                base,
                shape * am.equilateral_triangle(),
                "absdiff", 1.e-4
            );
            if (not test_transformation)
            {   
                std::cout << "Failed to get transformation matrix: \n"
                          << base << std::endl
                          << shape * am.equilateral_triangle() << std::endl;
            }
            BOOST_TEST(test_transformation);

            // test how 
            auto [tr, sym, asym] = am.decompose(shape);
            arma::mat22 decomposition = sym + asym + arma::eye(2, 2) * tr / 2.;
            auto test_decomposition = arma::approx_equal(
                shape, decomposition, "absdiff", 1.e-10
            );
            if (not test_decomposition)
            {   
                std::cout << "Failed to decompose matrix: \n"
                          << shape << std::endl
                          << "Decomposition is \n"
                          << decomposition << std::endl
                          << "with trace " << tr << std::endl
                          << "symmetric \n"
                          << sym << std::endl
                          << "and anti-symmetric \n"
                          << asym << std::endl
                          << "components.\n";
            }
            BOOST_TEST(test_decomposition);

            auto [A, q, r, Theta] = am.interpretation(shape);
            arma::mat22 interpretation = (
                sqrt(A) * arma::expmat(q)
                * arma::mat22({{cos(Theta), -sin(Theta)},
                               {sin(Theta),  cos(Theta)}}));
            bool test_interpretation = arma::approx_equal(
                    shape, interpretation,
                    "absdiff", 1.e-10
            );
            if (not test_interpretation) {
                std::cout << "Shape interpretation failed! \n Shape is \n"
                          << shape << std::endl
                          << " and interpretation is \n"
                          << interpretation << std::endl
                          << "with area " << A
                          << " axis ratio " << r
                          << " and rotation angle " << Theta;

            }
            BOOST_TEST(test_interpretation);
            BOOST_CHECK_CLOSE(A, area, 1.e-6);
            BOOST_CHECK_CLOSE(r, ratio, 1.e-6);
            BOOST_CHECK_SMALL(
                (  std::fmod(theta - Theta + 7 * M_PI / 3., 2. * M_PI / 3.)
                 - M_PI / 3.),
                1.e-8);
            // NOTE shift both to avoid comparison 0 == 2 * M_PI

            return true;
        };

        std::function<bool(SpaceVec, SpaceVec, SpaceVec, double, double, double,
                           std::string)>
        test_triangle = [am, _test_triangle](SpaceVec a, SpaceVec b, SpaceVec c,
                             double area, double ratio, double theta,
                             std::string message = "")
        {
            return _test_triangle(a, b, c, area, ratio, theta, message,
                                  a, b, c);
        };

        // the equilateral triangle
        SpaceVec a({0., 0.}); 
        SpaceVec b({1., 0.});
        SpaceVec c({0.5, sqrt(3) / 2.});
        

        // the expected outcome
        double _A = 1.;
        double _r = 1.;
        double _theta = 0.;

        BOOST_TEST(test_triangle(a, b, c, _A, _r, _theta, "equilateral"));
        BOOST_TEST(_test_triangle(b, c, a, _A, _r, _theta, "equilateral swap",
                                  a, b, c));
        BOOST_TEST(_test_triangle(c, a, b, _A, _r, _theta, "equilateral swap 2",
                                  a, b, c));
        BOOST_TEST(_test_triangle(a, c, b, _A, _r, _theta, "equilateral swap 3",
                                  a, b, c));


        // shift that barycenter is at origin
        SpaceVec shift = (a + b + c) / 3.;
        a -= shift;
        b -= shift;
        c -= shift;

        BOOST_TEST(test_triangle(a, b, c, _A, _r, _theta, "shift"));


        // rotate triangle
        std::size_t N_rot = 13;
        for (std::size_t i = 0; i < N_rot; i++) {
            double rotate = 2 * M_PI / N_rot;
            arma::mat22 rot({{cos(rotate), -sin(rotate)},
                             {sin(rotate),  cos(rotate)}});

            a = rot * a;    
            b = rot * b;
            c = rot * c;

            _theta += rotate;
            
            // ensure c uppermost vertex
            if (_theta < M_PI / 3. or _theta > 2 * M_PI - M_PI / 3.) {
                BOOST_TEST(test_triangle(a, b, c, _A, _r, _theta,
                                         fmt::format("rotate {}", _theta)));
            }
            else if (_theta < 3 * M_PI / 3.) {
                BOOST_TEST(_test_triangle(a, b, c, _A, _r, _theta,
                                          fmt::format("rotate {}", _theta),
                                          c, a, b));
            }
            else {
                BOOST_TEST(_test_triangle(a, b, c, _A, _r, _theta,
                                          fmt::format("rotate {}", _theta),
                                          b, c, a));
            }
        }


        // scale by 2
        a = 2 * a;
        b = 2 * b;
        c = 2 * c;

        _A *= 4;

        BOOST_TEST(test_triangle(a, b, c, _A, _r, _theta, "scale"));


        // reset 
        a /= 2.;
        b /= 2.;
        c /= 2.;
        _A = 1.;
        _r = 1.;
        _theta = 0.;
        BOOST_TEST(test_triangle(a, b, c, _A, _r, _theta, "reset"));


        // squeeze
        a = SpaceVec({a[0], a[1] / 2.});
        b = SpaceVec({b[0], b[1] / 2.});
        c = SpaceVec({c[0], c[1] / 2.});
        _A /= 2.;
        _r *= 2.;
        _theta = 0.;


        // rotate
        for (std::size_t i = 0; i < N_rot; i++) {
            double rotate = 2 * M_PI / N_rot;
            arma::mat22 rot({{cos(rotate), -sin(rotate)},
                             {sin(rotate),  cos(rotate)}});

            a = rot * a;    
            b = rot * b;
            c = rot * c;

            _theta += rotate;
            
            // ensure c uppermost vertex
            if (c[1] > a[1] and c[1] > b[1]) {
                BOOST_TEST(test_triangle(a, b, c, _A, _r, _theta,
                                         fmt::format("rotate {}", _theta)));
            }
            else if (b[1] > a[1] and b[1] > c[1]) {
                BOOST_TEST(_test_triangle(a, b, c, _A, _r, _theta,
                                          fmt::format("rotate {}", _theta),
                                          c, a, b));
            }
            else {
                BOOST_TEST(_test_triangle(a, b, c, _A, _r, _theta,
                                          fmt::format("rotate {}", _theta),
                                          b, c, a));
            }
        }
    }

    BOOST_AUTO_TEST_CASE(test_cell_elongation) {
        using SpaceVec = PCPVertex::SpaceVec;

        Fixture<true> fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        const auto& am = model.get_am();
        SpaceVec domain0 = model.get_space()->get_domain_size();

        std::function<void(arma::mat22, double, double, double)>
        check_elongation = [](arma::mat22 q, double qxx, double qxy,
                              double p)
        {
            if (fabs(qxx) > 0.1 * p) {
                BOOST_CHECK_CLOSE(q[0],  qxx, p * 100);
                BOOST_CHECK_CLOSE(q[3], -qxx, p * 100);
            }
            else {
                BOOST_CHECK_SMALL(q[0], p);
                BOOST_CHECK_SMALL(q[3], p);
            }
            if (fabs(qxy) > 0.1 * p) {
                BOOST_CHECK_CLOSE(q[1], qxy, p * 100);
                BOOST_CHECK_CLOSE(q[2], qxy, p * 100);
            }
            else {
                BOOST_CHECK_SMALL(q[1], p);
                BOOST_CHECK_SMALL(q[2], p);
            }
        };


        // for a cell
        auto cell = am.cells()[5];
        // NOTE cell is from initialisation a hexagon

        am.get_logger()->debug("Testing shape of cell {} ..", cell->id());
        check_elongation(am.elongation_of(cell), 0., 0., 1.e-5);
        std::map<std::size_t, SpaceVec> initial_positions;
        for (const auto& v : am.vertices()) {
            initial_positions[v->id()] = am.position_of(v);
        }
        std::map<std::size_t, SpaceVec> initial_edges;
        for (const auto& e : am.edges()) {
            initial_edges[e->id()] = am.displacement(e);
        }


        // pure shear domain
        double rho = 0.;        
        for (std::size_t i = 0; i < 5; i++) {
            rho += 0.2;
            SpaceVec _stretch = domain0 % SpaceVec({exp(rho), exp(-rho)});
            _stretch -= model.get_space()->get_domain_size();
            model.stretch_domain(_stretch, false, false, false, true);
            // NOTE deform plastic
            
            am.get_logger()->debug("Pure shear reached {}", rho);
            check_elongation(am.elongation_of(cell), rho, 0., 1.e-3);
        }
        model.stretch_domain(domain0 - model.get_space()->get_domain_size(),
                             false, false, false, true);

        am.get_logger()->debug("Reset ..");
        check_elongation(am.elongation_of(cell), 0., 0., 1.e-5);


        // simple shear domain
        std::vector<SpaceVec> simple_shear_experiments({
               SpaceVec({1., 0.0}),  // simple shear
               SpaceVec({0.0, 1.})   // simple shear
        });
        for (SpaceVec add_skew : simple_shear_experiments) {
            model.skew_domain(add_skew, true, true);
            SpaceVec skew = model.get_space()->get_skew();
            for (const auto& c : am.cells()) {
                BOOST_CHECK_CLOSE(am.area_of(c), 1., 1.e-2);
            }

            am.get_logger()->debug("Simple shear is ({}, {})",
                                skew[0], skew[1]);
            arma::mat22 q = am.elongation_of<false>(cell);
            double strain = (skew[0]/domain0[1] + skew[1]/domain0[0]) / 2.;
            double abs_q = sqrt(arma::trace(q * q.t()) / 2.);
            BOOST_CHECK_CLOSE(abs_q, strain, 0.5); // precise to 5 permille


            am.get_logger()->debug("Reset ..");
            model.skew_domain(-1. * model.get_space()->get_skew(), true, true);
            // NOTE undo skew does not work if skew_x != 0  and skew_y != 0
            //  ->  use brute force
            for (const auto& v : am.vertices()) {
                am.move_to(v, initial_positions[v->id()]);
            }
            BOOST_CHECK_SMALL(arma::norm(model.get_space()->get_skew()), 1.e-5);
            
            check_elongation(am.elongation_of(cell), 0., 0., 1.e-5);
            double max_displ = 0;
            for (const auto& v : am.vertices()) {
                SpaceVec _pos = initial_positions[v->id()];
                SpaceVec pos = am.position_of(v);
                SpaceVec displ = model.get_space()->displacement(pos, _pos);
                max_displ = std::max(max_displ, arma::norm(displ));
            }
            BOOST_CHECK_SMALL(max_displ, 1.e-5);
        }


        // isotropic expansion
        SpaceVec stretch(arma::fill::zeros);
        for (std::size_t i = 0; i < 4; i++) {
            SpaceVec _stretch = 0.01 * model.get_space()->get_domain_size();
            model.stretch_domain(_stretch, false, false, false);
            stretch += _stretch;

            am.get_logger()->debug("Isotropic expansion ..");
            check_elongation(am.elongation_of(cell), 0., 0., 1.e-5);
        }

        model.stretch_domain(domain0 - model.get_space()->get_domain_size(),
                             false, false, false);
        check_elongation(am.elongation_of(cell), 0., 0., 1.e-5);
    }

BOOST_AUTO_TEST_SUITE_END()
