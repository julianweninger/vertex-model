#define BOOST_TEST_MODULE PCPVertexTestTransitions

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <numeric>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "utils.hh"
#include "../PCPVertex.hh"
#include "../energy.hh"
#include "../algorithm.hh"
#include "../operations.hh"
#include "../PCPVertex_write_tasks.hh"

using namespace Utopia;
using namespace Utopia::Models::PCPVertex;
using Utopia::DataIO::Config;

using Vertex = Utopia::Models::PCPVertex::PCPVertex::Vertex;
using Edge = Utopia::Models::PCPVertex::PCPVertex::Edge;
using Cell = Utopia::Models::PCPVertex::PCPVertex::Cell;

template<bool periodic>
struct Fixture {    
    Models::PCPVertex::PCPVertex vertex_model;

    Fixture ()
    :
        vertex_model(model_factory())
    { }

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
        cell->state._area_preferential = 0.314;
        cell->state.type = PCPVertex::CellType::support;
        cell->state.shape_index_preferential = 4.;
        cell->state.contractility = 0.114;
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
            BOOST_TEST(new_cell->state.area_preferential()
                       ==  cell->state.area_preferential());
            BOOST_TEST(new_cell->state.type
                       ==  cell->state.type);
            BOOST_TEST(new_cell->state.shape_index_preferential
                       ==  cell->state.shape_index_preferential);
            BOOST_TEST(new_cell->state.contractility
                       ==  cell->state.contractility);
        }

        test_custom_links(model);

        BOOST_TEST_MESSAGE("Iterating model after cell division.");
        
        // check that it does not fail somewhere ...
        model.iterate();
        test_custom_links(model);
    } // test divide cell

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

        edge->state._linetension = 200.;

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

        edge->state._linetension = 200.;

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

        edge->state._linetension = 200.;

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

        // decrement area in small steps to avoid numerical errors!
        for (int i = 0; i < 9; i++) {
            if (not cell) { break; }

            cell->state._area_preferential -= 0.1;
            
            for (int i = 0; i < 250; i++) {
                model.iterate();
            }
        }

        if (cell) {
            cell->state._area_preferential -= 0.075;
            
            for (int i = 0; i < 10; i++) {
                model.iterate();
            }
        }

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

        // decrement area in small steps to avoid numerical errors!
        for (int i = 0; i < 9; i++) {
            if (not cell) { break; }

            cell->state._area_preferential -= 0.1;
            
            for (int i = 0; i < 250; i++) {
                model.iterate();
            }
        }

        if (cell) {
            cell->state._area_preferential -= 0.075;
            
            for (int i = 0; i < 10; i++) {
                model.iterate();
            }
        }

        BOOST_TEST((   std::find(cells.begin(), cells.end(), cell)
                    == cells.end()));

        BOOST_TEST(cell.use_count() == 1);

        BOOST_TEST(cells.size() == num_cells - 1);
        BOOST_TEST(edges.size() <= num_edges - 3);
        BOOST_TEST(vertices.size() <= num_vertices - 2);
        // NOTE T1 boundary edge removal possible

        test_custom_links(model);
    }

    
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

        // for a cell
        auto cell = am.cells()[0];
        // NOTE cell is from initialisation a hexagon

        double _theta = 0.;
        double _r = 1.;

        am.get_logger()->debug("Testing shape of cell {} ..", cell->id());

        for (const auto& [edge, flip] : cell->custom_links().edges) {
            auto a = edge->custom_links().a;
            auto b = edge->custom_links().b;
            if (flip) { std::swap(a, b); }
        }

        std::function<bool(
            const std::shared_ptr<typename PCPVertex::Cell>&,
            double, double,
            std::string)>
        test_elongation_of =
        [am](const auto& cell, double _r, double _theta,
             std::string message = "")
        {
            am.get_logger()->debug(message);

            arma::mat22 elongation = am.elongation_of(cell);
            SpaceVec tmp = elongation.col(0);
            double theta;
            if (fabs(tmp[0]) < 1.e-12) {
                if (tmp[1] > 1.e-12) {
                    theta = M_PI_2;
                }
                else {
                    theta = 3 * M_PI_2;
                }
            }
            else if (tmp[0] > 1.e-12) {
                theta = std::fmod(atan(tmp[1] / tmp[0]) + 2 * M_PI, 2 * M_PI);
            }
            else {
                if (tmp[1] > 0.) {
                    theta = M_PI - atan(- tmp[1] / tmp[0]);
                }
                else  {
                    theta = M_PI + atan(tmp[1] / tmp[0]);
                }
            }

            double r = exp(2*sqrt(arma::trace(elongation * elongation.t())/2.));

            BOOST_CHECK_CLOSE(r, _r, 1);
            // TODO error of 1 % might be to big for this problem !!!!
            
            if (fabs(1. - r) > 1.e-4) {
                theta = std::fmod(theta + 3 * M_PI_2, M_PI) - M_PI_2;
                _theta = std::fmod(theta + 3 * M_PI_2, M_PI) - M_PI_2;
                if (  std::fmod(theta - _theta + 3 * M_PI_2, M_PI) - M_PI_2
                    > 1.e-3)
                {
                    am.get_logger()->error("angle of cell {} does not "
                                           "correspond to expected angle {}",
                                           theta, _theta);
                }
                BOOST_CHECK_SMALL(
                    std::fmod(theta - _theta + 3 * M_PI_2, M_PI) - M_PI_2,
                    1.e-3);
            }

            return true;
        };

        // pure shear domain
        SpaceVec stretch({0., 0.});
        SpaceVec domain0 = model.get_space()->get_domain_size();
        SpaceVec domain = model.get_space()->get_domain_size();
        test_elongation_of(cell, _r, _theta, fmt::format("stretch ({}, {})",
                                                     stretch[0], stretch[1]));
        for (std::size_t i = 0; i < 4; i++) {
            SpaceVec _stretch({0.2, -0.2});
            model.stretch_domain(_stretch, false, false, false);
            stretch += _stretch;
            domain = model.get_space()->get_domain_size();
            _r = domain[0] / domain[1] / (domain0[0] / domain0[1]);
            test_elongation_of(cell, _r, _theta, fmt::format("stretch ({}, {})",
                                                         stretch[0],
                                                         stretch[1]));
        }


        // isotropic expansion
        for (std::size_t i = 0; i < 4; i++) {
            SpaceVec _stretch = 0.01 * model.get_space()->get_domain_size();
            model.stretch_domain(_stretch, false, false, false);
            stretch += _stretch;
            domain = model.get_space()->get_domain_size();
            test_elongation_of(cell, _r, _theta, fmt::format("expand ({}, {})",
                                                             stretch[0],
                                                             stretch[1]));
        }

        model.get_space()->set_domain_size(domain0);
        for (const auto& vertex : am.vertices()) {
            SpaceVec pos = am.position_of(vertex);
            am.move_to(vertex, am.position_of(vertex) / domain % domain0);
        }
        _r = 1.;
        _theta = 0.;
        test_elongation_of(cell, _r, _theta, "reset.");


        // simple shear domain
        double skew = 0.;
        domain = model.get_space()->get_domain_size();
        for (std::size_t i = 0; i < 5; i++) {
            BOOST_CHECK_CLOSE(am.area_of(cell), 1., 1.e-2);
            double _skew = 0.25;
            skew += _skew;
            for (const auto& vertex : am.vertices()) {
                SpaceVec pos = am.position_of(vertex);
                SpaceVec displ({_skew * pos[1] / domain[1], 0.});
                am.move_by(vertex, displ);
            }
            model.get_space()->set_skew(SpaceVec({skew, 0.}));
            BOOST_CHECK_CLOSE(am.area_of(cell), 1., 1.e-2);

            _theta = asin(skew);
            _r = sqrt(  (std::pow(domain[0] + skew, 2) + std::pow(domain[1], 2))
                      / (std::pow(domain[0] - skew, 2) + std::pow(domain[1], 2))
                     );

            test_elongation_of(cell, _r, _theta, fmt::format("skew {}", skew));
        }


        // check on entire domain
        double dual_area = 0.;
        arma::mat22 av_q = arma::zeros(2, 2);

        arma::mat22 q_ref = am.elongation_of(cell);
        for (const auto& cell : am.cells()) {
            arma::mat22 q = am.elongation_of(cell);            
            double area = arma::det(q);

            dual_area += area;
            BOOST_TEST(arma::approx_equal(q_ref, q, "absdiff", 1.e-4));

            av_q += area * q;
        }
        av_q /= dual_area;

        
        SpaceVec tmp = av_q.col(0);
        double theta;
        if (fabs(tmp[0]) < 1.e-12) {
            if (tmp[1] > 1.e-12) {
                theta = M_PI_2;
            }
            else {
                theta = 3 * M_PI_2;
            }
        }
        else if (tmp[0] > 1.e-12) {
            theta = std::fmod(atan(tmp[1] / tmp[0]) + 2 * M_PI, 2 * M_PI);
        }
        else {
            if (tmp[1] > 0.) {
                theta = M_PI - atan(- tmp[1] / tmp[0]);
            }
            else  {
                theta = M_PI + atan(tmp[1] / tmp[0]);
            }
        }

        double r = exp(2*sqrt(arma::trace(av_q * av_q.t())/2.));

        BOOST_CHECK_CLOSE(r, _r, 1);
        // TODO error of 1 % might be to big for this problem !!!!
        
        if (fabs(1. - r) > 1.e-4) {
            theta = std::fmod(theta + 3 * M_PI_2, M_PI) - M_PI_2;
            _theta = std::fmod(theta + 3 * M_PI_2, M_PI) - M_PI_2;
            if (  std::fmod(theta - _theta + 3 * M_PI_2, M_PI) - M_PI_2
                > 1.e-3)
            {
                am.get_logger()->error("angle of cell {} does not "
                                        "correspond to expected angle {}",
                                        theta, _theta);
            }
            BOOST_CHECK_SMALL(
                std::fmod(theta - _theta + 3 * M_PI_2, M_PI) - M_PI_2, 
                1.e-5);
        }
    }

BOOST_AUTO_TEST_SUITE_END()
