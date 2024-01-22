#define BOOST_TEST_MODULE PCPVertexTest

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "utils.hh"

#include <utopia/core/apply.hh>
#include <utopia/core/testtools.hh>

#include "../PCPVertex.hh"
#include "../algorithm.hh"
#include "../operations.hh"
#include "../PCPVertex_write_tasks.hh"

#include "../work_function.hh"


using namespace Utopia;
using namespace Utopia::Models::PCPVertex;

using SpaceVec = typename PCPVertex::SpaceVec;

enum Case {
    periodic,
    non_periodic,
    columnar_periodic,
    columnar      // Columnar initialisation  
};

template<Case C>
struct Fixture {    
    const Case test_case;

    Models::PCPVertex::PCPVertex vertex_model;
    
    const Config wf_terms;

    Fixture ()
    :
        test_case(C),
        vertex_model(model_factory()),
        wf_terms(YAML::LoadFile("test_wf_terms.yml"))
    { }

    ~Fixture()
    {
        vertex_model.get_logger()->info("Tearing down ...");
        std::remove("test_data.h5");
        spdlog::drop_all();
    }
    
    PCPVertex model_factory() {
        using Utopia::Models::PCPVertex::DataIO::time_energy_adaptor;

        if constexpr (C == periodic) {
            Utopia::PseudoParent pp("test_periodic.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else if constexpr (C == non_periodic) {
            Utopia::PseudoParent pp("test.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else if constexpr (C == columnar) {
            Utopia::PseudoParent pp("test_columnar.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else if constexpr (C == columnar_periodic) {
            Utopia::PseudoParent pp("test_columnar_periodic.yml");
            return PCPVertex("PCPVertex", pp, {},
                            std::make_tuple(time_energy_adaptor));
        }
        else {
            throw std::runtime_error(fmt::format("Testcase {} not implemented",
                                                 C));
        }
    }
};

class TEST_PCPVertex_energy : public PCPVertex
{
public:
    template<class ParentModel>
    TEST_PCPVertex_energy (
        const std::string name,
        ParentModel &parent_model,
        const Utopia::DataIO::Config custom_cfg = {})
    :
        PCPVertex(name, parent_model, custom_cfg,
             std::make_tuple(Utopia::Models::PCPVertex::DataIO::time_energy_adaptor))
    { }

    /** Runs the full test */
    void run_test() {
        this->prolog();

        this->perform_test();

        this->epilog();
    }

    /** Performs the full test */
    void perform_test() {        
        // compare dE with energy
        this->jiggle_vertices(0.02);
        for (std::size_t i = 0; i < 100; i++) {
            test_energy_gradient();
        }
    }

    void test_energy_gradient() {
        const double precision = 1.e-2; // relative error in percent
        // NOTE only doing linear estimation here
        
        // chose a step size in steepest gradient update
        double dt = 0.00001;
        // NOTE error is expected to be smaller with small step size
        
        // set the gradient
        this->init_minimization();

        // the current energy
        auto energy = this->get_energy();

        const auto& am = this->get_am();

        // A [0,N]-range uniform distribution used for evaluating probabilities
        std::uniform_int_distribution<> prob_distr(0, am.vertices().size()-1);

        // move a random vertex
        const auto& random_vertex = am.vertices()[prob_distr(*this->_rng)];
        const SpaceVec pos = am.position_of(random_vertex);
        this->get_am().move_by(
            random_vertex,
            random_vertex->state.get_force() * dt
        );

        const SpaceVec d_pos = this->get_space()->displacement(
            pos, am.position_of(random_vertex));

        // integrate the gradient to obtain the energy change
        double dE = arma::dot(-1. * random_vertex->state.get_force(), d_pos);

        // compare to the change in energy
        BOOST_CHECK_CLOSE(dE, this->get_energy() - energy, precision);
    }
};

typedef boost::mpl::vector< Fixture<Case::periodic>,
                            Fixture<Case::non_periodic>,
                            Fixture<Case::columnar_periodic>,
                            Fixture<Case::columnar>
                          > Fixtures;

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex, ModelFixture)

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_model_initialization, F, Fixtures) {
        F fixture;
        test_custom_links(fixture.vertex_model);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_minimization, F, Fixtures) {
        F fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        Utopia::DataIO::Config minimization_cfg;

        std::cout << std::endl << "Beginning test of minimization in "
                                  "steepest gradient scheme .. \n\n";

        minimization_cfg["tolerance"] = 1e-6;
        minimization_cfg["dt"] = 1e-2;
        minimization_cfg["max_steps"] = 10000;
        minimization_cfg["num_repeat"] = 2;
        minimization_cfg["jiggle_tolerance"] = 1e-6;
        minimization_cfg["jiggle_intensity"] = 0.001;
        MinimizationParams minimization(minimization_cfg);

        model.minimize_energy(minimization);

        BOOST_TEST(model.get_time() >= 5);
        BOOST_TEST(model.get_energy_change() < minimization.tolerance);


        std::cout << std::endl << "Beginning test of minimization in "
                                  "steepest gradient noisy scheme .. \n\n";

        minimization_cfg["num_steps"] = 100;
        minimization_cfg["temperature"] = 0.0001;
        minimization_cfg["linetension_fluctuation"] = 0.0001;
        minimization_cfg["area_fluctuation"] = 0.0001;

        auto time_start = model.get_time();
        minimization = MinimizationParams(minimization_cfg);
        model.minimize_energy(minimization);

        BOOST_TEST(model.get_time() - time_start == 2 * (100 + 1) + 1);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE (test_WF_terms, F, Fixtures)
    {
        std::size_t test_term_index = 0;
        while (true) {
            F fixture;
            auto& model = fixture.vertex_model;

            Config work_function_terms = get_as<Config>(
                "test_work_function_terms", 
                fixture.wf_terms
            );
            
            if (test_term_index >= work_function_terms.size()) {
                BOOST_TEST_MESSAGE("Done.");
                break;
            }

            // this iterates below the level of defined work functions
            double cumulated_energy = 0.;
            bool test_skipped = true;
            for (const auto& term_pair : work_function_terms[test_term_index++])
            {
                auto term_name = term_pair.first.as<std::string>();
                auto params = term_pair.second;

                auto _test_cases = get_as<std::vector<std::string>>(
                    "skip_test_cases", params, {}
                );
                std::set<std::string> test_cases(_test_cases.begin(),
                                                 _test_cases.end());
                if (    fixture.test_case == Case::periodic
                    and test_cases.find("periodic") != test_cases.end())
                {
                    std::cout << "Skipping wf term " << term_name
                              << " on case " << "periodic";
                    continue;
                }
                if (    fixture.test_case == Case::non_periodic
                    and test_cases.find("non-periodic") != test_cases.end())
                {
                    std::cout << "Skipping wf term " << term_name
                              << " on case " << "non-periodic";
                    continue;
                }
                if (    fixture.test_case == Case::columnar_periodic
                    and test_cases.find("columnar-periodic") != test_cases.end())
                {
                    std::cout << "Skipping wf term " << term_name
                              << " on case " << "columnar-periodic";
                    continue;
                }
                if (    fixture.test_case == Case::columnar
                    and test_cases.find("columnar") != test_cases.end())
                {
                    std::cout << "Skipping wf term " << term_name
                              << " on case " << "columnar-non-periodic";
                    continue;
                }

                double precision = get_as<double>(
                    "test_precision", params, 1.e-3
                );

                const auto& am = model.get_am();
                const auto& cells = am.cells();

                // minimize the default work function terms
                model.prolog();


                // Differentiate 50 % type 1
                std::uniform_int_distribution<std::size_t> distr(0, cells.size()-1);
                std::size_t cnt_1 = 0;
                while(cnt_1 != cells.size() / 2) {
                    auto c_id = distr(*model.get_rng());
                    if (cells[c_id]->state.type == 0) {
                        cells[c_id]->state.type = 1;
                        cnt_1++;
                    }
                }
                BOOST_TEST(cnt_1 == cells.size() / 2);
                
                for (std::size_t i = 0; i < 50; i++) {
                    model.iterate();
                }
                
                // erase all work function terms
                while (model.get_work_function_terms().size() > 0) {
                    auto [name, term] = *model.get_work_function_terms().begin();
                    model.erase_work_function_term(name);
                }

                // register a single work function term
                auto name = get_as<std::string>("name", params, term_name);
                BOOST_TEST_MESSAGE("Work function term: " << name);
                model.register_work_function_term(term_name, name, params);

                const auto& term_function = model.get_work_function_term(name);
                cumulated_energy += fabs(term_function->compute_energy());
                test_skipped = false;

                term_function->test_constraints(model.get_logger());

                std::uniform_real_distribution<double> uniform_distr(0, 2*M_PI);
                std::size_t steps = 200;

                for (const auto& vertex : am.vertices()) {
                    const SpaceVec v_pos0 = am.position_of(vertex);

                    // chose a random direction of motion
                    double random_angle = uniform_distr(*model.get_rng());
                    SpaceVec displ = 0.01 * SpaceVec({cos(random_angle),
                                                      sin(random_angle)});

                    // calculate current force
                    vertex->state.reset_force();
                    term_function->compute_and_set_forces();
                    const double energy = term_function->compute_energy();

                    // integrate force along the displacement
                    double integral = 0.;
                    for (std::size_t i = 0; i < steps; i++) {
                        SpaceVec force = vertex->state.get_force();

                        SpaceVec tmp = am.position_of(vertex);
                        am.move_by(vertex, displ / double(steps));

                        SpaceVec delta = am.get_space()->displacement(
                            tmp,
                            am.position_of(vertex)
                        );
                        vertex->state.reset_force();
                        term_function->compute_and_set_forces();
                        SpaceVec mean_force = 0.5 *(force + vertex->state.get_force());

                        integral -= arma::dot(mean_force, delta);
                    }

                    double delta_E = term_function->compute_energy() - energy;
                    if (fabs(delta_E) < 1.e-12) {
                        BOOST_CHECK_SMALL(integral, 1.e-12);
                    }
                    else {
                        BOOST_CHECK_CLOSE(delta_E, integral, precision);
                    }

                    am.move_to(vertex, v_pos0);
                }

                term_function->test_constraints(model.get_logger());
            }

            if (not test_skipped){
                // the machine epsilon has to be scaled to the magnitude of the 
                // values used and multiplied by the desired precision in ULPs 
                // (units in the last place)
                bool test_non_zero = (
                    cumulated_energy > std::numeric_limits<double>::epsilon() 
                                       * cumulated_energy * 10
                );
                test_non_zero = (
                    test_non_zero 
                    // unless the result is subnormal
                    or (cumulated_energy >= std::numeric_limits<double>::min())
                );
                BOOST_TEST(test_non_zero);
                    
            }
        }
    }
    
    BOOST_AUTO_TEST_CASE (test_update_area_elasticity_heterotypic)
    {
        const double precision = 1.e-3;

        Fixture<Case::periodic> fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        const auto& am = model.get_am();
        const auto& cells = am.cells();

        std::uniform_int_distribution<std::size_t> distr(0, cells.size()-1);
        std::size_t cnt_1 = 0;
        while(cnt_1 != cells.size() / 2) {
            auto c_id = distr(*model.get_rng());
            if (cells[c_id]->state.type == 0) {
                cells[c_id]->state.type = 1;
                cnt_1++;
            }
        }
        BOOST_TEST(cnt_1 == cells.size() / 2);

        Config cfg = get_as<Config>(
            "test_update_area_elasticity_heterotypic", 
            fixture.wf_terms
        );

        std::string name = "area_elasticity_heterotypic";
        Config wft_cfg = get_as<Config>(name, cfg);
        std::string term_name = get_as<std::string>("term", wft_cfg, name);
        model.register_work_function_term(term_name, name, wft_cfg);
        auto term_fct = model.get_work_function_term(name);

        std::function<void(std::vector<double>,double)> test_cell_areas =
            [cells](std::vector<double> targets, double prec)
        {
            for (const auto& cell : cells) {
                BOOST_CHECK_CLOSE(
                    cell->state.area_preferential, 
                    targets[cell->state.type],
                    prec
                );
            }
        };
        test_cell_areas(std::vector<double>({1., 1}), precision);

        name = "area_elasticity_heterotypic__set_heterotypic";
        model.get_logger()->info("Updating WF ({})", name);
        term_fct->update_parameters(wft_cfg);
        term_fct->update_parameters(get_as<Config>(name, cfg));
        test_cell_areas(std::vector<double>({1.1, 0.9}), precision);

        name = "area_elasticity_heterotypic__incr_heterotypic";
        model.get_logger()->info("Updating WF ({})", name);
        term_fct->update_parameters(wft_cfg);
        term_fct->update_parameters(get_as<Config>(name, cfg));
        test_cell_areas(std::vector<double>({1.1, 0.9}), precision);

        name = "area_elasticity_heterotypic__incr_heterotypic__restrain";
        model.get_logger()->info("Updating WF ({})", name);
        term_fct->update_parameters(wft_cfg);
        term_fct->update_parameters(get_as<Config>(name, cfg));
        test_cell_areas(std::vector<double>({1.1, 0.9}), precision);
    }
    
    BOOST_AUTO_TEST_CASE (test_update_area_elasticity_gradient)
    {
        Fixture<Case::periodic> fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        const auto& am = model.get_am();
        const auto& cells = am.cells();

        std::uniform_int_distribution<std::size_t> distr(0, cells.size()-1);
        std::size_t cnt_1 = 0;
        while(cnt_1 != cells.size() / 3) {
            auto c_id = distr(*model.get_rng());
            if (cells[c_id]->state.type == 0) {
                cells[c_id]->state.type = 1;
                cnt_1++;
            }
        }
        BOOST_TEST(cnt_1 == cells.size() / 3);

        Config cfg = get_as<Config>(
            "test_update_area_elasticity_gradient", 
            fixture.wf_terms
        );

        std::string name = "area_elasticity_gradient";
        Config wft_cfg = get_as<Config>(name, cfg);
        std::string term_name = get_as<std::string>("term", wft_cfg, name);

        model.register_work_function_term(term_name, name, wft_cfg);
        auto term_fct = model.get_work_function_term(name);

        for (const auto& cell : cells) {
            BOOST_CHECK_CLOSE(cell->state.area_preferential, 1., 1.e-10);
        }

        name = "area_elasticity_gradient__increase";
        for (std::size_t it = 1; it <= 3; it++) {
            model.get_logger()->info("Updating WF ({})", name);
            term_fct->update_parameters(get_as<Config>(name, cfg));
            double max = std::numeric_limits<double>::min();
            double min = std::numeric_limits<double>::max();
            double grad_mean = 0.;
            double max_ = std::numeric_limits<double>::min();
            double min_ = std::numeric_limits<double>::max();
            double global_mean = 0.;
            double mean = 1.;
            std::size_t cnt = 0;
            for (const auto& cell : cells) {
                BOOST_TEST(cell->state.area_preferential > 0.1999999999);
                global_mean += cell->state.area_preferential;
                if (cell->state.type == 0) {
                    max = std::max(max, cell->state.area_preferential);
                    min = std::min(min, cell->state.area_preferential);
                    grad_mean += cell->state.area_preferential;
                    cnt++;
                }
                else {
                    max_ = std::max(max_, cell->state.area_preferential);
                    min_ = std::min(min_, cell->state.area_preferential);
                }
            }
            BOOST_CHECK_CLOSE(global_mean / double(cells.size()), 
                              mean + 0.05 * cnt / double(cells.size()),
                              0.5);
            BOOST_CHECK_CLOSE(max, 1. + it * 0.2, 0.001);
            BOOST_CHECK_CLOSE(min, 1 - it * 0.1, 0.001);
            BOOST_CHECK_CLOSE(grad_mean / double(cnt), 1. + it * 0.05, 1);
            if (it == 1) {
                BOOST_CHECK_CLOSE(max_, 1., 1.e-10);
                BOOST_CHECK_CLOSE(min_, 1., 1.e-10);
            }

            // update to relax total area
            for (std::size_t i = 0; i < 1000; i++) {
                term_fct->update(0.);
            }
            max = std::numeric_limits<double>::min();
            min = std::numeric_limits<double>::max();
            grad_mean = 0.;
            global_mean = 0.;
            max_ = std::numeric_limits<double>::min();
            min_ = std::numeric_limits<double>::max();
            for (const auto& cell : cells) {
                BOOST_TEST(cell->state.area_preferential > 0.1999999999);
                global_mean += cell->state.area_preferential;
                if (cell->state.type == 0) {
                    max = std::max(max, cell->state.area_preferential);
                    min = std::min(min, cell->state.area_preferential);
                    grad_mean += cell->state.area_preferential;
                }
                else {
                    max_ = std::max(max_, cell->state.area_preferential);
                    min_ = std::min(min_, cell->state.area_preferential);
                }
            }
            BOOST_CHECK_CLOSE(global_mean / double(cells.size()), 1., 1);
            BOOST_CHECK_CLOSE(max, 1. + it * 0.2, 0.001);
            BOOST_CHECK_CLOSE(min, 1 - it * 0.1, 0.001);
            BOOST_CHECK_CLOSE(grad_mean / double(cnt), 1. + it * 0.05, 1);

            mean = global_mean / double(cells.size());
        }
    }
    
    BOOST_AUTO_TEST_CASE (test_update_edge_contractility_polar)
    {
        const double precision = 1.e-4;

        Fixture<Case::non_periodic> fixture;
        auto& model = fixture.vertex_model;
        model.prolog();

        const auto& am = model.get_am();
        const auto& cells = am.cells();

        std::uniform_int_distribution<std::size_t> distr(0, cells.size()-1);
        std::size_t cnt_1 = 0;
        while(cnt_1 != cells.size() / 2) {
            auto c_id = distr(*model.get_rng());
            if (cells[c_id]->state.type == 0) {
                cells[c_id]->state.type = 1;
                cnt_1++;
            }
        }
        BOOST_TEST(cnt_1 == cells.size() / 2);

        Config cfg = get_as<Config>(
            "test_update_edge_contractility_polar", 
            fixture.wf_terms
        );

        std::string name = "edge_contractility_polar";
        Config wft_cfg = get_as<Config>(
            name, 
            cfg
        );
        std::string term_name = get_as<std::string>("term", wft_cfg, name);

        model.register_work_function_term(term_name, name, wft_cfg);

        auto term_fct = model.get_work_function_term(name);

        double energy = model.get_energy();
        double dE = 0.;
        for (std::size_t i = 0; i < 100; i++) {
            double E = model.get_energy();
            
            term_fct->update(0.);

            dE += model.get_energy() - E;
            BOOST_TEST(model.get_energy() - E < -1.e-12);
        }

        BOOST_TEST(model.get_energy() - energy < -1.e-12);
        BOOST_CHECK_CLOSE(model.get_energy() - energy, dE, precision);
    }
    
    BOOST_AUTO_TEST_CASE (test_update_volume_elasticity)
    {
        const double precision = 1.e-6;

        Fixture<Case::periodic> fixture;
        auto& model = fixture.vertex_model;
        const auto& logger = model.get_logger();
        model.prolog();

        // erase all work function terms
        while (model.get_work_function_terms().size() > 0) {
            auto [name, term] = *model.get_work_function_terms().begin();
            model.erase_work_function_term(name);
        }


        const auto& am = model.get_am();
        const auto& cells = am.cells();

        Config cfg = get_as<Config>(
            "test_update_volume_elasticity", 
            fixture.wf_terms
        );

        std::string name = "volume_elasticity";
        Config wft_cfg = get_as<Config>(name, cfg);
        std::string term_name = get_as<std::string>("term", wft_cfg, name);
        model.register_work_function_term(term_name, name, wft_cfg);
        auto term_fct = model.get_work_function_term(name);

        BOOST_TEST(term_fct->test_constraints(logger));

        const auto& cell = cells[cells.size() / 2];
        cell->state.type = 1;

        std::string _height = term_fct->get_name() 
            + "_" + get_as<std::string>("height_parameter_name", wft_cfg);
        std::string _height_derivative = _height + "_derivative";

        double E0 = term_fct->compute_energy(am.vertices(), am.edges(), cells);

        double integral = 0.;
        for (std::size_t i = 0; i < 100; i++) {
            cell->state.update_parameter(_height_derivative, 0.);
            term_fct->compute_and_set_forces();

            double dW_dH = cell->state.get_parameter(_height_derivative);

            double tmp_H = cell->state.get_parameter(_height);
            cell->state.update_parameter(_height, tmp_H - 0.1 / 100.);

            double dH = cell->state.get_parameter(_height) - tmp_H;
            cell->state.update_parameter(_height_derivative, 0.);
            term_fct->compute_and_set_forces();
            double mean_deriv = 0.5 *(dW_dH + cell->state.get_parameter(_height_derivative));

            integral += mean_deriv * dH;
        }

        double E1 = term_fct->compute_energy(am.vertices(), am.edges(), cells);
        BOOST_CHECK_CLOSE(E1 - E0, integral, precision);

        BOOST_TEST(term_fct->test_constraints(logger));


    }

    BOOST_AUTO_TEST_CASE(test_energy_periodic)
    {
        Utopia::PseudoParent pp("test_periodic.yml");
        TEST_PCPVertex_energy test_model("PCPVertex", pp);

        test_model.run_test();
    }

    BOOST_AUTO_TEST_CASE(test_energy)
    {
        Utopia::PseudoParent pp("test.yml");
        TEST_PCPVertex_energy test_model("PCPVertex", pp);

        test_model.prolog();

        test_model.perform_test();

        test_model.epilog();
    }


BOOST_AUTO_TEST_SUITE_END()
