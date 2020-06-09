#define BOOST_TEST_MODULE NotchDelta test

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "../NotchDelta.hh"
#include "../NotchDelta_write_tasks.hh"

using namespace Utopia;
using namespace Utopia::Models::NotchDelta;


/// Factory for NotchDelta model
NotchDelta model_factory() {
    // initialize pp inside to avoid problems with datamanager and memory
    // access
    PseudoParent pp("./test.yml");

    using Utopia::Models::NotchDelta::DataIO::density_time;

    return NotchDelta("NotchDelta", pp, {}, std::make_tuple(density_time));
}

/// A fixture used in the test_PCPVertex test suite
/** Besides destructing the model it is also necessary to destruct the logger
 *  and to free the datapath to construct a new model from the same config
 *
 *  NOTE the models output_path needs to be "test_data.h5"
 *
 */
struct ModelFixture {
    std::shared_ptr<spdlog::logger> log;

    // Construct
    ModelFixture ()
        : log([](){
            auto logger = spdlog::get("test");

            // Create it only if it does not already exist
            if (not logger) {
                logger = spdlog::stdout_color_mt("test");
            }

            // Set level and global logging pattern
            logger->set_level(spdlog::level::debug);
            spdlog::set_pattern("[%T.%e] [%^%l%$] [%n]  %v");
            // "[HH:MM:SS.mmm] [level(colored)] [logger]  <message>"

            return logger;
        }())
    {}

    // Teardown, invoked after each test
    ~ModelFixture () {
        log->info("Tearing down ...");
        spdlog::drop_all();
        std::remove("test_data.h5");
    }
};

// ++ Tests +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
BOOST_FIXTURE_TEST_SUITE (test_NotchDelta, ModelFixture)

/// Test the usage of the custom_neighborhood
BOOST_AUTO_TEST_CASE(test_custom_neighborhood)
{
    auto model = model_factory();
    auto cm = model.get_cm();

    // set up custom neighborhood with 10 inactive cells
    std::uniform_int_distribution<> prob_distr(0, cm.cells().size()-11);

    unsigned int it;
    for (it = 0; it < cm.cells().size() - 10; it++) {
        auto c = cm.cells()[it];
        c->custom_links().neighbors.clear();
        // add 8 random neighbors
        for (int i = 0; i < 8; i++) {
            auto n = cm.cells()[prob_distr(*model.get_rng())];
            while (n == c) {
                n = cm.cells()[prob_distr(*model.get_rng())];
            }
            c->custom_links().neighbors.push_back(n);
        }
    }
    // the inactive cells
    for (void(); it < cm.cells().size(); it++) {
        auto c = cm.cells()[it];
        c->custom_links().neighbors.clear();
        c->state.cell_type = NotchDelta::CellType::inactive;
    }

    model.run();
    
    for (it = 0; it < cm.cells().size() - 10; it++) {
        auto c = cm.cells()[it];
        BOOST_TEST(c->state.cell_type != NotchDelta::CellType::inactive);
        BOOST_TEST(c->state.cell_type != NotchDelta::CellType::progenitor);
        BOOST_TEST(c->custom_links().neighbors.size() == 8);
        for (auto n : c->custom_links().neighbors) {
            BOOST_TEST(n->state.cell_type != NotchDelta::CellType::inactive);
        }
    }
    // the inactive cells
    for (void(); it < cm.cells().size(); it++) {
        auto c = cm.cells()[it];
        BOOST_TEST(c->state.cell_type == NotchDelta::CellType::inactive);
    }

    // compare the number of contacts
    unsigned int num_hh_contacts = 0;
    for (const auto& cell : cm.cells()) {
        if (cell->state.cell_type != NotchDelta::CellType::hair) {
            continue;
        }

        bool has_contact = false;
        for (auto n : cell->custom_links().neighbors) {
            if (n->state.cell_type == NotchDelta::CellType::hair) {
                has_contact = true;
                break;
            }
        }
        if (has_contact) {
            num_hh_contacts++;
        }
    }
    BOOST_TEST(num_hh_contacts == model.get_hh_contacts());
}

BOOST_AUTO_TEST_SUITE_END()
