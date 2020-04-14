#define BOOST_TEST_MODULE PCPVertexTestNew

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "../PCPVertex.hh"
#include "../energy.hh"
#include "../algorithm.hh"
#include "../operations.hh"
#include "../PCPVertex_write_tasks.hh"

using namespace Utopia::Models::PCPVertex;
using Utopia::Models::PCPVertex::DataIO::time_energy_adaptor;

auto model_factory(bool periodic) {
    using namespace Utopia::Models::PCPVertex::DataIO;

    if (periodic) {
        Utopia::PseudoParent pp("./test_periodic.yml");
        return PCPVertex("PCPVertex", pp,
            // the energy adaptors
            time_energy_adaptor);
    }
    else {
        Utopia::PseudoParent pp("./test.yml");
        return PCPVertex("PCPVertex", pp,
            // the energy adaptors
            time_energy_adaptor);
    }
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
    :
        log([]() {
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
    { }

    // Teardown, invoked after each test
    ~ModelFixture () {
        log->info("Tearing down ...");
        spdlog::drop_all();
        std::remove("test_data.h5");
    }
};

BOOST_FIXTURE_TEST_SUITE (test_PCPVertex_transitions, ModelFixture)

    void test_divide_cell(bool periodic) {
        auto model = model_factory(periodic);
        model.run();
    }

    BOOST_AUTO_TEST_CASE(test_non_periodic) {
        test_divide_cell(false);
    }

    BOOST_AUTO_TEST_CASE(test_periodic) {
        test_divide_cell(true);
    }

BOOST_AUTO_TEST_SUITE_END()
