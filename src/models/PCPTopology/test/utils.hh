#ifndef UTOPIA_MODELS_PCPTOPOLOGY_TEST_UTILS_HH
#define UTOPIA_MODELS_PCPTOPOLOGY_TEST_UTILS_HH

#include <assert.h>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <numeric>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

/// A fixture used in the test_PCPVertex test suite
/** Besides destructing the model it is also necessary to destruct the logger
 *  and to free the datapath to construct a new model from the same config
 *
 *  \note the models output_path needs to be "test_data.h5"
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
        std::remove("test_data.h5");
        spdlog::drop_all();
    }
};

std::pair<double, double> get_statistics (std::vector<double>& values)
{
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    double mean = sum / values.size();

    double sq_sum = std::inner_product(values.begin(), values.end(),
                                       values.begin(), 0.0);
    double stddev = std::sqrt(sq_sum / values.size() - mean * mean);

    return std::make_pair(mean, stddev);
}

#endif // UTOPIA_MODELS_PCPVERTEX_TEST_UTILS_HH
