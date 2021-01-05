# Add the model target
add_model(PCPTopology PCPTopology.cc)
# NOTE The target should have the same name as the model folder and the *.cc

# Add compile options to make release builds run faster
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    # Disable all Armadillo run-time checks (e.g. bounds checking)
    target_compile_definitions(PCPTopology PUBLIC "ARMA_NO_DEBUG")
    message(STATUS "(armadillo runtime checks disabled for this model)")
endif()

# Add test directories
add_subdirectory(test EXCLUDE_FROM_ALL)