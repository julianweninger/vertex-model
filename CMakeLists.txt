# Set minimum CMake version to the one supplied with Ubuntu 18.04
cmake_minimum_required(VERSION 3.13)

# Start the project
project(UtopiaVertexModel
    DESCRIPTION "The vertex model and associated models within the Utopia framework"
    LANGUAGES C CXX
    VERSION 0.1
)

# Parse CMake options
include(CMakeDependentOption)
# Option is only presented if CMAKE_BUILD_TYPE is not set to 'Release'.
# It always defaults to OFF and will always be OFF for 'Release' builds.
cmake_dependent_option(CPP_COVERAGE "Compile C++ code with coverage flags" OFF
                       "NOT CMAKE_BUILD_TYPE STREQUAL Release" OFF)

# --- CMake Modules ---

# Insert Utopia macros on top of the module path list
set(module_path ${PROJECT_SOURCE_DIR}/cmake/modules)
list(INSERT CMAKE_MODULE_PATH 0 ${module_path})

# Load the Utopia macros and execute them here
# (This checks dependencies, enables CMake functions, etc.)
include(UtopiaModelsMacros)

# add extra flags to debug compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
# NOTE The ${CMAKE_CXX_FLAGS_DEBUG} variable can be set in the cache. This is,
#      for example, used in the CI/CD pipeline to add additional flags to DEBUG


# Add this compile definition for debug builds, this same logic works for
# target_compile_options, target_link_options, etc.
target_compile_definitions(Utopia::utopia INTERFACE
    $<$<CONFIG:Debug>:
        UTOPIA_DEBUG=1
    >
)

# Create symlinks to the Utopia Python virtual environment
include(CreateSymlink)
message(STATUS "Creating symlinks to the Utopia Python virtual environment "
               "located at: ${UTOPIA_ENV_DIR}")
create_symlink(${RUN_IN_UTOPIA_ENV} ${PROJECT_BINARY_DIR}/run-in-utopia-env)
create_symlink(${UTOPIA_ENV_DIR}/bin/activate ${PROJECT_BINARY_DIR}/activate)

# --- Include Config Tree ---

# Enable testing via CTest engine
enable_testing()

# Include subdirectories
add_subdirectory(doc)
add_subdirectory(src/models)
add_subdirectory(python)

# Add a test target to rule them all
add_custom_target(test_all)
add_dependencies(test_all test_models)
