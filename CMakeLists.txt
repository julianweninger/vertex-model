# Set minimum CMake version to the one supplied with Ubuntu 18.04
cmake_minimum_required(VERSION 3.10)

# Start the project
project(UtopiaModels
    DESCRIPTION "Additional models for the Utopia framework"
    LANGUAGES C CXX
    VERSION 0.1
)

# --- CMake Modules ---

# Insert Utopia macros on top of the module path list
set(module_path ${PROJECT_SOURCE_DIR}/cmake/modules)
list(INSERT CMAKE_MODULE_PATH 0 ${module_path})

# Load the Utopia macros and execute them here
# (This checks dependencies, enables CMake functions, etc.)
include(UtopiaModelsMacros)

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
