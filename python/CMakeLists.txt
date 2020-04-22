# Handle configuration and installation of python packages and tests here

# Pass info to utopya  ........................................................
message(STATUS "")
message(STATUS "Providing model information to Utopia frontend ...")

# Declare the directories containing the models tests and plots
set(UTOPIA_PYTHON_MODEL_TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/model_tests
    CACHE STRING "the directory containing all python-based model tests" FORCE)
set(UTOPIA_PYTHON_MODEL_PLOTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/model_plots
    CACHE STRING "the directory containing all python-based model plots" FORCE)

# Register the existing models and associated frontend-related information
register_models_with_frontend()


# Tests .......................................................................
message(STATUS "")
message(STATUS "Adding python test targets ...")

# Add tests for models, accessible via `make test_models_python`
add_custom_target(test_models_python
                  COMMAND ${RUN_IN_UTOPIA_ENV}
                          python -m pytest -v
                          ${CMAKE_CURRENT_SOURCE_DIR}/model_tests/
)

# Add the python model tests as a dependency for the test_models target
add_dependencies(test_models test_models_python)


# Additional python packages ..................................................
# ... from requirements file
message(STATUS "")
message(STATUS "Installing custom python dependencies ...")
execute_process(COMMAND
                    ${RUN_IN_UTOPIA_ENV}
                    pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
                RESULT_VARIABLE exit_code
                OUTPUT_VARIABLE output
                ERROR_VARIABLE output
                )
if (NOT ${exit_code} STREQUAL "0")
    message(WARNING "Error ${exit_code} in installing additional python "
                    "requirements! Install the packages defined in the "
                    "python/requirements.txt file manually to remove this "
                    "warning.\nOutput of pip install was:\n${output}")
endif()


message(STATUS "Python setup complete.\n")
