# Handle configuration and installation of python packages and tests here

# Pass info to utopya  ........................................................
message(STATUS "")
message(STATUS "Providing model information to Utopia frontend ...")

# Register the existing models
register_models_with_frontend()

# Register paths of the additional python plot modules with the frontend.
# The modules given here are then available as imports under module strings
# starting with `model_plots.<model_name>`.
register_plot_module(NAME ccees_group_models
                     MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/model_plots")


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
    message(SEND_ERROR "Error ${exit_code} in installing additional "
                       "python requirements!\n"
                       "Output of pip install was:\n${output}")
endif()


message(STATUS "Python setup complete.\n")
