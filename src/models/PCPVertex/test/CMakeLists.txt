# Add model tests.
# They will be available as a custom target:  test_model_${MODEL_NAME}

add_model_tests(# Use consistent capitalization for the model name!
                MODEL_NAME PCPVertex
                # The sources of the model tests to carry out. Each of these
                # will become a test target with the same name.
                SOURCES
                    # "test.cc"  # FIXME Reactivate after fixing
                    "test_transition.cc"
                    "test_space.cc"
                # Optional: Files to be copied to the build directory
                AUX_FILES
                    "space_test.yml"
                    "test_periodic.yml"
                    "test_non_periodic.yml")

