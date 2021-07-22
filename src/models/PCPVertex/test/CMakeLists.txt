# Add model tests.
# They will be available as a custom target:  test_model_${MODEL_NAME}

add_model_tests(# Use consistent capitalization for the model name!
                MODEL_NAME PCPVertex
                # The sources of the model tests to carry out. Each of these
                # will become a test target with the same name.
                SOURCES
                    "test_space.cc"
                    "test.cc"
                    "test_transition.cc"
                # Optional: Files to be copied to the build directory
                AUX_FILES
                    "space_test.yml"
                    
                    # test cases
                    "test.yml"
                    "test_periodic.yml"
                    "test_shape_elastic.yml"
                    "test_columnar.yml"
                    "test_columnar_periodic.yml"
                    
                    # test transitions
                    "test_transitions.yml"
                    "test_transitions_periodic.yml")

