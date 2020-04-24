# Add the model target
add_model(NotchDelta NotchDelta.cc)
# NOTE The target should have the same name as the model folder and the *.cc

# Add test directories
add_subdirectory(test EXCLUDE_FROM_ALL)