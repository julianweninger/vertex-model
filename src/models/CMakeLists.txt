# Instructions for adding a model to Utopia:
#   1. Create a subfolder with the name of your model, <name>
#   2. Add the corresponding `add_subdirectory(<name>)` below
#   3. Add a CMakeLists.txt in the model directory
#   4. Call `add_model(<name> <name>.cc)` there, where <name> is the name of
#      the model directory. For consistency and identification of the models,
#      it is important that the names of folders and targets match.
#   5. Re-configure using `cmake ..`

# Declare this directory as containing all models of this project
set(UTOPIA_MODELS_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE STRING "the directory containing all model source directories" FORCE)

message(STATUS "Starting Utopia model registration ...")

# Register the model directories
add_subdirectory(ActiveSLS)
add_subdirectory(BendSplay)
add_subdirectory(Collier)
add_subdirectory(CopyMeVertex)
add_subdirectory(NotchDelta)
add_subdirectory(PCPTopology)
add_subdirectory(PCPVertex)
# ...
# NOTE Add your subdirectory above (alphabetically sorted)

message(STATUS "Model registration finished.")
