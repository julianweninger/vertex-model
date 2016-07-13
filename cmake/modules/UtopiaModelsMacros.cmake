# New CMake policy stack
cmake_policy(PUSH)

# Use <Package>_ROOT variables in respective find_package functions/modules
if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif ()

# Utopia, for sure!
find_package(Utopia REQUIRED)

# FFTW3 for SavannaHeterogeneous model
find_package(FFTW3 3.3 REQUIRED)

# Doxygen for the C++ documentation
find_package(Doxygen
             OPTIONAL_COMPONENTS dot)

# Remove latest policy stack
cmake_policy(POP)
