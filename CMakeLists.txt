cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(mini-coevol CXX)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic")

include_directories("include")

# Utils static library
set(UTILS src/utils/Chrono.cpp src/utils/linalg.cpp)
add_library (coevol_utils STATIC ${UTILS})

# Core static library
file(GLOB CORE "src/core/*.cpp")
add_library (coevol_core STATIC ${CORE})
target_compile_features(coevol_core INTERFACE cxx_auto_type)
target_link_libraries(coevol_core coevol_utils)

# PoissonGamma executable
add_executable(poisson_gamma "src/PoissonGamma.cpp")
target_link_libraries(poisson_gamma coevol_core)

# Diffsel things
file (GLOB DIFFSEL "src/diffsel/*.cpp")
add_executable(diffsel "src/DirichletNormalCodonUsageSelectionMS.cpp" ${DIFFSEL})
target_link_libraries(diffsel coevol_core)
target_compile_features(diffsel PRIVATE cxx_auto_type)
