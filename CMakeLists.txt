cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(mini-coevol CXX)

# Compilation options
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic")
set(CXX_FEATURES cxx_auto_type cxx_nullptr cxx_final cxx_override cxx_range_for)
include_directories("src")

# Utils static library
set(UTILS src/utils/Chrono.cpp src/utils/linalg.cpp src/utils/Random.cpp)
add_library (coevol_utils STATIC ${UTILS})
target_compile_features(coevol_utils PUBLIC ${CXX_FEATURES})

# Core static library
set(CORE
  src/core/BaseType.cpp
  src/core/DAGnode.cpp
  src/core/Move.cpp
  src/core/ProbModel.cpp
  src/core/RandomTypes.cpp)
add_library (coevol_core STATIC ${CORE})
target_compile_features(coevol_core PUBLIC ${CXX_FEATURES})
target_link_libraries(coevol_core coevol_utils)

# PoissonGamma executable
add_executable(poisson_gamma "src/PoissonGamma.cpp")
target_link_libraries(poisson_gamma coevol_core)
target_compile_features(poisson_gamma PRIVATE ${CXX_FEATURES})

# Diffsel things
file (GLOB DIFFSEL "src/diffsel/*.cpp")
add_executable(diffsel "src/DirichletNormalCodonUsageSelectionMS.cpp" ${DIFFSEL})
target_link_libraries(diffsel coevol_core)
target_compile_features(diffsel PRIVATE ${CXX_FEATURES})

# move 1
add_executable(move_test "src/Test.cpp")
target_link_libraries(move_test coevol_core)
target_compile_features(move_test PRIVATE ${CXX_FEATURES})

# move 2
add_executable(move_test2 "src/Test2.cpp")
target_link_libraries(move_test2 coevol_core)
target_compile_features(move_test2 PRIVATE ${CXX_FEATURES})

# TESTS
enable_testing()
add_test(move1 move_test)
add_test(move2 move_test2)
add_test(poissonGamma poisson_gamma ../data/test.data test.out)
