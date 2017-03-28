cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(mini-coevol CXX)

# Compilation options
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic")
set(CXX_FEATURES cxx_auto_type cxx_nullptr cxx_final cxx_override cxx_range_for)
include_directories("src")


########################
#     STATIC LIBS      #
########################

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


########################
#     EXECUTABLES      #
########################

# Diffsel
file (GLOB DIFFSEL "src/diffsel/*.cpp")
add_executable(diffsel "app/DirichletNormalCodonUsageSelectionMS.cpp" ${DIFFSEL})
target_link_libraries(diffsel coevol_core)
target_compile_features(diffsel PRIVATE ${CXX_FEATURES})


########################
#        TESTS         #
########################
enable_testing()

file(GLOB APP_SOURCES test/*.cpp) # all files in test/
foreach(testsourcefile ${APP_SOURCES})
  # get filename withtout extension to use as executable name
  get_filename_component(testname ${testsourcefile} NAME_WE)
  message(STATUS "Found test " ${testname} ": " ${testsourcefile}) # pretty message :)

  # add executable
  add_executable(${testname} ${testsourcefile})
  target_link_libraries(${testname} coevol_core)
  target_compile_features(${testname} PRIVATE ${CXX_FEATURES})

  # add to tests
  add_test(${testname} ${testname})
endforeach(testsourcefile ${APP_SOURCES})
