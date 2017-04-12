cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(mini-coevol CXX)

# Compilation options
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic") # normal mode
if (CMAKE_BUILD_TYPE MATCHES PERF)
  message(STATUS "Compiling in performance measurement mode")
  set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic -g -fno-omit-frame-pointer")
endif()
if (CMAKE_BUILD_TYPE MATCHES PERFFULL)
  message(STATUS "Compiling in fine-grain performance measurement mode")
  set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic -g -fno-omit-frame-pointer -fno-inline-functions -fno-optimize-sibling-calls")
endif()
if (CMAKE_BUILD_TYPE MATCHES COVERAGE)
  message(STATUS "Compiling in coverage measurement mode")
  set(CMAKE_CXX_FLAGS "-O0 -fprofile-arcs -fprofile-instr-generate -fcoverage-mapping -g")
endif()

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
  src/core/Chain.cpp
  src/core/Correlation.cpp
  src/core/DAGnode.cpp
  src/core/Move.cpp
  src/core/ProbModel.cpp
  src/core/RandomTypes.cpp)
add_library (coevol_core STATIC ${CORE})
target_compile_features(coevol_core PUBLIC ${CXX_FEATURES})
target_link_libraries(coevol_core coevol_utils)

# Phylogeny
set(PHYLOGENY
  src/phylogeny/BranchSitePath.cpp
  src/phylogeny/CodonSequenceAlignment.cpp
  src/phylogeny/CodonStateSpace.cpp
  src/phylogeny/CodonSubMatrix.cpp
  src/phylogeny/GTRSubMatrix.cpp
  src/phylogeny/MSCodonSubMatrix.cpp
  src/phylogeny/PhyloProcess.cpp
  src/phylogeny/RandomBranchSitePath.cpp
  src/phylogeny/SequenceAlignment.cpp
  src/phylogeny/SiteMapping.cpp
  src/phylogeny/StateSpace.cpp
  src/phylogeny/SubMatrix.cpp
  src/phylogeny/Tree.cpp
  src/phylogeny/TaxonSet.cpp)
add_library (coevol_phylogeny STATIC ${PHYLOGENY})
target_compile_features(coevol_phylogeny PUBLIC ${CXX_FEATURES})
target_link_libraries(coevol_phylogeny coevol_core)


########################
#     EXECUTABLES      #
########################

# Diffsel
add_executable(diffsel "app/DirichletNormalCodonUsageSelectionMS.cpp")
target_link_libraries(diffsel coevol_phylogeny)
target_compile_features(diffsel PRIVATE ${CXX_FEATURES})

# Tracecomp
add_executable(tracecomp "app/TraceComp.cpp")
target_link_libraries(tracecomp coevol_phylogeny)
target_compile_features(tracecomp PRIVATE ${CXX_FEATURES})


########################
#        TESTS         #
########################
enable_testing()

file(GLOB APP_SOURCES test/*.cpp) # all files in test/
foreach(testsourcefile ${APP_SOURCES})
  # get filename without extension to use as executable name
  get_filename_component(testname ${testsourcefile} NAME_WE)
  message(STATUS "Found test " ${testname} ": " ${testsourcefile}) # pretty message :)

  # add executable
  add_executable(${testname} ${testsourcefile})
  target_link_libraries(${testname} coevol_core)
  target_compile_features(${testname} PRIVATE ${CXX_FEATURES})

  # add to tests
  add_test(${testname} ${testname})
endforeach(testsourcefile ${APP_SOURCES})
