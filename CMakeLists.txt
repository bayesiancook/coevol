cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(mini-coevol CXX)

include_directories("include")

file(GLOB CORE
  "src/core/*.cpp"
  )

file (GLOB DIFFSEL
  "src/diffsel/*.cpp"
  )

add_library (coevol_core STATIC ${CORE})
target_compile_features(coevol_core INTERFACE cxx_auto_type)

add_executable(poisson_gamma "src/PoissonGamma.cpp")
target_link_libraries(poisson_gamma coevol_core)
# set(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic")

add_executable(diffsel "src/DirichletNormalCodonUsageSelectionMS.cpp" ${DIFFSEL})
target_link_libraries(diffsel coevol_core)
target_compile_features(diffsel PRIVATE cxx_auto_type)
