cmake_minimum_required(VERSION 2.8)

include_directories("include")

file(GLOB CORE
  "src/core/*.cpp"
  )

add_library (coevol_core STATIC ${CORE})

add_executable(poisson_gamma "src/PoissonGamma.cpp")
target_link_libraries(poisson_gamma coevol_core)
# set(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic")

file (GLOB DIFFSEL
  "src/diffsel/*.cpp"
  )

add_executable(diffsel "src/DirichletNormalCodonUsageSelectionMS.cpp" ${DIFFSEL})
target_link_libraries(diffsel coevol_core)
