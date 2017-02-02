cmake_minimum_required(VERSION 2.8)

include_directories("src/core")

file(GLOB CORE
  "src/core/*.cpp"
  )

add_executable(PoissonGamma "src/PoissonGamma.cpp" ${CORE})
set(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic")

file (GLOB DIFFSEL
  "src/diffsel/*.cpp"
  )

# add_executable(DiffSel ${DIFFSEL} ${CORE})
