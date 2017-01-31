cmake_minimum_required(VERSION 2.8)

file(GLOB SRC
    "src/*.cpp"
)

add_executable(PoissonGamma ${SRC})
set(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic")
