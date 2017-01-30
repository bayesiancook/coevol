cmake_minimum_required(VERSION 2.8)

file(GLOB SRC
    "src/*.hpp"
    "src/*.cpp"
)

add_executable(PoissonGamma ${SRC})
