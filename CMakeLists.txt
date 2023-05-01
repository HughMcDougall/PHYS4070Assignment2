cmake_minimum_required(VERSION 3.24)
project(Assignment_2)

set(CMAKE_CXX_STANDARD 17)

#add_executable(Assignment_2 Part_1/simulation.cpp Part_1/vector_utils.cpp Part_1/params.hpp Part_1/forces_and_integrators.cpp Part_1/forces_and_integrators.hpp Part_1/sysvec_utils.cpp Part_1/sysvec_utils.hpp Part_1/manybody_example.hpp)
add_executable(Assignment_2 Part_2/singlegrid_parallel.cpp Part_2/monte_carlo.cpp Part_2/rand_utils.cpp Part_2/vector_utils.cpp Part_2/vector_utils.hpp Part_2/monte_carlo.cpp)