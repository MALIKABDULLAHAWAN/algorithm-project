cmake_minimum_required(VERSION 3.15)
project(ccp_graph_sssp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
  target_compile_definitions(OpenMP::OpenMP_CXX INTERFACE HAS_OPENMP)
endif()

# Sources
add_library(graph STATIC
  src/graph.cpp
  src/graph_io.cpp
)

target_include_directories(graph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(bench
  src/bench.cpp
  src/dijkstra_seq.cpp
  src/delta_stepping_omp.cpp
)

target_link_libraries(bench PRIVATE graph)

if(OpenMP_CXX_FOUND)
  target_link_libraries(bench PRIVATE OpenMP::OpenMP_CXX)
endif()

# Windows MSVC warnings level
if(MSVC)
  target_compile_options(graph PRIVATE /W4)
  target_compile_options(bench PRIVATE /W4)
else()
  target_compile_options(graph PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(bench PRIVATE -Wall -Wextra -Wpedantic)
endif()
