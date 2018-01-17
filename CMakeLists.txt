cmake_minimum_required(VERSION 3.8)
project(Labb3)

set(CMAKE_CXX_STANDARD 11)
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
SET(GCC_COVERAGE_COMPILE_FLAGS "-Ofast")
SET(files_Labb3 ${Labb3_headers} ${Labb3_sources})


set(SOURCE_FILES main.cpp halfEdgeMesh.cpp halfEdgeMesh.h tdMathLibNew.h mem_block.cpp mem_block.h memHalfEdgeMesh.cpp memHalfEdgeMesh.h)
add_executable(Labb3 ${SOURCE_FILES})