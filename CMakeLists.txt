cmake_minimum_required(VERSION 2.8.4)

project(TSP)

set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -flto -Wall -fPIC -fopenmp")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX   "-isystem ")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(LINKER_OPTIONS                  -flto -Wl,--no-as-needed)

FILE(GLOB SOURCE_FILES src/*.cpp src/*/*.c src/*/*.cpp src/*/*/*.c)

find_package(Cplex)
include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})

add_executable( TSP ${SOURCE_FILES})


target_include_directories(TSP PRIVATE
						"${CMAKE_CURRENT_SOURCE_DIR}/src/svm"
						"${CMAKE_CURRENT_SOURCE_DIR}/src/linear_svm"
						"${CMAKE_CURRENT_SOURCE_DIR}/src/linear_svm/blas"
						)

target_link_libraries(  TSP PRIVATE ${CPLEX_LIBRARIES})
target_link_libraries(  TSP PRIVATE dl)

