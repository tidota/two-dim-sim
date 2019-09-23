SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS_DEF} -g -Wall")

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TWO_DIM_SIM)

find_package(yaml-cpp REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(twodimsim main.cpp)
target_link_libraries(twodimsim yaml-cpp Eigen3::Eigen)

#configure_file(settings.yaml settings.yaml COPYONLY)
add_custom_target(settings ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/settings.yaml settings.yaml)
