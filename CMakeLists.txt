cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)

project(heat_equation_mpi)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS system mpi filesystem REQUIRED)
find_package(PNG REQUIRED)
find_package(nlohmann_json REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PNG_INCLUDE_DIR})


add_executable(main main.cpp)

target_link_libraries(main ${PNG_LIBRARY} Boost::filesystem Boost::mpi Boost::system ${MPI_LIBRARIES} nlohmann_json::nlohmann_json)