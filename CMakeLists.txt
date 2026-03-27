# cmake_minimum_required (VERSION 3.8)

# set(PROJECT_NAME "APPROX")
#     project(${PROJECT_NAME})

# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # file(GLOB SOURCES "*.cpp")
# file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
# add_executable(${PROJECT_NAME} ${SOURCES})

# if (CMAKE_VERSION VERSION_GREATER 3.12)
#   set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 20)
# endif()

# set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/CPET_SEAL/install/lib/cmake/SEAL-4.1")

# find_package(SEAL 4.1 REQUIRED)
# target_link_libraries(${PROJECT_NAME} PRIVATE SEAL::seal) 

# set(MY_MODULE "EBModule")
# find_package(pybind11 REQUIRED)
# pybind11_add_module(${MY_MODULE} src/bindings.cpp src/error_bound.cpp)
# target_include_directories(${MY_MODULE} PRIVATE ${CMAKE_SOURCE_DIR}/include)

cmake_minimum_required(VERSION 3.15)
project(APPROX LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(FILTER APP_SOURCES EXCLUDE REGEX ".*/bindings\\.cpp$")
add_executable(${PROJECT_NAME} ${APP_SOURCES})

set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/CPET_SEAL/install/lib/cmake/SEAL-4.1")
find_package(SEAL 4.1 REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE SEAL::seal)

set(MY_MODULE "EBModule")

set(pybind11_DIR
    "${CMAKE_SOURCE_DIR}/.venv/lib/python3.12/site-packages/pybind11/share/cmake/pybind11"
)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(${MY_MODULE}
    src/eval_plain.cpp
    src/bindings.cpp
    src/parse.cpp
    src/error_bound.cpp
    src/polyEval_class.cpp
)

target_include_directories(${MY_MODULE} PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(${MY_MODULE} PRIVATE SEAL::seal)