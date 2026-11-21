cmake_minimum_required(VERSION 3.16)
project(sciddica_ss3_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python 3.15 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

option(SCIDDICA_PGO_GENERATE "PGO generate phase" OFF)
option(SCIDDICA_PGO_USE "PGO use phase" OFF)

pybind11_add_module(sciddica_ss3_cpp cpp/sciddica_ss3_cpp/module.cpp)
target_include_directories(sciddica_ss3_cpp PRIVATE cpp)
target_compile_options(sciddica_ss3_cpp PRIVATE -O3 -march=native -funroll-loops)

if(SCIDDICA_PGO_GENERATE)
  target_compile_options(sciddica_ss3_cpp PRIVATE -fprofile-generate)
  target_link_options(sciddica_ss3_cpp PRIVATE -fprofile-generate)
elseif(SCIDDICA_PGO_USE)
  target_compile_options(sciddica_ss3_cpp PRIVATE -fprofile-use -fprofile-correction)
  target_link_options(sciddica_ss3_cpp PRIVATE -fprofile-use)
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(sciddica_ss3_cpp PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS sciddica_ss3_cpp LIBRARY DESTINATION .)
