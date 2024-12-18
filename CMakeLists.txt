cmake_minimum_required(VERSION 3.20.2)
project(VoroCrustImplementation)
enable_language(CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

message("CXX Compiler: " ${CMAKE_CXX_COMPILER})
message("C Compiler: " ${CMAKE_C_COMPILER})
message("Fortran Compiler: " ${CMAKE_Fortran_COMPILER})
message("CXX Standard: " ${CMAKE_CXX_STANDARD})
message("Build Type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -O0 -g3 -gdwarf-3 -Wall -Wextra -Wpedantic -Wcast-align -Wunused -Wold-style-cast -Wpointer-arith -Wcast-qual -Wno-missing-braces -Wmissing-declarations -Wno-long-long -Wfatal-errors")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

message("CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
message("CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE})
message("CMAKE_CXX_FLAGS_DEBUG = " ${CMAKE_CXX_FLAGS_DEBUG})
message("Root_dir = " ${PROJECT_ROOT_DIR})

set(VTK_DIRECTORY "/software/x86_64/5.14.0/vtk/9.3.0/Intel/OneApi/2024.0.1/with_mesa")

find_package(VTK COMPONENTS
    CommonColor
    CommonCore
    CommonDataModel
    CommonTransforms
    FiltersGeneral
    FiltersSources
    IOXML 
    IOParallelXML
    # ParallelMPI 
    InteractionStyle
    NO_MODULE
    PATHS ${VTK_DIRECTORY} NO_DEFAULT_PATH
)
if (NOT VTK_FOUND)
    message(FATAL_ERROR "NO VTK")
endif()

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${BOOST_INCLUDE})

message("Boost Path ${BOOST_INCLUDE}")

message("project src directory" ${PROJECT_SOURCE_DIR})
message("pybind11 directory: " ${pybind11_DIR})
# pybind modules 


find_package(pybind11 REQUIRED)
set(CPP_PYBIND_MODULES_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")

add_subdirectory(src)