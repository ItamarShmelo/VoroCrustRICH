cmake_minimum_required(VERSION 3.20.2)
project(VoroCrustImplementation)
enable_language(CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
message("CXX Compiler: " ${CMAKE_CXX_COMPILER})
message("C Compiler: " ${CMAKE_C_COMPILER})
message("Fortran Compiler: " ${CMAKE_Fortran_COMPILER})
message("CXX Standard: " ${CMAKE_CXX_STANDARD})
message("Build Type: " ${CMAKE_BUILD_TYPE})

message("CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
message("CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE})
message("CMAKE_CXX_FLAGS_DEBUG = " ${CMAKE_CXX_FLAGS_DEBUG})
message("CMAKE_Fortran_FLAGS = " ${CMAKE_Fortran_FLAGS})
message("CMAKE_Fortran_FLAGS_RELEASE = " ${CMAKE_Fortran_FLAGS_RELEASE})
message("CMAKE_Fortran_FLAGS_DEBUG = " ${CMAKE_Fortran_FLAGS_DEBUG})
message("Root_dir = " ${PROJECT_ROOT_DIR})

set(proj_dir "${PROJECT_SOURCE_DIR}")

file(GLOB_RECURSE SRC_CPP "src/*.cpp")
file(GLOB R3D "../source/3D/GeometryCommon/Vector3D.cpp")
file(GLOB MISC "../source/misc/utils.cpp" "../source/misc/serializable.cpp" "../source/misc/universal_error.cpp")

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${BOOST_INCLUDE})

message("proj_dir ${proj_dir}")

add_executable("VoroCrust.out" ${SRC_CPP} ${R3D} ${MISC})
