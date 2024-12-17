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

# set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -O0 -g3 -gdwarf-3 -Wall -Wextra -Wpedantic -Wcast-align -Wunused -Wold-style-cast -Wpointer-arith -Wcast-qual -Wno-missing-braces -Wmissing-declarations -Wno-long-long -Wfatal-errors")
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -O0 -g3 -gdwarf-3")
# set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g -O0")
# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp")

message("CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
message("CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE})
message("CMAKE_CXX_FLAGS_DEBUG = " ${CMAKE_CXX_FLAGS_DEBUG})
message("Root_dir = " ${PROJECT_ROOT_DIR})

file(GLOB SRC_CPP "${PROJECT_SOURCE_DIR}/PLC/PL_Complex.cpp" 
                  "${PROJECT_SOURCE_DIR}/PLC/VoroCrustEdge.cpp" 
                  "${PROJECT_SOURCE_DIR}/PLC/VoroCrustFace.cpp" 
                  "${PROJECT_SOURCE_DIR}/PLC/VoroCrustVertex.cpp" 
                  "${PROJECT_SOURCE_DIR}/RMPS/SliverDriver.cpp" 
                  "${PROJECT_SOURCE_DIR}/RMPS/CornersRMPS.cpp" 
                  "${PROJECT_SOURCE_DIR}/RMPS/EdgesRMPS.cpp" 
                  "${PROJECT_SOURCE_DIR}/RMPS/FacesRMPS.cpp" 
                  "${PROJECT_SOURCE_DIR}/VoroCrust_kd_tree/VoroCrust_kd_tree.cpp" 
                  "main.cpp" "${PROJECT_SOURCE_DIR}/runVoroCrust.cpp" 
                  "${PROJECT_SOURCE_DIR}/trees/trees.cpp" 
                  "${PROJECT_SOURCE_DIR}/vorocrust_vtk/vorocrust_vtk.cpp" 
                  "${PROJECT_SOURCE_DIR}/VoroCrustAlgorithm.cpp")

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

set(EXE_NAME "VoroCrust.out")
message("executable name ${EXE_NAME}")

add_executable(${EXE_NAME} ${SRC_CPP})
target_link_libraries(${EXE_NAME} PRIVATE ${VTK_LIBRARIES} gomp pthread)

vtk_module_autoinit(
    TARGETS ${EXE_NAME}
    MODULES ${VTK_LIBRARIES}
    )
    

message("SRC FILES:")
foreach(file IN LISTS SRC_CPP)
    message(${file})
endforeach()

message("project src directory" ${PROJECT_SOURCE_DIR})
message("pybind11 directory: " ${pybind11_DIR})
# pybind modules 


find_package(pybind11 REQUIRED)
set(CPP_PYBIND_MODULES_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")

add_subdirectory(src)