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

set(VOROCRUST_VTK_SRC "vorocrust_vtk.cpp"
                      "../VoroCrust_kd_tree/VoroCrust_kd_tree.cpp" 
                      "../PLC/PL_Complex.cpp"
                      "../PLC/VoroCrustVertex.cpp"
                      "../PLC/VoroCrustEdge.cpp"
                      "../PLC/VoroCrustFace.cpp"
                      "../trees/trees.cpp"
                    )

pybind11_add_module(_vorocrust_vtk "_vorocrust_vtk.cpp" ${VOROCRUST_VTK_SRC})
target_link_libraries(_vorocrust_vtk PRIVATE pybind11::module ${VTK_LIBRARIES})
set_target_properties(_vorocrust_vtk PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}" LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")