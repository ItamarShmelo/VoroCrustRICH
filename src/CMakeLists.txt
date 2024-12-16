message("In src file ${CPP_PYBIND_MODULES_DIRECTORY}")

add_subdirectory(miscellaneous)
add_subdirectory(PLC)
add_subdirectory(vorocrust_vtk)
add_subdirectory(VoroCrust_kd_tree)
add_subdirectory(trees)

set(VOROCRUST_ALGORITHM_SRC "VoroCrustAlgorithm.cpp" 
                            "VoroCrust_kd_tree/VoroCrust_kd_tree.cpp" 
                            "RMPS/CornersRMPS.cpp"
                            "RMPS/EdgesRMPS.cpp"
                            "RMPS/FacesRMPS.cpp"
                            "RMPS/SliverDriver.cpp"
                            "PLC/PL_Complex.cpp"
                            "PLC/VoroCrustVertex.cpp"
                            "PLC/VoroCrustEdge.cpp"
                            "PLC/VoroCrustFace.cpp"
                            "trees/trees.cpp"
                            )
                            
pybind11_add_module(_VoroCrustAlgorithm "_VoroCrustAlgorithm.cpp" ${VOROCRUST_ALGORITHM_SRC})
target_link_libraries(_VoroCrustAlgorithm PRIVATE pybind11::module)
set_target_properties(_VoroCrustAlgorithm PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}" LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")