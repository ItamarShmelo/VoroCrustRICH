set(RMPS_SRC "SliverDriver.cpp"
             "CornersRMPS.cpp"
             "EdgesRMPS.cpp"
             "FacesRMPS.cpp"
             "../trees/trees.cpp"
             "../PLC/PL_Complex.cpp"
             "../PLC/VoroCrustVertex.cpp"
             "../PLC/VoroCrustEdge.cpp"
             "../PLC/VoroCrustFace.cpp"
             "../VoroCrust_kd_tree/VoroCrust_kd_tree.cpp"
)

pybind11_add_module(_RMPS "_RMPS.cpp" ${RMPS_SRC})
target_link_libraries(_RMPS PRIVATE pybind11::module)
set_target_properties(_RMPS PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}" LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")