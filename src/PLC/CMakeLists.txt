set(PL_COMPLEX_SRC "PL_Complex.cpp"
                   "VoroCrustVertex.cpp"
                   "VoroCrustEdge.cpp"
                   "VoroCrustFace.cpp"
                   "../VoroCrust_kd_tree/VoroCrust_kd_tree.cpp"
)

pybind11_add_module(_PL_Complex "_PL_Complex.cpp" ${PL_COMPLEX_SRC})
target_link_libraries(_PL_Complex PRIVATE pybind11::module)
set_target_properties(_PL_Complex PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}" LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")