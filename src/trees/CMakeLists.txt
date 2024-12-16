set(TREES_SRC "trees.cpp"
              "../PLC/PL_Complex.cpp"
              "../PLC/VoroCrustVertex.cpp"
              "../PLC/VoroCrustEdge.cpp"
              "../PLC/VoroCrustFace.cpp"
              "../VoroCrust_kd_tree/VoroCrust_kd_tree.cpp"
)

pybind11_add_module(_trees "_trees.cpp" ${TREES_SRC})
target_link_libraries(_trees PRIVATE pybind11::module)
set_target_properties(_trees PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}" LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")