set(VOROCRUST_KD_TREE_SRC "VoroCrust_kd_tree.cpp")

pybind11_add_module(_VoroCrust_kd_tree "_VoroCrust_kd_tree.cpp" ${VOROCRUST_KD_TREE_SRC})
target_link_libraries(_VoroCrust_kd_tree PRIVATE pybind11::module)
set_target_properties(_VoroCrust_kd_tree PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}" LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")