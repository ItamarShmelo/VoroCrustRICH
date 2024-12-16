pybind11_add_module(_Vector3D "_Vector3D.cpp")
target_link_libraries(_Vector3D PRIVATE pybind11::module)
set_target_properties(_Vector3D PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}" LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/../cpp_modules")