#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>

#include "Vector3D.hpp"

void bind_vector_3d(pybind11::module& m);

void bind_vector_3d(pybind11::module& m) {
    using namespace pybind11::literals;
    using Vector3DPtr = std::shared_ptr<Vector3D>;
    
    pybind11::class_<Vector3D, Vector3DPtr>(m, "Vector3D")
        .def(pybind11::init<>())
        .def(pybind11::init<double, double, double>(), pybind11::kw_only(), "x"_a, "y"_a, "z"_a)
        .def(pybind11::init<Vector3D const&>(), pybind11::kw_only(), "other"_a)
        .def_readonly("x", &Vector3D::x)
        .def_readonly("y", &Vector3D::y)
        .def_readonly("z", &Vector3D::z)
        .def("Set", &Vector3D::Set, pybind11::kw_only(), "x"_a, "y"_a, "z"_a);
}

PYBIND11_MODULE(_Vector3D, m) {
    m.doc() = "Vector3D class";

    bind_vector_3d(m);
}