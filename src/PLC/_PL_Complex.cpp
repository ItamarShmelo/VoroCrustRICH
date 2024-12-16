#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "PL_Complex.hpp"

void bind_pl_complex(pybind11::module& m) {
    using namespace pybind11::literals;

    pybind11::class_<PL_Complex, std::shared_ptr<PL_Complex>>(m, "PL_Complex")
        .def(pybind11::init<std::vector<Vector3D> const&>(), pybind11::kw_only(), "vertices"_a)
        .def("addFace", &PL_Complex::addFace, pybind11::kw_only(), "indices"_a)
        .def("checkAllVerticesAreUnique", &PL_Complex::checkAllVerticesAreUnique)
        .def("checkAllVerticesAreOnFace", &PL_Complex::checkAllVerticesAreOnFace)
        .def("detectSharpFeatures", &PL_Complex::detectFeatures, pybind11::kw_only(), "sharpTheta"_a)
        .def("getBoundingBox", &PL_Complex::getBoundingBox);
}

PYBIND11_MODULE(_PL_Complex, m) {
    m.doc() = "PL_Complex class";

    bind_pl_complex(m);
}