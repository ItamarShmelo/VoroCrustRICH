#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SliverDriver.hpp"

void bind_rmps(pybind11::module& m){
    using namespace pybind11::literals;

    pybind11::class_<Seed, SeedPtr>(m, "Seed")
        .def(pybind11::init<>())
        .def(pybind11::init<Vector3D const&, double const>(), pybind11::kw_only(), "position"_a, "radius"_a)
        .def_readonly("p", &Seed::p)
        .def_readonly("radius", &Seed::radius);
}

PYBIND11_MODULE(_RMPS, m){
    m.doc() = "RMPS module";

    bind_rmps(m);
}