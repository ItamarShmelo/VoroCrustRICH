#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "trees.hpp"

void bind_trees(pybind11::module& m);

void bind_trees(pybind11::module& m){
    using namespace pybind11::literals;
    
    pybind11::class_<Trees, TreesPtr>(m, "Trees")
        .def(pybind11::init<>())
        .def("loadPLC", &Trees::loadPLC, pybind11::kw_only(), "plc"_a, "Nsample_edges"_a, "Nsample_faces"_a)
        .def_readonly("VC_kd_sharp_corners", &Trees::VC_kd_sharp_corners)
        .def_readonly("VC_kd_sharp_edges", &Trees::VC_kd_sharp_edges)
        .def_readonly("VC_kd_faces", &Trees::VC_kd_faces)
        .def_readonly("ball_kd_vertices", &Trees::ball_kd_vertices)
        .def_readonly("ball_kd_edges", &Trees::ball_kd_edges)
        .def_readonly("ball_kd_faces", &Trees::ball_kd_faces);
}

PYBIND11_MODULE(_trees, m){
    m.doc() = "Trees class";

    bind_trees(m);
}