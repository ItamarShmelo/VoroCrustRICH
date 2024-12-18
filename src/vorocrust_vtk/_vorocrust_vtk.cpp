#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "vorocrust_vtk.hpp"

namespace vorocrust_vtk {
void bind_vorocrust_vtk(pybind11::module& m);

void bind_vorocrust_vtk(pybind11::module& m){
    using namespace pybind11::literals;

    m.def("write_vtu_PL_Complex", &write_vtu_PL_Complex, pybind11::kw_only(), "filename"_a, "plc"_a);

    m.def("write_vtu_trees", &write_vtu_trees, pybind11::kw_only(), "filename"_a, "trees"_a);
    
    m.def("write_ballTree", &write_ballTree, pybind11::kw_only(), "filename"_a, "ball_tree"_a);

    m.def("write_points", &write_points, pybind11::kw_only(), "filename"_a, "points"_a);
}

} // namespace vorocrust_vtk

PYBIND11_MODULE(_vorocrust_vtk, m) {
    m.doc() = "VoroCrust VTK utilities";

    vorocrust_vtk::bind_vorocrust_vtk(m);
}

