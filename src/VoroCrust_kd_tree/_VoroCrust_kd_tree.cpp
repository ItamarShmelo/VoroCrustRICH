#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "VoroCrust_kd_tree.hpp"

void bind_vorocrust_kd_tree(pybind11::module& m);

void bind_vorocrust_kd_tree(pybind11::module& m) {
    using namespace pybind11::literals;

    pybind11::class_<VoroCrust_KD_Tree, VoroCrust_KD_TreePtr>(m, "VoroCrust_KD_Tree")
        .def(pybind11::init<>())
        .def(pybind11::init<std::vector<Vector3D> const&>(), pybind11::kw_only(), "points"_a)
        .def("clear", &VoroCrust_KD_Tree::clear)
        .def("makeTree", &VoroCrust_KD_Tree::makeTree, pybind11::kw_only(), "points"_a)
        .def("remakeTree", &VoroCrust_KD_Tree::remakeTree)
        .def("insert", &VoroCrust_KD_Tree::insert, pybind11::kw_only(), "point"_a)
        .def("empty", &VoroCrust_KD_Tree::empty)
        .def("size", &VoroCrust_KD_Tree::size);
    
    pybind11::class_<VoroCrust_KD_Tree_Boundary, VoroCrust_KD_Tree_BoundaryPtr, VoroCrust_KD_Tree>(m, "VoroCrust_KD_Tree_Boundary")
        .def(pybind11::init<>())
        .def(pybind11::init<std::vector<Vector3D> const&>(), pybind11::kw_only(), "points"_a)
        .def(pybind11::init<std::vector<Vector3D> const&,
                            std::vector<Vector3D> const&,
                            std::vector<std::size_t> const&,
                            std::vector<std::size_t> const&>(), 
                                pybind11::kw_only(),
                                 "points"_a, 
                                 "vecs"_a, 
                                 "feature_index"_a, 
                                 "plc_index"_a);

    pybind11::class_<VoroCrust_KD_Tree_Ball, VoroCrust_KD_Tree_BallPtr, VoroCrust_KD_Tree_Boundary>(m, "VoroCrust_KD_Tree_Ball")
        .def(pybind11::init<>())
        .def(pybind11::init<std::vector<Vector3D> const&,
                            std::vector<Vector3D> const&,
                            std::vector<std::size_t> const&,
                            std::vector<std::size_t> const&,
                            std::vector<double> const&>(), 
                                pybind11::kw_only(),
                                 "points"_a, 
                                 "vecs"_a, 
                                 "feature_index"_a, 
                                 "plc_index"_a,
                                 "radii"_a);
    
    
}

PYBIND11_MODULE(_VoroCrust_kd_tree, m) {
    m.doc() = "VoroCrust KD-Tree class";

    bind_vorocrust_kd_tree(m);
}