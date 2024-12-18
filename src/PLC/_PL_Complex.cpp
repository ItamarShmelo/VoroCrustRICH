#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "PL_Complex.hpp"
#include "VoroCrustVertex.hpp"
#include "VoroCrustEdge.hpp"
#include "VoroCrustFace.hpp"

void bind_pl_complex(pybind11::module& m);

void bind_pl_complex(pybind11::module& m) {
    using namespace pybind11::literals;

    pybind11::class_<PL_Complex, PL_ComplexPtr>(m, "PL_Complex")
        .def(pybind11::init<std::vector<Vector3D> const&>(), pybind11::kw_only(), "vertices"_a)
        .def("addFace", &PL_Complex::addFace, pybind11::kw_only(), "indices"_a)
        .def("checkAllVerticesAreUnique", &PL_Complex::checkAllVerticesAreUnique)
        .def("checkAllVerticesAreOnFace", &PL_Complex::checkAllVerticesAreOnFace)
        .def("detectSharpFeatures", &PL_Complex::detectFeatures, pybind11::kw_only(), "sharpTheta"_a)
        .def("getBoundingBox", &PL_Complex::getBoundingBox)
        .def_readonly("sharp_corners", &PL_Complex::sharp_corners)
        .def_readonly("sharp_edges", &PL_Complex::sharp_edges)
        .def_readonly("faces", &PL_Complex::faces);

    pybind11::class_<VoroCrustVertex, VoroCrustVertexPtr>(m, "VoroCrustVertex")
        .def(pybind11::init<Vector3D const&,
                            std::size_t const>(),
                                pybind11::kw_only(), 
                                "vertex"_a,
                                "index"_a)
        .def_readonly("index", &VoroCrustVertex::index)
        .def_readonly("vertex", &VoroCrustVertex::vertex);
    
    pybind11::class_<VoroCrustEdge, VoroCrustEdgePtr>(m, "VoroCrustEdge")
        .def(pybind11::init<VoroCrust::Vertex const&,
                            VoroCrust::Vertex const&,
                            std::size_t const>(),
                                pybind11::kw_only(), 
                                "vertex1"_a, 
                                "vertex2"_a, 
                                "index"_a)
        .def_readonly("index",   &VoroCrustEdge::index)
        .def_readonly("vertex1", &VoroCrustEdge::vertex1)
        .def_readonly("vertex2", &VoroCrustEdge::vertex2);

    pybind11::class_<VoroCrustFace, VoroCrustFacePtr>(m, "VoroCrustFace")
        .def(pybind11::init<std::vector<VoroCrust::Vertex> const&,
                            std::size_t const>(),
                                pybind11::kw_only(), 
                                "vertices"_a, 
                                "index"_a)
        .def_readonly("index",    &VoroCrustFace::index)
        .def_readonly("vertices", &VoroCrustFace::vertices)
        .def_readonly("edges",    &VoroCrustFace::edges);

}

PYBIND11_MODULE(_PL_Complex, m) {
    m.doc() = "PL_Complex class";

    bind_pl_complex(m);
}