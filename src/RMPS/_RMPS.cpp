#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "CornersRMPS.hpp"
#include "EdgesRMPS.hpp"
#include "FacesRMPS.hpp"
#include "SliverDriver.hpp"

void bind_rmps(pybind11::module& m){
    using namespace pybind11::literals;

    pybind11::class_<Seed, SeedPtr>(m, "Seed")
        .def(pybind11::init<>())
        .def(pybind11::init<Vector3D const&, 
                            double const>(), pybind11::kw_only(), "position"_a, "radius"_a)
        .def_readonly("p", &Seed::p)
        .def_readonly("radius", &Seed::radius);

    pybind11::class_<CornersRMPS, CornersRMPSPtr>(m, "CornersRMPS")
        .def(pybind11::init<double const,
                            double const,
                            double const,
                            PL_ComplexPtr>(),
                                pybind11::kw_only(), 
                                "maxRadius"_a, 
                                "L_Lipschitz"_a, 
                                "sharpTheta"_a, 
                                "plc"_a)
        .def("loadCorners", &CornersRMPS::loadCorners, pybind11::kw_only(), "sharp_corners"_a)
        .def("doSampling", &CornersRMPS::doSampling, pybind11::kw_only(), "corner_ball_tree"_a, "trees"_a);

    pybind11::class_<EdgesRMPS, EdgesRMPSPtr>(m, "EdgesRMPS")
        .def(pybind11::init<double const,
                            double const,
                            double const,
                            double const,
                            PL_ComplexPtr>(),
                                pybind11::kw_only(), 
                                "maxRadius"_a,
                                "L_Lipschitz"_a, 
                                "alpha"_a,
                                "sharpTheta"_a, 
                                "plc"_a)
        .def("loadEdges", &EdgesRMPS::loadEdges, pybind11::kw_only(), "sharp_edges"_a)
        .def("doSampling", &EdgesRMPS::doSampling, pybind11::kw_only(), "edges_ball_tree"_a, "trees"_a);

    pybind11::class_<FacesRMPS, FacesRMPSPtr>(m, "FacesRMPS")
        .def(pybind11::init<double const,
                            double const,
                            double const,
                            double const,
                            PL_ComplexPtr>(),
                                pybind11::kw_only(), 
                                "maxRadius"_a, 
                                "L_Lipschitz"_a, 
                                "alpha"_a,
                                "sharpTheta"_a,
                                "plc"_a)
        .def("loadFaces", &FacesRMPS::loadFaces, pybind11::kw_only(), "sharp_faces"_a)
        .def("doSampling", &FacesRMPS::doSampling, pybind11::kw_only(), "faces_ball_tree"_a, "trees"_a);
    
    pybind11::class_<SliverDriver, SliverDriverPtr>(m, "SliverDriver")
        .def(pybind11::init<double const>(), pybind11::kw_only(), "L_Lipschitz"_a)
        .def("eliminateSlivers", &SliverDriver::eliminateSlivers, pybind11::kw_only(), "trees"_a)
        .def("getSeeds", &SliverDriver::getSeeds, pybind11::kw_only(), "trees"_a);
}

PYBIND11_MODULE(_RMPS, m){
    m.doc() = "RMPS module";

    bind_rmps(m);
}