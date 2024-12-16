#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "VoroCrustAlgorithm.hpp"

void bind_vorocrust_algorithm(pybind11::module& m) {
    using namespace pybind11::literals;
    
    pybind11::class_<VoroCrustAlgorithm, std::shared_ptr<VoroCrustAlgorithm>>(m, "VoroCrustAlgorithm")
        .def(pybind11::init<PL_Complex const&, 
                            double const,
                            double const, 
                            double const, 
                            double const, 
                            std::size_t const,
                            std::size_t const,
                            std::size_t const>(),
                            pybind11::kw_only(),
                            "plc"_a,
                            "sharpTheta"_a,
                            "maxRadius"_a,
                            "L_Lipschitz"_a,
                            "alpha"_a,
                            "maximal_num_iter"_a,
                            "num_of_samples_edges"_a,
                            "num_of_samples_faces"_a)
        .def("run", &VoroCrustAlgorithm::run)
        .def_readonly("plc", &VoroCrustAlgorithm::plc);
}

PYBIND11_MODULE(_VoroCrustAlgorithm, m) {
    m.doc() = "VoroCrustAlgorithm class";

    bind_vorocrust_algorithm(m);
}