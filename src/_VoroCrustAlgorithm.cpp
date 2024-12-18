#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "VoroCrustAlgorithm.hpp"

void bind_vorocrust_algorithm(pybind11::module& m) {
    using namespace pybind11::literals;
    
    m.def("dumpSeeds", &dumpSeeds, pybind11::kw_only(), "dirname"_a, "seeds"_a);
   
    m.def("determineZoneOfSeeds", &determineZoneOfSeeds, pybind11::kw_only(), "seeds"_a, "zone_plcs"_a);
    
    m.def("enforceLipschitzness", &enforceLipschitzness, pybind11::kw_only(), "ball_tree"_a, "L_Lipschitz"_a);
    
    m.def("randomSampleVolumeSeeds", &randomSampleVolumeSeeds, pybind11::kw_only(), "zones_plcs"_a, "zones_boundary_seeds"_a, "maxSize"_a, "trees"_a, "L_Lipschitz"_a);
}

PYBIND11_MODULE(_VoroCrustAlgorithm, m) {
    m.doc() = "VoroCrustAlgorithm class";

    bind_vorocrust_algorithm(m);
}