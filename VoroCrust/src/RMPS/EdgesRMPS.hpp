#ifndef EDGES_RMPS
#define EDGES_RMPS

#include "../PLC/PL_Complex.hpp"
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"
#include "../VoroCrust_kd_tree/VoroCrust_kd_tree.hpp"
#include "../trees.hpp"
#include <vector>
#include <array>
#include <boost/random.hpp>

struct EligbleEdge {
    std::array<Vector3D, 2> edge;
    std::size_t crease_index;

    EligbleEdge() : edge(), crease_index(0) {}
    EligbleEdge(Vector3D const& v1, Vector3D const& v2, std::size_t const crease_index_) : edge({v1, v2}), crease_index(crease_index_) {}

    Vector3D const& operator [] (std::size_t const index) const {
        return edge[index];
    }

    Vector3D& operator [] (std::size_t const index)  {
        return edge[index];
    }
};

class EdgesRMPS {
};

#endif // EDGES_RMPS