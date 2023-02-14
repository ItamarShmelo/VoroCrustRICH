#ifndef FACES_RMPS
#define FACES_RMPS

#include <vector>
#include <boost/random.hpp>

#include "../PLC/PL_Complex.hpp"
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"
#include "../VoroCrust_kd_tree/VoroCrust_kd_tree.hpp"
#include "../trees.hpp"

struct EligbleFace {
    std::vector<Vector3D> face;
    std::size_t patch_index;
    std::size_t plc_index;

    EligbleFace() :  face(), patch_index(0), plc_index(0) {}
    EligbleFace(std::vector<Vector3D> const& face_, std::size_t const patch_index_, std::size_t const plc_index_) : face(face_), patch_index(patch_index_), plc_index(plc_index_) {}

    Vector3D const& operator [] (std::size_t const index) const {
        return face[index];
    }

    Vector3D& operator [] (std::size_t const index) {
        return face[index];
    }  
};


#endif // FACES_RMPS
