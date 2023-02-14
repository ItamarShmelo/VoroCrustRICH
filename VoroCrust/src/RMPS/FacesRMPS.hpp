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

class FacesRMPS {
    public:
        double const maxRadius;
        double const L_Lipschitz;
        double const alpha;
        double const sharpTheta;
        double const rejection_probability = 0.1;

        boost::variate_generator<boost::mt19937, boost::uniform_01<>> uni01_gen;

        std::shared_ptr<PL_Complex const> plc;
        std::vector<EligbleFace> eligble_faces;

        FacesRMPS(double const maxRadius_, double const L_Lipschitz_, double const alpha_, double const sharpTheta_, std::shared_ptr<PL_Complex> const& plc_);
        ~FacesRMPS() = default;
};
#endif // FACES_RMPS
