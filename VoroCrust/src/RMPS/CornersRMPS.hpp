#ifndef CORNER_RMPS
#define CORNER_RMPS

#include "../PLC/PL_Complex.hpp"
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"
#include "../VoroCrust_kd_tree/VoroCrust_kd_tree.hpp"

using EligbleCorner = Vector3D;

class CornersRMPS {
    public:
        double const maxRadius;
        double const L_Lipschitz;
        double const sharpTheta;
        std::shared_ptr<PL_Complex const> plc;
        std::vector<EligbleCorner> eligble_corners;

        CornersRMPS(double const maxRadius_, double const L_Lipschitz_);
        ~CornersRMPS() = default;

        void loadCorners(std::vector<Vertex> const& sharp_corners);

        void doSampling(VoroCrust_KD_Tree_Ball &corner_ball_tree, VoroCrust_KD_Tree_Boundary &corner_boundary_tree);

    private:
        std::pair<std::size_t, EligbleCorner> newSample();

        double calculateInitialRadius(EligbleCorner const& corner, VoroCrust_KD_Tree_Ball const& corner_ball_tree, VoroCrust_KD_Tree_Boundary const& corner_boundary_tree);
};


#endif // CORNER_RMPS