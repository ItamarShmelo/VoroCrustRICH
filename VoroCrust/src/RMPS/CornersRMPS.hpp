#ifndef CORNER_RMPS
#define CORNER_RMPS

#include "../PLC/PL_Complex.hpp"
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"
#include "../trees.hpp"
#include "../VoroCrust_kd_tree/VoroCrust_kd_tree.hpp"

using EligbleCorner = Vertex;

class CornersRMPS {
    public:
        double const maxRadius;
        double const L_Lipschitz;
        double const sharpTheta;
        std::shared_ptr<PL_Complex const> plc;
        std::vector<EligbleCorner> eligble_corners; // should be private

        CornersRMPS(double const maxRadius_, double const L_Lipschitz_, double const sharpTheta, std::shared_ptr<PL_Complex> const& plc_);
        ~CornersRMPS() = default;

        //! \brief loads the sharp corners from the plc
        //! should probably use the plc shared ptr instead of the explicit external input
        void loadCorners(std::vector<Vertex> const& sharp_corners);

        //! run the corner ball sampling
        void doSampling(VoroCrust_KD_Tree_Ball &corner_ball_tree, Trees const& trees);

    private:
        //! \brief get a new sample 
        EligbleCorner newSample();

        //! \brief calculate the initial radius given to a ball centered at `corner`
        double calculateInitialRadius(EligbleCorner const& corner, Trees const& trees) const;

        //! \brief calculate the limitation on the radius from co-smoothness condition
        double calculateSmoothnessLimitation(EligbleCorner const& corner, Trees const& trees) const;
        
};


#endif // CORNER_RMPS