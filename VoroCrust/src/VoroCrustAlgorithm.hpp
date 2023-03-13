#ifndef VOROCRUST_ALGORITHM
#define VOROCRUST_ALGORITHM

#include "PLC/PL_Complex.hpp"
#include "trees.hpp"
#include "RMPS/CornersRMPS.hpp"
#include "RMPS/EdgesRMPS.hpp"
#include "RMPS/FacesRMPS.hpp"
#include "RMPS/SliverDriver.hpp"
#include <sstream>
#include <memory>


class VoroCrustAlgorithm {
    public:
        std::shared_ptr<PL_Complex> plc;
        Trees trees;

        double const sharpTheta;
        double const flatTheta;
        double const maxRadius;
        double const L_Lipschitz;
        double const alpha;
        std::size_t const maximal_num_iter;

        CornersRMPS cornersDriver;
        EdgesRMPS edgesDriver;
        FacesRMPS facesDriver;
        SliverDriver sliverDriver;

        VoroCrustAlgorithm( PL_Complex const& plc_,
                            double const sharpTheta_,
                            double const flatTheta_,
                            double const maxRadius_,
                            double const L_Lipschitz_,
                            double const alpha_);

        ~VoroCrustAlgorithm() = default;

        /*! \brief runs the VoroCrust Algorithm*/
        void run();

        std::vector<Seed> getSeeds() const;

        std::string repr() const;

        std::pair<std::vector<Seed>, std::vector<Seed>> determineIfSeedsAreInsideOrOutside(std::vector<Seed> const& seeds) const;

        std::pair<std::vector<Vector3D>, std::vector<Vector3D>> calcVolumeSeedsUniform(std::vector<Seed> const& seeds, std::size_t const num_points_x, std::size_t const num_points_y, std::size_t const num_points_z) const;
    
    private:
        /*! \brief enforces the Lipschitzness for a strata ball_tree
            \return true if some ball shrunk
        */
        bool enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree);

        bool sliverElimination();
};

VoroCrust_KD_Tree_Ball makeSeedBallTree(std::vector<Seed> const& seeds);


#endif /* VOROCRUST_ALGORITHM */