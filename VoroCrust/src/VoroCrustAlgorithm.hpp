#ifndef VOROCRUST_ALGORITHM
#define VOROCRUST_ALGORITHM

#include "PLC/PL_Complex.hpp"
#include "trees.hpp"
#include "RMPS/CornersRMPS.hpp"
#include "RMPS/EdgesRMPS.hpp"
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

        VoroCrustAlgorithm( PL_Complex const& plc_,
                            double const sharpTheta_,
                            double const flatTheta_,
                            double const maxRadius_,
                            double const L_Lipschitz_,
                            double const alpha_);

        ~VoroCrustAlgorithm() = default;

        /*! \brief runs the VoroCrust Algorithm*/
        void run();

        /*! \brief enforces the Lipschitzness for a strata ball_tree
            \return true if some ball shrunk
        */
        bool enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree);

        std::string repr() const;
};


#endif /* VOROCRUST_ALGORITHM */