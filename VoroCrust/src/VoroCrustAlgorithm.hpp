#ifndef VOROCRUST_ALGORITHM
#define VOROCRUST_ALGORITHM

#include "PLC/PL_Complex.hpp"
#include "trees.hpp"
#include <sstream>
#include <boost/random.hpp>
using EligbleVertex = Vector3D;

class VoroCrustAlgorithm {
    public:
        PL_Complex plc;
        Trees trees;

        double const sharpTheta;
        double const flatTheta;
        double const maxRadius;
        double const L_Lipschitz;
        std::size_t const maximal_num_iter;

        boost::variate_generator<boost::mt19937, boost::uniform_01<>> uni01_gen;
        std::vector<EligbleVertex> eligble_vertices;

        VoroCrustAlgorithm( PL_Complex const& plc_,
                            double const sharpTheta_,
                            double const flatTheta_,
                            double const maxRadius_,
                            double const L_Lipschitz_);

        ~VoroCrustAlgorithm() = default;

        /*! \brief runs the VoroCrust Algorithm*/
        void run();

        /*! \brief enforces the Lipschitzness for a strata ball_tree*/
        void enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree);

        std::string repr() const;
};


#endif /* VOROCRUST_ALGORITHM */