#ifndef VOROCRUST_ALGORITHM
#define VOROCRUST_ALGORITHM

#include "PLC/PL_Complex.hpp"
#include "trees.hpp"
#include <sstream>

class VoroCrustAlgorithm {
    public:
        PL_Complex plc;
        Trees trees;

        double const sharpTheta;
        double const flatTheta;
        double const maxRadius;
        double const L_Lipschitz;

        std::vector<Vertex> eligble_vertices;

        VoroCrustAlgorithm( PL_Complex const& plc_,
                            double const sharpTheta_,
                            double const flatTheta_,
                            double const maxRadius_,
                            double const L_Lipschitz_);

        ~VoroCrustAlgorithm() = default;


        void run();

        std::pair<unsigned int, Vertex> sampleEligbleVertices();

        void RMPS_Vertices();

        double calculateInitialRadiusOfVertex(Vertex const& vertex);

        void enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree);

        std::string repr() const;
};


#endif /* VOROCRUST_ALGORITHM */