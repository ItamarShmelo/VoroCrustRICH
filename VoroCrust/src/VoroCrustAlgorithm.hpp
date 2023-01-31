#ifndef VOROCRUST_ALGORITHM
#define VOROCRUST_ALGORITHM

#include "PLC/PL_Complex.hpp"
#include "trees.hpp"
#include <sstream>
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

        std::vector<EligbleVertex> eligble_vertices;

        VoroCrustAlgorithm( PL_Complex const& plc_,
                            double const sharpTheta_,
                            double const flatTheta_,
                            double const maxRadius_,
                            double const L_Lipschitz_);

        ~VoroCrustAlgorithm() = default;

        /*! \brief runs the VoroCrust Algorithm*/
        void run();

        /*! \brief samples a vertex from the `eligble_vertices` and returns it 
            together with its index then erase it from the vector */
        std::pair<unsigned int, EligbleVertex> sampleEligbleVertices();

        /*! \brief RMPS on the vertices, creates the initial balls on the sharp corners*/
        void RMPS_Vertices();

        /*! \brief calculates the initial radius of a vertex ball */
        double calculateInitialRadiusOfVertex(EligbleVertex const& vertex);

        /*! \brief enforces the Lipschitzness for a strata ball_tree*/
        void enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree);

        std::string repr() const;
};


#endif /* VOROCRUST_ALGORITHM */