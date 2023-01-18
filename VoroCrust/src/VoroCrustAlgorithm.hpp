#ifndef VOROCRUST_ALGORITHM
#define VOROCRUST_ALGORITHM

#include "PLC/PL_Complex.hpp"
#include <sstream>

class VoroCrustAlgorithm {
    public:
        PL_Complex plc;

        double const sharpTheta;
        double const flatTheta;
        double const maxRadius;
        double const L_Lipschitz;

        VoroCrustAlgorithm( PL_Complex const& plc_,
                            double const sharpTheta_,
                            double const flatTheta_,
                            double const maxRadius_,
                            double const L_Lipschitz_);

        ~VoroCrustAlgorithm() = default;

        void run();

        std::string repr() const;
};


#endif /* VOROCRUST_ALGORITHM */