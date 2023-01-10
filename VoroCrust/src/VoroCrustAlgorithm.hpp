#ifndef VOROCRUST_ALGORITHM
#define VOROCRUST_ALGORITHM

#include "PL_Complex.hpp"
#include <sstream>

class VoroCrustAlgorithm {
    public:
        PL_Complex plc;

        double const sharpTheta;
        double const flatTheta;
        double const maxRadius;

        VoroCrustAlgorithm( PL_Complex const& plc_,
                            double const sharpTheta_,
                            double const flatTheta_,
                            double const maxRadius_);

        void run();

        std::string repr() const;
};


#endif /* VOROCRUST_ALGORITHM */