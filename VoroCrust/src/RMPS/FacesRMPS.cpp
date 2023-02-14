#include "FacesRMPS.hpp"

FacesRMPS::FacesRMPS(double const maxRadius_, 
                     double const L_Lipschitz_, 
                     double const alpha_, 
                     double const sharpTheta_, 
                     std::shared_ptr<PL_Complex> const& plc_) : 
                                    maxRadius(maxRadius_), 
                                    L_Lipschitz(L_Lipschitz_), 
                                    alpha(alpha_), 
                                    sharpTheta(sharpTheta_), 
                                    uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>()), 
                                    plc(plc_), 
                                    eligble_faces() {}