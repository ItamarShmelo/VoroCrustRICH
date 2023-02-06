#include "EdgesRMPS.hpp"

EdgesRMPS::EdgesRMPS(double const maxRadius_, double const L_Lipschitz_, double const alpha_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), alpha(alpha_), uni01_gen(boost::mt19937(std::time(nullptr)), boost::random::uniform_01<>()) {}
