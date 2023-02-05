#include "CornersRMPS.hpp"
#include <boost/random.hpp>

CornersRMPS::CornersRMPS(double const maxRadius_, double const L_Lipschitz_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), eligble_corners() {}


void CornersRMPS::loadCorners(std::vector<Vertex> const& sharp_corners){
    eligble_corners.clear();
    for(Vertex const& corner_ptr : sharp_corners)
        eligble_corners.push_back(corner_ptr->vertex);
}
