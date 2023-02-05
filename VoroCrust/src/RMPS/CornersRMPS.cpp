#include "CornersRMPS.hpp"
#include <boost/random.hpp>

CornersRMPS::CornersRMPS(double const maxRadius_, double const L_Lipschitz_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), eligble_corners() {}


void CornersRMPS::loadCorners(std::vector<Vertex> const& sharp_corners){
    eligble_corners.clear();
    for(Vertex const& corner_ptr : sharp_corners)
        eligble_corners.push_back(corner_ptr->vertex);
}

double CornersRMPS::calculateInitialRadius(EligbleCorner const& corner, VoroCrust_KD_Tree_Ball const& corner_ball_tree, VoroCrust_KD_Tree_Boundary const& corner_boundary_tree){
    // find nearest ball center
    int nearestBall_index = corner_ball_tree.nearestNeighbor(corner);

    Vector3D const& nearestBallCenter = corner_ball_tree.points[nearestBall_index];
    double const dist_q = distance(corner, nearestBallCenter); //||p-q||
    double const r_q = corner_ball_tree.ball_radii[nearestBall_index];

    // find nearest sharp corner
    int nearestSharpCorner_index = corner_boundary_tree.kNearestNeighbors(corner, 2)[1];
    Vector3D const& nearestSharpCorner = corner_boundary_tree.points[nearestSharpCorner_index];

    double const dist_q_prime = distance(corner, nearestSharpCorner); //||p-q^*||

    return std::min<double>({maxRadius, 0.49*dist_q_prime, r_q + L_Lipschitz*dist_q});
}