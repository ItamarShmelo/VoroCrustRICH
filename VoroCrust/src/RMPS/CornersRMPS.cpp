#include "CornersRMPS.hpp"
#include <boost/random.hpp>

CornersRMPS::CornersRMPS(double const maxRadius_, double const L_Lipschitz_, double const sharpTheta_, std::shared_ptr<PL_Complex> const& plc_) : maxRadius(maxRadius_), L_Lipschitz(L_Lipschitz_), sharpTheta(sharpTheta_), plc(plc_), eligble_corners() {}


void CornersRMPS::loadCorners(std::vector<Vertex> const& sharp_corners){
    eligble_corners.clear();
    for(Vertex const& corner_ptr : sharp_corners)
        eligble_corners.push_back(corner_ptr);
}

void CornersRMPS::doSampling(VoroCrust_KD_Tree_Ball &corner_ball_tree, Trees const& trees){

    Vector3D const empty_vector(0, 0, 0);

    while(not eligble_corners.empty()){
        std::cout << "corner sample " << corner_ball_tree.points.size() << std::endl;

        // sample a vertex
        EligbleCorner const sample = newSample();

        double const radius = calculateInitialRadius(sample, trees);
        
        corner_ball_tree.insert(sample->vertex, empty_vector, radius, corner_ball_tree.points.size(), sample->index);
    }
}

EligbleCorner CornersRMPS::newSample(){
    // create a random number generator to sample from 0 - (eligble_corners.size()-1)
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_int_distribution<> int_distribution(0, eligble_corners.size() - 1);
    boost::variate_generator<boost::mt19937, boost::random::uniform_int_distribution<>> rand_gen(rng, int_distribution);

    // sample a random index;
    std::size_t const index = rand_gen();
    EligbleCorner sample = eligble_corners[index];

    eligble_corners.erase(eligble_corners.begin() + index);

    return sample;
}

double CornersRMPS::calculateInitialRadius(EligbleCorner const& corner, Trees const& trees) const {
    VoroCrust_KD_Tree_Ball const& corner_ball_tree = trees.ball_kd_vertices;

    double r_q = std::numeric_limits<double>::max();
    double dist_q = 0.0;
    
    if(not corner_ball_tree.points.empty()){
        // find nearest ball center
        int nearestBall_index = corner_ball_tree.nearestNeighbor(corner->vertex);

        Vector3D const& nearestBallCenter = corner_ball_tree.points[nearestBall_index];
        dist_q = distance(corner->vertex, nearestBallCenter); //||p-q||
        r_q = corner_ball_tree.ball_radii[nearestBall_index];
    }
    
    double const dist_q_prime = calculateSmoothnessLimitation(corner, trees);

    return std::min<double>({maxRadius, 0.49*dist_q_prime, r_q + L_Lipschitz*dist_q});
}

    // find nearest sharp corner
    int nearestSharpCorner_index = corner_boundary_tree.kNearestNeighbors(corner, 2)[1];
    Vector3D const& nearestSharpCorner = corner_boundary_tree.points[nearestSharpCorner_index];

    double const dist_q_prime = distance(corner, nearestSharpCorner); //||p-q^*||

    return std::min<double>({maxRadius, 0.49*dist_q_prime, r_q + L_Lipschitz*dist_q});
}