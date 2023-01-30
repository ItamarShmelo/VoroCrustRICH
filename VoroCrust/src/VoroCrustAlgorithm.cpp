#include "VoroCrustAlgorithm.hpp"
#include <cmath>
#include <iostream>
#include <boost/random.hpp>

VoroCrustAlgorithm::VoroCrustAlgorithm( PL_Complex const& plc_, 
                                        double const sharpTheta_, 
                                        double const flatTheta_, 
                                        double const maxRadius_,
                                        double const L_Lipschitz_): plc(plc_), 
                                                                    sharpTheta(sharpTheta_),
                                                                    flatTheta(flatTheta_),
                                                                    maxRadius(maxRadius_),
                                                                    L_Lipschitz(L_Lipschitz_),
                                                                    trees(){

    if(sharpTheta > M_PI_2){
        std::cout << "ERROR: sharpTheta > pi/2" << std::endl;
        exit(1);
    }           

    if(L_Lipschitz >= 1){
        std::cout << "ERROR: L_Lipschitz >= 1" << std::endl;
        exit(1);
    }

    if(L_Lipschitz <= 0){
        std::cout << "ERROR: L_Lipschitz <= 0 " << std::endl;
        exit(1);
    }

}

void VoroCrustAlgorithm::run() {
    
    if(not plc.checkAllVerticesAreOnFace()) exit(1);

    if(not plc.checkIfALLFacesAreFlat()) exit(1);

    plc.detectFeatures(sharpTheta, flatTheta);

    //! TODO: make sampling size a user input!
    trees.loadPLC(plc, 1e5, 1e6);    

    //! TODO: init eligable edges vertices and faces
    eligble_vertices = plc.sharp_corners;
    RMPS_Vertices();
    enforceLipschitzness(trees.ball_kd_vertices);
    std::cout << "\nRun enforceLipschitzness again\n--------------\n" << std::endl;
    enforceLipschitzness(trees.ball_kd_vertices);
    trees.ball_kd_vertices.remakeTree();
    for(std::size_t iteration = 0; iteration < maximal_num_iter; ++iteration){
        
    }

}

std::pair<unsigned int, Vertex> VoroCrustAlgorithm::sampleEligbleVertices(){
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_int_distribution<> int_distribution(0, eligble_vertices.size());
    boost::variate_generator<boost::mt19937, boost::random::uniform_int_distribution<>> rand_gen(rng, int_distribution);

    unsigned int index = rand_gen();
    std::pair<unsigned int, Vertex> pair(index, eligble_vertices[index]);
    
    // erase vertex from eligable vertices after it was sampled.
    eligble_vertices.erase(eligble_vertices.begin()+index);
    
    return pair;
}

void VoroCrustAlgorithm::RMPS_Vertices(){
    while(not eligble_vertices.empty()){
        std::pair<unsigned int, Vertex> const sample = sampleEligbleVertices();
        double radius;
        if(not trees.ball_kd_vertices.points.empty())
            radius = calculateInitialRadiusOfVertex(sample.second);
        else
            radius = maxRadius;

        trees.ball_kd_vertices.insert(sample.second->vertex, radius);
    }

    trees.ball_kd_vertices.remakeTree();
}

double VoroCrustAlgorithm::calculateInitialRadiusOfVertex(Vertex const& vertex){

    int nearsetBall_index = trees.ball_kd_vertices.nearestNeighbor(vertex->vertex);

    Vector3D const& nearsetBallCenter = trees.ball_kd_vertices.points[nearsetBall_index];
    double const dist_q = distance(vertex->vertex, nearsetBallCenter);
    double const r_q = trees.ball_kd_vertices.ball_radii[nearsetBall_index];

    int nearsetSharpCorner_index = trees.VC_kd_sharp_corners.kNearestNeighbors(vertex->vertex, 2)[1];
    Vector3D const& nearestSharpCorner = trees.VC_kd_sharp_corners.points[nearsetSharpCorner_index];

    double const dist_q_prime = distance(vertex->vertex, nearestSharpCorner);
    
    return std::min<double>({maxRadius, 0.49*dist_q_prime, r_q + L_Lipschitz*dist_q});
}

void VoroCrustAlgorithm::enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree){
    std::size_t num_of_points = ball_tree.points.size();
    for(std::size_t i = 0; i<num_of_points; ++i){
        for(std::size_t j = 0; j<num_of_points; ++j){
            if(i == j) continue;
            double const dist = distance(ball_tree.points[i], ball_tree.points[j]);
            if(ball_tree.ball_radii[i] > ball_tree.ball_radii[j] + L_Lipschitz*dist){
                std::cout << "Enforce Lipschitzness ball " << i << std::endl;
            }
            ball_tree.ball_radii[i] = std::min({ball_tree.ball_radii[i], ball_tree.ball_radii[j] + L_Lipschitz*dist}); 
        }
    }
}

std::string VoroCrustAlgorithm::repr() const {
    std::ostringstream s;
    
    s << "VoroCrustAlgorithm : \n--------------------------------\n\n";
    s << "PLC : \n------------\n" << plc.repr() << std::endl;

    return s.str();
}