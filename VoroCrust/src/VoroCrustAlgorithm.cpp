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
                                                                    trees(),
                                                                    maximal_num_iter(15),
                                                                    uni01_gen(  boost::mt19937(std::time(nullptr)), 
                                                                                boost::random::uniform_01<>())  {

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
    for(std::size_t i=0; i<plc.sharp_corners.size(); ++i){
        eligble_vertices.push_back(plc.sharp_corners[i]->vertex);
    }
    RMPS_Vertices();
    enforceLipschitzness(trees.ball_kd_vertices);
    
    std::cout << "\nRun enforceLipschitzness again\n--------------\n" << std::endl;
    enforceLipschitzness(trees.ball_kd_vertices);
    
    trees.ball_kd_vertices.remakeTree();
    for(std::size_t iteration = 0; iteration < maximal_num_iter; ++iteration){
        
    }

}

std::pair<unsigned int, EligbleVertex> VoroCrustAlgorithm::sampleEligbleVertices(){
    // create a random number generator
    boost::mt19937 rng(std::time(nullptr));
    boost::random::uniform_int_distribution<> int_distribution(0, eligble_vertices.size());
    boost::variate_generator<boost::mt19937, boost::random::uniform_int_distribution<>> rand_gen(rng, int_distribution);

    // sample a random index
    unsigned int index = rand_gen();
    std::pair<unsigned int, EligbleVertex> pair(index, eligble_vertices[index]);
    
    // erase vertex from eligable vertices after it was sampled.
    eligble_vertices.erase(eligble_vertices.begin()+index);
    
    return pair;
}

void VoroCrustAlgorithm::RMPS_Vertices(){
    // while there are eligble vertices
    while(not eligble_vertices.empty()){
        // sample a vertex
        std::pair<unsigned int, EligbleVertex> const sample = sampleEligbleVertices();
        
        double radius;
        // determine a radius
        if(not trees.ball_kd_vertices.points.empty())
            radius = calculateInitialRadiusOfVertex(sample.second);
        else
            radius = maxRadius;

        trees.ball_kd_vertices.insert(sample.second, radius);
    }
}

double VoroCrustAlgorithm::calculateInitialRadiusOfVertex(EligbleVertex const& vertex){

    // fined nearest ball center
    int nearsetBall_index = trees.ball_kd_vertices.nearestNeighbor(vertex);

    Vector3D const& nearsetBallCenter = trees.ball_kd_vertices.points[nearsetBall_index];
    double const dist_q = distance(vertex, nearsetBallCenter); // ||p-q||
    double const r_q = trees.ball_kd_vertices.ball_radii[nearsetBall_index]; 

    // find nearest sharp corner
    int nearsetSharpCorner_index = trees.VC_kd_sharp_corners.kNearestNeighbors(vertex, 2)[1];
    Vector3D const& nearestSharpCorner = trees.VC_kd_sharp_corners.points[nearsetSharpCorner_index];

    double const dist_q_prime = distance(vertex, nearestSharpCorner); // ||p-q^*||
    
    return std::min<double>({maxRadius, 0.49*dist_q_prime, r_q + L_Lipschitz*dist_q}); // min {sz, 0.49*||r-q^*||, r_q + L * ||p - q||}
}

void VoroCrustAlgorithm::enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree){
    std::size_t num_of_points = ball_tree.points.size();

    // go through all pairs i != j and enforce r_i <= r_j + L * ||p_i - p_j||
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