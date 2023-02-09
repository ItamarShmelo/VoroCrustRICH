#include "VoroCrustAlgorithm.hpp"
#include <cmath>
#include <iostream>

VoroCrustAlgorithm::VoroCrustAlgorithm( PL_Complex const& plc_, 
                                        double const sharpTheta_, 
                                        double const flatTheta_, 
                                        double const maxRadius_,
                                        double const L_Lipschitz_,
                                        double const alpha_): plc(plc_), 
                                                                    sharpTheta(sharpTheta_),
                                                                    flatTheta(flatTheta_),
                                                                    maxRadius(maxRadius_),
                                                                    L_Lipschitz(L_Lipschitz_),
                                                                    alpha(alpha_),
                                                                    trees(),
                                                                    maximal_num_iter(15),
                                                                    cornersDriver(maxRadius_, L_Lipschitz_),
                                                                    edgesDriver(maxRadius, L_Lipschitz_, alpha_, sharpTheta_) {

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
    //! IMPORTANT: sampling size can effect the convergence of the algorithm because the radius is determined using proximity to the sampled points on different features. Make sure that the sampling size is compatible to the size of the smallest polygon in the data. One wants the sampling to be "dense" in the edges and faces.
    trees.loadPLC(plc, 1e6, 1e6);    

    //! TODO: init eligable edges vertices and faces
    cornersDriver.loadCorners(plc.sharp_corners);
    edgesDriver.loadEdges(plc.sharp_edges);

    cornersDriver.doSampling(trees.ball_kd_vertices, trees.VC_kd_sharp_corners);
    
    enforceLipschitzness(trees.ball_kd_vertices);
    
    std::cout << "\nRun enforceLipschitzness again\n--------------\n" << std::endl;
    enforceLipschitzness(trees.ball_kd_vertices);
    
    trees.ball_kd_vertices.remakeTree();

    for(std::size_t iteration = 0; iteration < maximal_num_iter; ++iteration){
        edgesDriver.doSampling(trees.ball_kd_edges, trees);
        break;
    }

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