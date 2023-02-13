#include "VoroCrustAlgorithm.hpp"
#include <cmath>
#include <iostream>

VoroCrustAlgorithm::VoroCrustAlgorithm( PL_Complex const& plc_, 
                                        double const sharpTheta_, 
                                        double const flatTheta_, 
                                        double const maxRadius_,
                                        double const L_Lipschitz_,
                                        double const alpha_): plc(std::make_shared<PL_Complex>(plc_)), 
                                                              trees(),
                                                              sharpTheta(sharpTheta_),
                                                              flatTheta(flatTheta_),
                                                              maxRadius(maxRadius_),
                                                              L_Lipschitz(L_Lipschitz_),
                                                              alpha(alpha_),
                                                              maximal_num_iter(15),
                                                              cornersDriver(maxRadius_, L_Lipschitz_, sharpTheta_, plc),
                                                              edgesDriver(maxRadius, L_Lipschitz_, alpha_, sharpTheta_, plc) {

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
    
    //! TODO: Maybe all this need to be in the plc under detect features?
    if(not plc->checkAllVerticesAreOnFace()) exit(1);

    if(not plc->checkIfALLFacesAreFlat()) exit(1);

    plc->detectFeatures(sharpTheta, flatTheta);

    std::cout << plc->repr() << std::endl;
    //! TODO: make sampling size a user input!
    //! IMPORTANT: sampling size can effect the convergence of the algorithm because the radius is determined using proximity to the sampled points on different features. Make sure that the sampling size is compatible to the size of the smallest polygon in the data. One wants the sampling to be "dense" in the edges and faces.
    trees.loadPLC(*plc, 1e5, 1e6);    

    //! TODO: init eligable edges vertices and faces
    cornersDriver.loadCorners(plc->sharp_corners);
    cornersDriver.doSampling(trees.ball_kd_vertices, trees);
    trees.ball_kd_vertices.remakeTree();
    // sliver elimination loop
    bool redoVertices = false;
    for(std::size_t iteration = 0; iteration < maximal_num_iter; ++iteration){
        while(true){
            redoVertices = false;
            // enfore lipchitzness on vertices
            std::cout << "\nCorners Lipchitzness\n--------------\n" << std::endl;
            enforceLipschitzness(trees.ball_kd_vertices);    
            
            do {
                std::cout << "\nEdgesRMPS\n--------------\n" << std::endl;
                edgesDriver.loadEdges(plc->sharp_edges);
                redoVertices = edgesDriver.doSampling(trees.ball_kd_edges, trees);
                trees.ball_kd_edges.remakeTree();
            } while(enforceLipschitzness(trees.ball_kd_edges) && !redoVertices);
            
            if(redoVertices) continue; // if vertices enforce lipschitzness again

            break;
        }
        break;
    }

}

bool VoroCrustAlgorithm::enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree){
    std::size_t num_of_points = ball_tree.points.size();

    bool isBallsShrunk = false;
    // go through all pairs i != j and enforce r_i <= r_j + L * ||p_i - p_j||
    for(std::size_t i = 0; i<num_of_points; ++i){
        for(std::size_t j = 0; j<num_of_points; ++j){
            if(i == j) continue;
            double const dist = distance(ball_tree.points[i], ball_tree.points[j]);
            if(ball_tree.ball_radii[i] > ball_tree.ball_radii[j] + L_Lipschitz*dist){
                std::cout << "Enforce Lipschitzness ball " << i << std::endl;
                isBallsShrunk = true;
            }
            ball_tree.ball_radii[i] = std::min({ball_tree.ball_radii[i], ball_tree.ball_radii[j] + L_Lipschitz*dist}); 
        }
    }

    return isBallsShrunk;
}

std::string VoroCrustAlgorithm::repr() const {
    std::ostringstream s;
    
    s << "VoroCrustAlgorithm : \n--------------------------------\n\n";
    s << "PLC : \n------------\n" << plc->repr() << std::endl;

    return s.str();
}