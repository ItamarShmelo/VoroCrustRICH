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
                                                              edgesDriver(maxRadius_, L_Lipschitz_, alpha_, sharpTheta_, plc),
                                                              facesDriver(maxRadius_, L_Lipschitz_, alpha_, sharpTheta_, plc),
                                                              sliverDriver(L_Lipschitz_) {

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
    for(std::size_t iteration = 0; iteration < maximal_num_iter; ++iteration){
        // enfore lipchitzness on vertices
        std::cout << "\nCorners Lipchitzness\n--------------\n" << std::endl;
        enforceLipschitzness(trees.ball_kd_vertices);    
        
        do {
            std::cout << "\nEdgesRMPS\n--------------\n" << std::endl;
            edgesDriver.loadEdges(plc->sharp_edges);
            edgesDriver.doSampling(trees.ball_kd_edges, trees);
            trees.ball_kd_edges.remakeTree();
        } while(enforceLipschitzness(trees.ball_kd_edges));
            
        do {
            std::cout << "\nFacesRMPS\n--------------\n" << std::endl;
            facesDriver.loadFaces(plc->faces);
            facesDriver.doSampling(trees.ball_kd_faces, trees);
            trees.ball_kd_faces.remakeTree();
        } while(enforceLipschitzness(trees.ball_kd_faces));
        
        if(not sliverDriver.eliminateSlivers(trees)) break;

        enforceLipschitzness(trees.ball_kd_edges);
        enforceLipschitzness(trees.ball_kd_faces);
    }
}

bool VoroCrustAlgorithm::enforceLipschitzness(VoroCrust_KD_Tree_Ball& ball_tree){
    std::size_t num_of_points = ball_tree.points.size();
    
    bool isBallsShrunk = false;

    std::cout << "Running enforceLipschitzness " << std::endl;
    // go through all pairs i != j and enforce r_i <= r_j + L * ||p_i - p_j||
    //! MAYBE: only consider the overlapping balls or balls up to radius r_i?
    for(std::size_t i = 0; i<num_of_points; ++i){
        for(std::size_t j = 0; j<num_of_points; ++j){
            if(i == j) continue;
            double const dist = distance(ball_tree.points[i], ball_tree.points[j]);
            
            if(ball_tree.ball_radii[i] > ball_tree.ball_radii[j] + L_Lipschitz*dist){
                std::cout << "Enforce Lipschitzness ball " << i << ", r_old " << ball_tree.ball_radii[i] << ", r_new " << ball_tree.ball_radii[j] + L_Lipschitz*dist << std::endl;
                isBallsShrunk = true;
            }
            
            ball_tree.ball_radii[i] = std::min(ball_tree.ball_radii[i], ball_tree.ball_radii[j] + L_Lipschitz*dist); 
        }
    }

    return isBallsShrunk;
}

std::vector<Vector3D> VoroCrustAlgorithm::getSeeds() const {
    return sliverDriver.getSeeds(trees);
}

std::pair<std::vector<Vector3D>, std::vector<Vector3D>> VoroCrustAlgorithm::determineIfSeedsAreInsideOrOutside(std::vector<Vector3D> const& seeds) const {

    // if(seeds.size() % 2 != 0){
    //     std::cout << "ERROR: seeds is not even" << std::endl;
    //     exit(1);
    // }

    if(seeds.empty()){
        std::cout << "ERROR: seeds is empty" << std::endl;
        exit(1);
    }

    std::vector<Vector3D> in_seeds, out_seeds;

    int i = 0;
    for(Vector3D const& seed : seeds) {
        // determine if a seed is in or out using the ray casting algorithm
        std::cout << "seed num " << ++i << std::endl;
        
        int count = 0;
        for(Face const& face : plc->faces) {
            auto const& [success, p_inter] = face->pointXYaxisRayIntersectsAt(seed);

            if(not success) continue;

            if(p_inter.z > seed.z && face->pointIsInsideFace(p_inter)){
                count++;
            }
        }

        if(count % 2 == 0){
            out_seeds.push_back(seed);
        } else {
            in_seeds.push_back(seed);
        }
    }

    return std::pair<std::vector<Vector3D>, std::vector<Vector3D>>(in_seeds, out_seeds);
}

std::string VoroCrustAlgorithm::repr() const {
    std::ostringstream s;
    
    s << "VoroCrustAlgorithm : \n--------------------------------\n\n";
    s << "PLC : \n------------\n" << plc->repr() << std::endl;

    return s.str();
}