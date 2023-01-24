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

std::string VoroCrustAlgorithm::repr() const {
    std::ostringstream s;
    
    s << "VoroCrustAlgorithm : \n--------------------------------\n\n";
    s << "PLC : \n------------\n" << plc.repr() << std::endl;

    return s.str();
}