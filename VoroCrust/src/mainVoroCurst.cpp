#include "PL_Complex.hpp"
#include "VoroCrustAlgorithm.hpp"
#include <iostream>

int main(int argc, char *argv[]){
    
    std::vector<Vector3D> vertices{Vector3D(0, 0, 0), Vector3D(1, 0, 0), Vector3D(1, 1, 0), Vector3D(0, 1, 0), Vector3D(0, 0, 1)};
    PL_Complex plc = PL_Complex(vertices);
    plc.addFace(std::vector<unsigned int>{0,1,2,3});
    plc.addFace(std::vector<unsigned int>{0,1,4});

    std::cout << plc.repr() << std::endl;

    std::cout << "\n\nchange one vertex" << std::endl;
    plc.vertices[2]->vertex.x = 2;
    std::cout << plc.repr() << std::endl;

    std::cout << "\n\nCheck VoroCrustAlgorithm " << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    VoroCrustAlgorithm alg(plc, 15., 20., 1.);

    std::cout << alg.repr() << std::endl;

    return 0;
}