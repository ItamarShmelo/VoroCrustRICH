#include "PL_Complex.hpp"

PL_Complex::PL_Complex(std::vector<Vector3D> const& vertices_) : vertices(), faces(){
    
    for(auto& vertex : vertices_) 
        vertices.push_back(std::make_shared<VoroCrustVertex>(vertex));
}


std::string PL_Complex::repr() const{
    std::ostringstream s;
    int face_index = 0;
    for(auto& face : faces){
        face_index++;
        s << "face " << face_index << ": " << face->repr() << "\n";
    }

    return s.str();
}