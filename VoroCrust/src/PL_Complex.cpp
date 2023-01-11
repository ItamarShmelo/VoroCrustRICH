#include "PL_Complex.hpp"

PL_Complex::PL_Complex(std::vector<Vector3D> const& vertices_) : vertices(), faces(){
    std::size_t index = 0;
    for(auto& vertex : vertices_) {
        vertices.push_back(std::make_shared<VoroCrustVertex>(vertex, index));
        index++;
    }
}

void PL_Complex::addFace(std::vector<unsigned int> const& indices){

    std::vector<std::shared_ptr<VoroCrustVertex>> face_vertices = std::vector<std::shared_ptr<VoroCrustVertex>>();

    for(auto& index : indices)
        face_vertices.push_back(vertices[index]);

    std::shared_ptr<VoroCrustFace> new_face_ptr = std::make_shared<VoroCrustFace>(VoroCrustFace(face_vertices));

    for (std::shared_ptr<VoroCrustVertex> vertex_ptr : face_vertices){
        vertex_ptr->addFace(new_face_ptr);
    }

    faces.push_back(new_face_ptr);
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