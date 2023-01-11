#include "PL_Complex.hpp"
#include <iostream>
#include "../../source/misc/utils.hpp"

PL_Complex::PL_Complex(std::vector<Vector3D> const& vertices_) : vertices(), faces(){
    std::size_t index = 0;
    for(auto& vertex : vertices_) {
        vertices.push_back(std::make_shared<VoroCrustVertex>(vertex, index));
        index++;
    }
}

std::shared_ptr<VoroCrustEdge> PL_Complex::addEdge(std::shared_ptr<VoroCrustVertex> const& v1, std::shared_ptr<VoroCrustVertex> const& v2){
    for (auto& edge : edges){
        if (edge->checkIfEqual(v1, v2)){
            return edge;
        }
    }

    auto new_edge_ptr = std::make_shared<VoroCrustEdge>(v1, v2, edges.size());
    
    new_edge_ptr->vertex1->addEdge(new_edge_ptr);
    new_edge_ptr->vertex2->addEdge(new_edge_ptr);

    edges.push_back(new_edge_ptr);
    return new_edge_ptr;
}

void PL_Complex::addFace(std::vector<unsigned int> const& indices){

    if (indices.size() != unique(indices).size()){
        std::cout << "ERROR: Repeated indices in PL_Complex::addFace!" << std::endl;
        exit(1);
    }

    if(indices.size() < 3){
        std::cout << "ERROR: Face has to have at least three vertices!" << std::endl;
        exit(1);
    }

    std::vector<std::shared_ptr<VoroCrustVertex>> face_vertices = std::vector<std::shared_ptr<VoroCrustVertex>>();

    for(auto& index : indices)
        face_vertices.push_back(vertices[index]);

    std::shared_ptr<VoroCrustFace> new_face_ptr = std::make_shared<VoroCrustFace>(face_vertices, faces.size());

    for (std::shared_ptr<VoroCrustVertex> vertex_ptr : face_vertices){
        vertex_ptr->addFace(new_face_ptr);
    }

    faces.push_back(new_face_ptr);
}

std::string PL_Complex::repr() const{
    std::ostringstream s;


    for(auto& face : faces){
        s << "Face " << face->index << ": \n" << face->repr() << "\n";
    }

    return s.str();
}