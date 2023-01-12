#include "VoroCrustFace.hpp"
#include <iostream>

VoroCrustFace::VoroCrustFace(std::vector<std::shared_ptr<VoroCrustVertex>> const& vertices_, 
                             std::size_t const index_) : vertices(vertices_), edges(), neighbors(), index(index_), current_normal(), orientationFixed(false){}

void VoroCrustFace::addEdge(std::shared_ptr<VoroCrustEdge> edge){
    edges.push_back(edge);
}


Vector3D VoroCrustFace::calcNormal(){
    Vector3D const& v1 = vertices[1]->vertex - vertices[0]->vertex;
    Vector3D const& v2 = vertices[2]->vertex - vertices[1]->vertex;

    Vector3D const& cross_product = CrossProduct(v2, v1);
    current_normal = cross_product / abs(cross_product); 
    
    return current_normal;
}

std::string VoroCrustFace::repr() const {
    std::ostringstream s;
    s << "\tEdges: ";
    for (auto& edge : edges){
        s << "\n\t\tedge " << edge->index << ": " << edge->repr() << ", ";
    }

    return s.str();
}