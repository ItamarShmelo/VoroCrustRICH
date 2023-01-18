#include "VoroCrustFace.hpp"
#include <iostream>
#include <sstream>

                                                         current_normal(),
                                                         isPatched(false) {}

void VoroCrustFace::addEdge(Edge edge){
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