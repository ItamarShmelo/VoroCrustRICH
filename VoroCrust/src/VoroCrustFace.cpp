#include "VoroCrustFace.hpp"

VoroCrustFace::VoroCrustFace(std::vector<std::shared_ptr<VoroCrustVertex>> const& vertices_, 
                             std::size_t const index_) : vertices(vertices_), index(index_) {}


Vector3D VoroCrustFace::calcNormal(){
    Vector3D const& v1 = vertices[1]->vertex - vertices[0]->vertex;
    Vector3D const& v2 = vertices[2]->vertex - vertices[1]->vertex;

    Vector3D const& cross_product = CrossProduct(v1, v2);
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