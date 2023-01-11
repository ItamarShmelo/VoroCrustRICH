#include "VoroCrustFace.hpp"

VoroCrustFace::VoroCrustFace(std::vector<std::shared_ptr<VoroCrustVertex>> const& vertices_, 
                             std::size_t const index_) : vertices(vertices_), index(index_) {}

std::string VoroCrustFace::repr() const {
    std::ostringstream s;
    s << "\tEdges: ";
    for (auto& edge : edges){
        s << "\n\t\tedge " << edge->index << ": " << edge->repr() << ", ";
    }

    return s.str();
}