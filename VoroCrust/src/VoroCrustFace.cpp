#include "VoroCrustFace.hpp"

VoroCrustFace::VoroCrustFace(std::vector<std::shared_ptr<VoroCrustVertex>> const& vertices_) : vertices(vertices_) {}

std::string VoroCrustFace::repr() const {
    std::ostringstream s;
    int i = 0;
    for (auto& vertex : vertices){
        i++;
        s << "Vertex " << i << ": " << vertex->repr() << ", ";
    }

    return s.str();
}