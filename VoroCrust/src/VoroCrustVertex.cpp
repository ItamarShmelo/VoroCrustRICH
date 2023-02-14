#include "VoroCrustVertex.hpp"

VoroCrustVertex::VoroCrustVertex(Vector3D const& vertex_) : vertex(vertex_), faces(){}

void VoroCrustVertex::addFace(std::shared_ptr<VoroCrustFace> new_face){
    faces.push_back(new_face);
}

std::string VoroCrustVertex::repr() const {
    std::ostringstream s;
    s << "x = " << vertex.x << ", y = " << vertex.y << ", z = " << vertex.z;
    return s.str();
}
