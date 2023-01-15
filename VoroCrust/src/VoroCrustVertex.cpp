#include "VoroCrustVertex.hpp"
#include <sstream>

VoroCrustVertex::VoroCrustVertex(Vector3D const& vertex_, 
                                 std::size_t const index_) : vertex(vertex_), 
                                                             faces(), 
                                                             index(index_),
                                                             isSharp(false){}

void VoroCrustVertex::addFace(Face new_face){
    faces.push_back(new_face);
}

std::string VoroCrustVertex::repr() const {
    std::ostringstream s;
    s << "x = " << vertex.x << ", y = " << vertex.y << ", z = " << vertex.z;
    return s.str();
}
