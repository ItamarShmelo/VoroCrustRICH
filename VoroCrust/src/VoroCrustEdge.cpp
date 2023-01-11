#include "VoroCrustEdge.hpp"
#include <iostream>
#include <sstream>

VoroCrustEdge::VoroCrustEdge(std::shared_ptr<VoroCrustVertex> const& v1, 
                             std::shared_ptr<VoroCrustVertex> const& v2, 
                             std::size_t const index_) : vertex1(v1), vertex2(v2), faces(), index(index_){}

bool VoroCrustEdge::checkIfEqual(std::shared_ptr<VoroCrustVertex> const& v1, std::shared_ptr<VoroCrustVertex> const& v2){
    if (v1->index == vertex1->index && v2->index == vertex2->index)
        return true;

    if (v2->index == vertex1->index && v1->index == vertex2->index)
        return true;

    return false;
}

void VoroCrustEdge::addFace(std::shared_ptr<VoroCrustFace> new_face){
    if(faces.size() == 2){
        std::cout << "ERROR : can't have more than 2 faces sharing an edge" << std::endl;
        exit(1);
    }

    faces.push_back(new_face);
}

std::string VoroCrustEdge::repr() {
    std::ostringstream s;

    s << "\n\t\t\tvertex " << vertex1->index << ": " << vertex1->repr() << ", vertex " << vertex2->index << " : " << vertex2->repr();

    return s.str();
}