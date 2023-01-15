#include "VoroCrustEdge.hpp"
#include <iostream>
#include <sstream>
#include <cmath>

VoroCrustEdge::VoroCrustEdge(Vertex const& v1, 
                             Vertex const& v2, 
                             std::size_t const index_) : vertex1(v1), 
                                                         vertex2(v2), 
                                                         faces(), 
                                                         index(index_),
                                                         isSharp(false){}

bool VoroCrustEdge::checkIfEqual(Vertex const& v1, Vertex const& v2){
    if (v1->index == vertex1->index && v2->index == vertex2->index)
        return true;

    if (v2->index == vertex1->index && v1->index == vertex2->index)
        return true;

    return false;
}

void VoroCrustEdge::addFace(Face new_face){
    if(faces.size() == 2){
        std::cout << "ERROR : can't have more than 2 faces sharing an edge" << std::endl;
        exit(1);
    }

    faces.push_back(new_face);
}


double VoroCrustEdge::calcDihedralAngle(){
    if(faces.size() != 2){
        std::cout << "ERROR: can't calculate dihedral angle if number of faces " << std::endl;
        exit(1);
    }

    Vector3D const& b0 = vertex2->vertex - vertex1->vertex;
    Vector3D b1, b2;
    for(auto& edge : faces[0]->edges){
        if(edge->index == index) continue;

        if(edge->vertex1->index == vertex1->index){
            b1 = edge->vertex2->vertex - edge->vertex1->vertex;
            break;        
        }

        if(edge->vertex2->index == vertex1->index){
            b1 = edge->vertex1->vertex - edge->vertex2->vertex;
            break;        
        }
    }

    for(auto& edge : faces[1]->edges){
        if(edge->index == index) continue;

        if(edge->vertex1->index == vertex1->index){
            b2 = edge->vertex2->vertex - edge->vertex1->vertex;
            break;        
        }

        if(edge->vertex2->index == vertex1->index){
            b2 = edge->vertex1->vertex - edge->vertex2->vertex;
            break;        
        }
    }

    Vector3D const& b0_X_b1 = CrossProduct(b0, b1);
    Vector3D const& b0_X_b2 = CrossProduct(b0, b2);
    double const cos_phi = ScalarProd(b0_X_b1, b0_X_b2) / (abs(b0_X_b1) * abs(b0_X_b2));

    return std::acos(cos_phi);
}

std::string VoroCrustEdge::repr() {
    std::ostringstream s;

    s << "\n\t\t\tvertex " << vertex1->index << ": " << vertex1->repr() << ", vertex " << vertex2->index << " : " << vertex2->repr();

    return s.str();
}