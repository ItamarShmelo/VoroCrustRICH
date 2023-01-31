#include "VoroCrustFace.hpp"
#include <iostream>
#include <sstream>
#include <algorithm>

VoroCrustFace::VoroCrustFace(std::vector<Vertex> const &vertices_,
                             std::size_t const index_) : vertices(vertices_),
                                                         edges(),
                                                         neighbors(),
                                                         index(index_),
                                                         current_normal(),
                                                         isPatched(false) {}

void VoroCrustFace::addEdge(Edge const& edge)
{
    //! TODO: CHECK if edge is really on the boundary.
    edges.push_back(edge);
}

Vector3D VoroCrustFace::calcNormal()
{
    Vector3D const &v1 = vertices[1]->vertex - vertices[0]->vertex;
    Vector3D const &v2 = vertices[2]->vertex - vertices[1]->vertex;

    Vector3D const &cross_product = CrossProduct(v2, v1);
    current_normal = cross_product / abs(cross_product);

    return current_normal;
}

double VoroCrustFace::calcSignedArea(){
    /*calculates the signed area using formula in https://math.stackexchange.com/questions/3207981/how-do-you-calculate-the-area-of-a-2d-polygon-in-3d*/
    //! WARNING: This is not really a signed area !!!!
    Vector3D const& v = vertices[0]->vertex;

    Vector3D signed_area_vec({0, 0, 0});
    for (unsigned int i=1; i<vertices.size()-1; ++i){
        Vector3D const& v1 = vertices[i]->vertex - v;
        Vector3D const& v2 = vertices[i+1]->vertex - v;

        signed_area_vec += CrossProduct(v1, v2);
    }

    return 0.5*abs(signed_area_vec);
}

double VoroCrustFace::calcArea(){
    return std::abs(calcSignedArea());
}


void VoroCrustFace::flipOrientation(){
    std::reverse(vertices.begin(), vertices.end());
}

void VoroCrustFace::orientWithRespectTo(Face const& face){
    Vector3D const& n1 = this->calcNormal();
    Vector3D const& n2 = face->calcNormal();
    
    if(ScalarProd(n1, n2) < 0){
        std::cout << "flips orientation of Face " << index << ", because of Face " << face->index << "\n";
        flipOrientation();
    }
}

Vector3D VoroCrustFace::calculateCenteroid() const {
    Vector3D res(0.0, 0.0, 0.0);
    for(Vertex const& vertex : vertices){
        res += vertex->vertex;
    }

    return res / 3.0;
}

std::string VoroCrustFace::repr() const
{
    std::ostringstream s;
    s << "\tEdges: ";
    for (auto &edge : edges)
    {
        s << "\n\t\tedge " << edge->index << ": " << edge->repr() << ", ";
    }

    return s.str();
}