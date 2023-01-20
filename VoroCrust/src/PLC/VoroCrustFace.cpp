#include "VoroCrustFace.hpp"
#include <iostream>
#include <sstream>

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
    Vector3D const& v = vertices[0]->vertex;

    double signed_area = 0;
    for (unsigned int i=1; i<vertices.size()-1; ++i){
        Vector3D const& v1 = vertices[i]->vertex - v;
        Vector3D const& v2 = vertices[i+1]->vertex - v;

        signed_area += abs(CrossProduct(v1, v2));
    }

    return 0.5*signed_area;
}

double VoroCrustFace::calcArea(){
    return std::abs(calcSignedArea());
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