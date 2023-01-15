#ifndef VOROCRUST_VERTEX_HPP
#define VOROCRUST_VERTEX_HPP 1


#include <vector>
#include <string>
#include "../../source/3D/GeometryCommon/Vector3D.hpp"

#include "VoroCrustUsing.hpp"
#include "VoroCrustFace.hpp"
#include "VoroCrustEdge.hpp"


class VoroCrustVertex
{
    public:
        Vector3D vertex;
        std::vector<Face> faces;
        std::vector<Edge> edges;

        std::size_t index;
        bool isSharp;

        VoroCrustVertex(Vector3D const& vertex_, std::size_t const index_);

        ~VoroCrustVertex() = default;
        
        void addFace(Face new_face);

        void addEdge(Edge const& new_edge) {edges.push_back(new_edge);}

        std::string repr() const;
};
#endif /* VOROCRUST_VERTEX_HPP */