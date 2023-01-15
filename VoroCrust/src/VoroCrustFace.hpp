#ifndef VOROCRUST_FACE_HPP
#define VOROCRUST_FACE_HPP 1

#include <vector>
#include "../../source/3D/GeometryCommon/Vector3D.hpp"

#include "VoroCrustUsing.hpp"
#include "VoroCrustVertex.hpp"

class VoroCrustFace
{
    public:
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
        std::vector<Face> neighbors;
        
        std::size_t index;
        Vector3D current_normal;
        
        VoroCrustFace(std::vector<Vertex> const& vertices_, std::size_t const index_);
        
        ~VoroCrustFace() = default;

        void addEdge(Edge edge);

        /*! \brief calculate the normal to the surface defined by the vectors `vertices[2] - vertices[1]`, `vertices[1]-vertices[0]`*/
        Vector3D calcNormal();

        std::string repr() const;

};

#endif /* VOROCRUST_FACE_HPP */