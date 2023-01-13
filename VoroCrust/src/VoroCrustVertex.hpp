#ifndef VOROCRUST_VERTEX_HPP
#define VOROCRUST_VERTEX_HPP 1

#include <vector>
#include <sstream>
#include <memory>
#include "../../source/3D/GeometryCommon/Vector3D.hpp"
#include "VoroCrustFace.hpp"
#include "VoroCrustEdge.hpp"
class VoroCrustFace;
class VoroCrustEdge;

class VoroCrustVertex
{
    public:
        Vector3D vertex;
        std::vector<std::shared_ptr<VoroCrustFace>> faces;
        std::vector<std::shared_ptr<VoroCrustEdge>> edges;

        std::size_t index;
        bool isSharp;

        VoroCrustVertex(Vector3D const& vertex_, std::size_t const index_);

        ~VoroCrustVertex() = default;
        
        void addFace(std::shared_ptr<VoroCrustFace> new_face);

        void addEdge(std::shared_ptr<VoroCrustEdge> const& new_edge) {edges.push_back(new_edge);}

        std::string repr() const;
};
#endif /* VOROCRUST_VERTEX_HPP */