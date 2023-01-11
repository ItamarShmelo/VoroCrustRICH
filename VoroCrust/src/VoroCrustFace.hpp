#ifndef VOROCRUST_FACE_HPP
#define VOROCRUST_FACE_HPP 1

#include <vector>
#include <sstream>
#include <memory>
#include "../../source/3D/GeometryCommon/Vector3D.hpp"
#include "VoroCrustVertex.hpp"
class VoroCrustVertex;
class VoroCrustEdge;

class VoroCrustFace
{
    public:
        std::vector<std::shared_ptr<VoroCrustVertex>> vertices;
        std::vector<std::shared_ptr<VoroCrustEdge>> edges;
        std::size_t index;
        
        VoroCrustFace(std::vector<std::shared_ptr<VoroCrustVertex>> const& vertices_, std::size_t const index_);
        
        ~VoroCrustFace() = default;

        void addEdge(std::shared_ptr<VoroCrustEdge> edge) {edges.push_back(edge);};

        std::string repr() const;

};

#endif /* VOROCRUST_FACE_HPP */