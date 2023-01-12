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
        std::vector<std::shared_ptr<VoroCrustFace>> neighbors;
        
        std::size_t index;
        Vector3D current_normal;
        
        VoroCrustFace(std::vector<std::shared_ptr<VoroCrustVertex>> const& vertices_, std::size_t const index_);
        
        ~VoroCrustFace() = default;

        void addEdge(std::shared_ptr<VoroCrustEdge> edge);

        /*! \brief calculate the normal to the surface defined by the vectors `vertices[2] - vertices[1]`, `vertices[1]-vertices[0]`*/
        Vector3D calcNormal();

        std::string repr() const;

};

#endif /* VOROCRUST_FACE_HPP */