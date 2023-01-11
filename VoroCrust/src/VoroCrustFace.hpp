#ifndef VOROCRUST_FACE_HPP
#define VOROCRUST_FACE_HPP 1

#include <vector>
#include <sstream>
#include <memory>
#include "../../source/3D/GeometryCommon/Vector3D.hpp"
#include "VoroCrustVertex.hpp"
class VoroCrustVertex;

class VoroCrustFace
{
    public:
        std::vector<std::shared_ptr<VoroCrustVertex>> vertices;
        
        VoroCrustFace(std::vector<std::shared_ptr<VoroCrustVertex>> const& vertices_);

        std::string repr() const;

};

#endif /* VOROCRUST_FACE_HPP */