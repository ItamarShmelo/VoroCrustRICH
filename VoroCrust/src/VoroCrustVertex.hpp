#ifndef VOROCRUST_VERTEX_HPP
#define VOROCRUST_VERTEX_HPP 1

#include <vector>
#include <sstream>
#include <memory>
#include "../../source/3D/GeometryCommon/Vector3D.hpp"
#include "VoroCrustFace.hpp"
class VoroCrustFace;

class VoroCrustVertex
{
    public:
        Vector3D vertex;
        std::vector<std::shared_ptr<VoroCrustFace>> faces;

        VoroCrustVertex(Vector3D const& vertex_);

        void addFace(std::shared_ptr<VoroCrustFace> new_face);
        
        std::string repr() const;
};
#endif /* VOROCRUST_VERTEX_HPP */