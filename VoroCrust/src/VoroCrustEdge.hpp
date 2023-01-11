#ifndef VOROCRUST_EDGE
#define VOROCRUST_EDGE

#include <vector>
#include <memory>
#include <string>
#include "VoroCrustVertex.hpp"
#include "VoroCrustFace.hpp"

class VoroCrustVertex;
class VoroCrustFace;

class VoroCrustEdge
{

public:
    std::shared_ptr<VoroCrustVertex> vertex1, vertex2;
    std::vector<std::shared_ptr<VoroCrustFace>> faces;
    
    std::size_t index;

    VoroCrustEdge(std::shared_ptr<VoroCrustVertex> const& v1, std::shared_ptr<VoroCrustVertex> const& v2, std::size_t const index_);
    
    ~VoroCrustEdge() = default;

    bool checkIfEqual(std::shared_ptr<VoroCrustVertex> const& v1, std::shared_ptr<VoroCrustVertex> const& v2);

    void addFace(std::shared_ptr<VoroCrustFace> new_face);

    std::string repr();
};

#endif /* VOROCRUST_EDGE */ 
