#ifndef VOROCRUST_EDGE
#define VOROCRUST_EDGE

#include <vector>
#include <string>

#include "VoroCrustUsing.hpp"
#include "VoroCrustVertex.hpp"
#include "VoroCrustFace.hpp"

class VoroCrustEdge
{

public:
    Vertex vertex1, vertex2;
    std::vector<Face> faces;
    
    std::size_t index;
    bool isSharp;
    bool isCreased;

    VoroCrustEdge(Vertex const& v1, Vertex const& v2, std::size_t const index_);
    
    ~VoroCrustEdge() = default;

    bool checkIfEqual(Vertex const& v1, Vertex const& v2);

    void addFace(Face new_face);

    double calcDihedralAngle();

    std::string repr();
};

#endif /* VOROCRUST_EDGE */ 
