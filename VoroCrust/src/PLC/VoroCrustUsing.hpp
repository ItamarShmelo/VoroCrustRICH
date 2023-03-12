#ifndef VOROCRUST_USING
#define VOROCRUST_USING

#include <memory>

class VoroCrustVertex;
class VoroCrustEdge;
class VoroCrustFace;

using Vertex = std::shared_ptr<VoroCrustVertex>;
using Edge = std::shared_ptr<VoroCrustEdge>;
using Face = std::shared_ptr<VoroCrustFace>;

using VertexWeakPtr = std::weak_ptr<VoroCrustVertex>;
using EdgeWeakPtr = std::weak_ptr<VoroCrustVertex>;
using FaceWeakPtr = std::weak_ptr<VoroCrustVertex>;

#endif /* VOROCRUST_USING */