#ifndef VOROCRUST_FACE_HPP
#define VOROCRUST_FACE_HPP 1

#include <vector>
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"

#include "VoroCrustUsing.hpp"
#include "VoroCrustVertex.hpp"

class VoroCrustFace
{
    public:
        //! \brief Vertices defining the face
        std::vector<Vertex> vertices;
        //! \brief Edges on the boundary of the Face
        std::vector<Edge> edges;
        //! \brief Neighboring Faces
        std::vector<Face> neighbors;
        
        //! \brief index in the PLC Faces vector
        std::size_t index;

        //! \brief the current calculated normal to the face.
        Vector3D current_normal;

        //! \brief true if Face was already assigned to a Patch
        bool isPatched;
        
        VoroCrustFace(std::vector<Vertex> const& vertices_, std::size_t const index_);
        
        ~VoroCrustFace() = default;

        //! \brief adds an Edge to the Face's Edges vector.
        void addEdge(Edge const& edge);

        /*! \brief calculate the normal to the surface defined by the vectors `vertices[2] - vertices[1]` X `vertices[1]-vertices[0]`*/
        Vector3D calcNormal();

        std::string repr() const;

};

#endif /* VOROCRUST_FACE_HPP */