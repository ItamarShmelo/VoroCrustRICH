#ifndef VOROCRUST_VERTEX_HPP
#define VOROCRUST_VERTEX_HPP 1


#include <vector>
#include <string>
#include "../../../source/3D/GeometryCommon/Vector3D.hpp"

#include "VoroCrustUsing.hpp"
#include "VoroCrustFace.hpp"
#include "VoroCrustEdge.hpp"


class VoroCrustVertex
{
    public:
        //! \brief Vector3D representing the Vertex in 3D space
        Vector3D vertex; 
        //! \brief Faces incident to the Vertex.
        std::vector<Face> faces;
        std::vector<std::vector<Face>> divided_faces;
        //! \brief Edges incident to the Vertex.
        std::vector<Edge> edges;

        //! \brief index of Vertex in PLC vertices vector.
        std::size_t index;

        //! \brief true if Vertex is a sharp corner.
        bool isSharp;

        VoroCrustVertex(Vector3D const& vertex_, std::size_t const index_);

        ~VoroCrustVertex() = default;
        
        //! \brief adds a Face to Vertex's `faces` vector.
        void addFace(Face const& new_face);
        
        //! \brief adds a Face to Vertex's `edges` vector.
        void addEdge(Edge const& new_edge);

        std::string repr() const;
};
#endif /* VOROCRUST_VERTEX_HPP */