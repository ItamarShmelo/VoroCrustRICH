#ifndef PL_Complex_HPP
#define PL_Complex_HPP 1

#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include "VoroCrustVertex.hpp"
#include "VoroCrustEdge.hpp"
#include "VoroCrustFace.hpp"

//! \brief a Crease is a chain of sharp edges, ends in a sharp corner or forms a cycle. 
using Crease = std::vector<Edge>;
//! \brief a Surface patch is the connected componnent containing no sharp features. A Surface patch is enveloped by Creases.
using SurfacePatch = std::vector<Face>;

/*! \brief Piecewise Linear Complex 
    \details Holds the mesh and does the preprocessing of the VoroCrust algorithm.
*/
class PL_Complex 
{
    public:

        //! \brief the PLC mesh vertices.
        std::vector<Vertex> vertices;
        //! \brief the PLC mesh edges.
        std::vector<Edge> edges;
        //! \brief the PLC mesh faces.
        std::vector<Face> faces;

        //! \brief holds the sharp edges of the PLC.
        std::vector<Edge> sharp_edges;
        //! \brief holds the sharp corners (sharp vertices) of the PLC.
        std::vector<Vertex> sharp_corners;

        //! \brief holds the Creases of the PLC.
        std::vector<Crease> creases;
        //! \brief holds the Surface Patches of the PLC.
        std::vector<SurfacePatch> patches;

        /*! \brief Constructor, consturcts a PLC from a given set of Vectors representing the vertices of the mesh.
            \param vertices a vector<Vector3D> of the initial vertices of the mesh.
        */
        PL_Complex(std::vector<Vector3D> const& vertices);
        
        ~PL_Complex() = default;

        /*! \brief adds an edge to the PLC starting at v1 and ending at v2.
            \param v1 first Vertex of Edge.
            \param v2 second Vertex of Edge.
        */
        Edge addEdge(Vertex const& v1, Vertex const& v2);
        
        /*! \brief adds a face to the PLC
            \param indices vector of indices of the `vertices` of the Face.
        */
        void addFace(std::vector<unsigned int> const& indices);

        /*! \brief Checks if all the vertices are assigned at least on face.*/
        bool checkAllVerticesAreOnFace();

        /*! \brief Checks if all defined faces are on planes */
        bool checkIfALLFacesAreFlat();
        
        /*! \brief check that all dihedral angles between faces are less than `\pi-sharpTheta` or above `\pi-flatTheta` 
            \param sharpTheta determins the sharp features.
            \param flatTheta constraint on the flatness of the non-sharp features.
        */
        void detectFeatures(double const sharpTheta, double const flatTheta);

        /*! \brief builds the Creases */
        void buildCreases();

        /*! \brief builds a Crease s.t. `edge` is in it.
            \param edge 
        */
        Crease createCrease(Edge const& edge);
        
        /*! \brief builds the Surface Patches */
        void buildSurfacePatches();
        
        /*! \brief create a Surface Patch s.t. face is in it 
            \param face
        */
        SurfacePatch createSurfacePatch(Face const& face);
        
        std::string repr() const;

};


#endif /* PL_Complex_HPP */