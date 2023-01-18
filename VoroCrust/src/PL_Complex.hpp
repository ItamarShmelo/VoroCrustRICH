#ifndef PL_Complex_HPP
#define PL_Complex_HPP 1

#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include "VoroCrustVertex.hpp"
#include "VoroCrustEdge.hpp"
#include "VoroCrustFace.hpp"

using Crease = std::vector<Edge>;
using SurfacePatch = std::vector<Face>;

class PL_Complex // Piecewise Linear Complex
{
    public:

        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
        std::vector<Face> faces;

        std::vector<Edge> sharp_edges;
        std::vector<Vertex> sharp_corners;

        std::vector<Crease> creases;
        std::vector<SurfacePatch> patches;

        PL_Complex(std::vector<Vector3D> const& vertices);
        ~PL_Complex() = default;

        Edge addEdge(Vertex const& v1, Vertex const& v2);
        
        void addFace(std::vector<unsigned int> const& indices);

        /*! \brief Checks if all the vertices are assigned at least on face.*/
        bool checkAllVerticesAreOnFace();

        /*! \brief Checks if all defined faces are on planes */
        bool checkIfALLFacesAreFlat();
        
        /*! \brief check that all dihedral angles between faces are less than `\pi-sharpTheta` or above `\pi-flatTheta` */
        void detectFeatures(double const sharpTheta, double const flatTheta);

        /*! \brief builds the creases */
        void buildCreases();

        /*! \brief builds a crease starting at `edge */
        Crease createCrease(Edge const& edge);

        std::string repr() const;

};


#endif /* PL_Complex_HPP */