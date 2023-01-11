#ifndef PL_Complex_HPP
#define PL_Complex_HPP 1

#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include "VoroCrustVertex.hpp"
#include "VoroCrustEdge.hpp"
#include "VoroCrustFace.hpp"

class PL_Complex // Piecewise Linear Complex
{
    public:

        std::vector<std::shared_ptr<VoroCrustVertex>> vertices;
        std::vector<std::shared_ptr<VoroCrustEdge>> edges;
        std::vector<std::shared_ptr<VoroCrustFace>> faces;

        PL_Complex(std::vector<Vector3D> const& vertices);
        ~PL_Complex() = default;

        std::shared_ptr<VoroCrustEdge> addEdge(std::shared_ptr<VoroCrustVertex> const& v1, std::shared_ptr<VoroCrustVertex> const& v2);
        
        void addFace(std::vector<unsigned int> const& indices);

        /*! \brief Checks if all the vertices are assigned at least on face.*/
        bool checkAllVerticesAreOnFace();

        /*! \brief Checks if all defined faces are on planes */
        bool checkIfALLFacesAreFlat();

        std::string repr() const;

};

#endif /* PL_Complex_HPP */