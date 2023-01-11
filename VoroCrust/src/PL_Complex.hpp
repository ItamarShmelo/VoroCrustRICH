#ifndef PL_Complex_HPP
#define PL_Complex_HPP 1

#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include "VoroCrustVertex.hpp"
#include "VoroCrustFace.hpp"

class PL_Complex // Piecewise Linear Complex
{
    public:

        std::vector<std::shared_ptr<VoroCrustVertex>> vertices;
        std::vector<std::shared_ptr<VoroCrustFace>> faces;

        PL_Complex(std::vector<Vector3D> const& vertices);

        void addFace(std::vector<unsigned int> const& indices);

        std::string repr() const;

};

#endif /* PL_Complex_HPP */