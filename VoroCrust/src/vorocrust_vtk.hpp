#ifndef WRITE_VTU_PL_COMPLEX
#define WRITE_VTU_PL_COMPLEX

#include "PLC/PL_Complex.hpp"
#include "VoroCrust_kd_tree/VoroCrust_kd_tree.hpp"
#include "trees.hpp"
#include <filesystem>
#include "../../source/3D/GeometryCommon/Vector3D.hpp"
#include <memory>

namespace vorocrust_vtk{

/*! \brief writes a PLC as vtu file.
    \param filename path of output file.
    \param plc  PL_Complex to be exported as vtu. 
*/
void write_vtu_PL_Complex(std::filesystem::path const& filename, PL_Complex const& plc);

void write_vtu_trees(std::filesystem::path const& filename, Trees const& trees);

void write_nearestNeighbor(std::filesystem::path const& filename, VoroCrust_KD_Tree const& tree, Vector3D const& query);

void write_kNearestNeighbors(std::filesystem::path const& filename, VoroCrust_KD_Tree const& tree, Vector3D const& query, int k);

void write_ballTree(std::filesystem::path const& filename, VoroCrust_KD_Tree_Ball const& b_tree);

} // namespace

#endif /* WRITE_VTU_PL_COMPLEX */