#ifndef WRITE_VTU_PL_COMPLEX
#define WRITE_VTU_PL_COMPLEX

#include "PL_Complex.hpp"
#include <filesystem>
#include "../../source/3D/GeometryCommon/Vector3D.hpp"
#include <memory>

namespace vorocrust_vtk{

void write_vtu_PL_Complex(std::filesystem::path const& filename, PL_Complex const& plc);

} // namespace

#endif /* WRITE_VTU_PL_COMPLEX */