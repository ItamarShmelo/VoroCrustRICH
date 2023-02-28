#ifndef UPDATE_BOX_HPP
#define UPDATE_BOX_HPP 1

#include "../../newtonian/three_dimensional/hdsim_3d.hpp"
#include "Voronoi3D.hpp"

void UpdateBox(Tessellation3D &tess, HDSim3D &sim, double const min_velocity, double const volume_fraction, ComputationalCell3D const& reference_cell
#ifdef RICH_MPI
    , Tessellation3D &tproc
#endif
);

#endif // UPDATE_BOX_HPP