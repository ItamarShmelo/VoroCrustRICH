#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#endif
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/tessellation/RoundGrid.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/two_dimensional/ghost_point_generators/RigidWallGenerator.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/misc/simple_io.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/two_dimensional/condition_action_sequence_2.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/two_dimensional/interpolations/LinearGaussImproved.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"
#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#include "source/mpi/ConstNumberPerProc.hpp"
#endif

#define cylindar_run 1

int main(void)
{
    int rank = 0;
#ifdef RICH_MPI
	MPI_Init(NULL,NULL);
    int ws = 0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&ws);
#endif
    // Number of linear points
    size_t const Np = 64;
    // End time of run
    double const tend = 0.05;
    // Create box and outer boundary conditions
    Vector2D lower_left(0, 0), upper_right(1, 1);
    SquareBox outer_boundary(lower_left, upper_right);
    

#ifdef RICH_MPI
    //Create domain decomposition
    std::vector<Vector2D> cpu_points = RandSquare(ws, lower_left.x, upper_right.x, lower_left.y, upper_right.y);
    VoronoiMesh cpu_tess;
    cpu_tess.Initialise(cpu_points, outer_boundary);
    // Create initial random mesh
    std::vector<Vector2D> mesh_points = RandSquare(Np * Np, cpu_tess, lower_left, upper_right);
#else
    // Create initial random mesh
    std::vector<Vector2D> mesh_points = RandSquare(Np * Np, lower_left.x, upper_right.x, lower_left.y, upper_right.y);
#endif
    // Make mesh nice and round
#ifdef RICH_MPI
    mesh_points = RoundGridV(mesh_points, outer_boundary, cpu_tess);
#else
    mesh_points = RoundGridV(mesh_points, outer_boundary);
#endif
    // Create voronoi mesh
    VoronoiMesh tess;
#ifdef RICH_MPI
    tess.Initialise(mesh_points, cpu_tess, outer_boundary);
#else
    tess.Initialise(mesh_points, outer_boundary);
#endif
    // Choose geometry
#ifdef cylindar_run 
    CylindricalSymmetry geometry(Vector2D(0, 0), Vector2D(0, 1));
#else
    SlabSymmetry geometry;
#endif
    // Set how outer interfaces move
    StationaryBox interface_velocity_calc;
    // Set source/force
#ifdef cylindar_run 
    CylindricalComplementary source(geometry.getAxis());
#else
    ZeroForce source;
#endif
    // Set time step calculator
    SimpleCFL time_step(0.3);
    // Riemann solver
    Hllc rs;
    //Equation of state
    IdealGas eos(5.0 / 3.0);
    //Ghost point generator
    RigidWallGenerator ghost;
    // Spatial interpolation scheme
    LinearGaussImproved interpolation(eos, ghost);
    // Set flux calculator
    std::vector<std::pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> > first_order;
    std::vector<std::pair<const ConditionActionSequence::Condition*, const ConditionActionSequence2::Action2*> > second_order;
    RegularFlux2 reg_flux(rs);
    RigidWallFlux2 rigid_flux(rs);
    second_order.push_back(std::pair<const ConditionActionSequence::Condition*, const ConditionActionSequence2::Action2*>(new IsBoundaryEdge(), &rigid_flux));
    second_order.push_back(std::pair<const ConditionActionSequence::Condition*, const ConditionActionSequence2::Action2*>(new IsBulkEdge(), &reg_flux));
    ConditionActionSequence2 flux_calc(first_order, second_order, interpolation);
    // Set mesh point motion
    Lagrangian base_point_motion;
    RoundCells point_motion(base_point_motion, eos);
    //How to update the extensive variables
    SimpleExtensiveUpdater extensive_update;
    //How to update the primitive variables
    SimpleCellUpdater cell_update;
    //Set the initial hydro conditions
#ifdef cylindar_run 
    double const pressure_ratio = 1e4;
#else
    double const pressure_ratio = 1e3;
#endif
    std::vector<ComputationalCell> cells(tess.GetPointNo());
    for(size_t i = 0; i < tess.GetPointNo(); ++i)
    {
        cells[i].density = 1;
        cells[i].pressure = abs(tess.GetMeshPoint(i)) < 0.1 ? pressure_ratio : 1;
    }
    // Create main sim class
#ifdef RICH_MPI
    // Create load balance scheme
    ConstNumberPerProc cpu_move(outer_boundary);
    hdsim sim(cpu_tess, tess, outer_boundary, geometry, cells, eos, point_motion, interface_velocity_calc, source, time_step, flux_calc, extensive_update, 
        cell_update, std::pair<std::vector<std::string>, std::vector<std::string> > (), &cpu_move);
#else
    hdsim sim(tess, outer_boundary, geometry, cells, eos, point_motion, interface_velocity_calc, source, time_step, flux_calc, extensive_update, cell_update);
#endif
    // Main loop
    while(sim.getTime() < tend)
    {
        try
        {
            if(rank == 0)
                std::cout<<"Cycle "<<sim.getCycle()<<" time "<<sim.getTime()<<std::endl;
           sim.TimeAdvance2Heun();
        }
        catch(UniversalError const& eo)
        {
            reportError(eo);
            throw;
        }  
    }
#ifdef RICH_MPI
    write_snapshot_to_hdf5(sim, "sedov_rank" + std::to_string(rank) + ".h5");
    MPI_Finalize();
#else
    write_snapshot_to_hdf5(sim, "sedov.h5");
#endif
    return 0;
}

