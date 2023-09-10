#include "source/3D/GeometryCommon/Voronoi3D.hpp"
#include "source/newtonian/three_dimensional/hdsim_3d.hpp"
#include "source/newtonian/three_dimensional/PCM3D.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/newtonian/common/TillotsonOrg.hpp"
#include "source/newtonian/three_dimensional/Hllc3D.hpp"
#include "source/misc/simple_io.hpp"
#include "source/misc/mesh_generator3D.hpp"
#include "source/newtonian/three_dimensional/Lagrangian3D.hpp"
#include "source/newtonian/three_dimensional/RoundCells3D.hpp"
#include "source/newtonian/three_dimensional/default_cell_updater.hpp"
#include "source/newtonian/three_dimensional/ConditionActionFlux1.hpp"
#include "source/newtonian/three_dimensional/ConditionExtensiveUpdater3D.hpp"
#include "source/newtonian/three_dimensional/CourantFriedrichsLewy.hpp"
#include "source/newtonian/three_dimensional/Ghost3D.hpp"
#include "source/3D/GeometryCommon/hdf_write.hpp"
#include <filesystem>
namespace fs = std::filesystem;
#include <sstream>
#include "source/VoroCrust/runVoroCrust.hpp"

int main(void)
{
    VoroCrustAlgorithm* alg = 0;
    std::vector<PL_Complex> zones_plcs;
    std::pair<std::vector<std::vector<Vector3D> >, std::vector<std::vector<Vector3D> >> points = runVoroCrust("/home/elads/RICH/source/VoroCrust/data/astroid", "/home/elads/source/VoroCrustDumps/astroid", alg, zones_plcs, true);  // volume_seeds, boundary_seeds
    std::cout<<"point size "<<points.first.size()<<" "<<points.second.size()<<" "<<points.second[0].size()<<std::endl;
    size_t const tot_size = points.second[0].size();
    double ll_x, ll_y, ll_z, ur_x, ur_y, ur_z;
    ll_x = ll_y = ll_z = std::numeric_limits<double>::max();
    ur_x = ur_y = ur_z = -std::numeric_limits<double>::max();
    
    for(size_t j = 0; j < points.second.size(); ++j)
    {
        std::cout<<"surface size "<<points.second[j].size()<<std::endl;
        for(size_t i=0; i<points.second[j].size(); ++i){
            Vector3D const& p = points.second[j][i];
            ll_x = std::min(p.x, ll_x);
            ll_y = std::min(p.y, ll_y);
            ll_z = std::min(p.z, ll_z);

            ur_x = std::max(p.x, ur_x);
            ur_y = std::max(p.y, ur_y);
            ur_z = std::max(p.z, ur_z);
        }
    }
    std::vector<Vector3D> all_points = points.second[0];
    double const off = 0.1;
    Vector3D ll(ll_x - off, ll_y - off, ll_z - off), ur(ur_x + off, ur_y + off, ur_z + off); // GIVE HERE VALUES
    std::vector<Vector3D> rand_points = RandRectangular(1e4, ll, ur);
    size_t const Nboundary_in = tot_size;
    size_t const Nboundary = all_points.size();
    std::cout<<"checking inside "<<zones_plcs.size()<<std::endl;
    for(size_t i = 0; i < rand_points.size(); ++i)
        if(alg->pointInSidePLC(zones_plcs[0], rand_points[i]))
            all_points.push_back(rand_points[i]);
    size_t const NInner = all_points.size();
    all_points.insert(all_points.end(), points.second[1].begin(), points.second[1].end());
    std::cout<<"checking outside"<<std::endl;
    for(size_t i = 0; i < rand_points.size(); ++i)
        if(alg->pointOutSidePLC(zones_plcs[0], rand_points[i]))
            all_points.push_back(rand_points[i]);
    size_t const Nsphere = all_points.size();
    std::cout<<Nboundary<<" "<<NInner<<" "<<Nsphere<<std::endl;
    rand_points = RandSphereR(1e4, ll, ur, 0, 0.05, Vector3D(0, 0.2, 0));
    std::cout<<"adding sphere"<<std::endl;
    all_points.insert(all_points.end(), rand_points.begin(), rand_points.end());
    size_t const N = all_points.size();
    for(size_t i = 0; i < N; ++i)
        all_points[i] *= 1e5;
    ll *= 1e5;
    ur *= 1e5;
	Voronoi3D tess(ll, ur);
    try{
        std::cout<<"starting voronoi with "<<all_points.size()<<" points"<<std::endl;
	    tess.Build(all_points);
    } catch (UniversalError const& eo){
        std::cout << eo.getErrorMessage() << std::endl;
        exit(1);
    }
  
    TillotsonOrg eos(0.5, 1.5, 2.67e11, 2.67e11, 2.7, 4.87e12, 4.72e10, 1.82e11, 5, 5);
    ComputationalCell3D basalt_cell, out_cell;
    basalt_cell.density = 2.7;
    basalt_cell.internal_energy = 1e8;
    basalt_cell.pressure = eos.de2p(basalt_cell.density, basalt_cell.internal_energy);
    out_cell.density = 1e-5;
    out_cell.internal_energy = 1e10;
    out_cell.pressure = eos.de2p(out_cell.density, out_cell.internal_energy);
	std::vector<ComputationalCell3D> cells(N);
    ComputationalCell3D::stickerNames.push_back("Inside");
    
    for(size_t i = 0; i < N; ++i)
    {  
        if(i < NInner || i >= Nsphere)
        {
            cells[i] = basalt_cell;
            if(i >= Nsphere)
                cells[i].velocity = Vector3D(0, 0, -2e5);
            else
                cells[i].stickers[0] = true;
        }
        else
            cells[i] = out_cell;
    }

	Hllc3D rs;
	RigidWallGenerator3D ghost;
	LinearGauss3D interp(eos, ghost, true, 0.2, 0.5, 0.7, false, std::vector<std::string>(), std::string(), false);

	Lagrangian3D bpm;
    RoundCells3D pm(bpm, eos, 0.5);
	DefaultCellUpdater cu;

	vector<pair<const ConditionActionFlux1::Condition3D *, const ConditionActionFlux1::Action3D *>> flux_vector;
	ConditionActionFlux1 fc(flux_vector, interp);

	vector<pair<const ConditionExtensiveUpdater3D::Condition3D *, const ConditionExtensiveUpdater3D::Action3D *>> eu_sequence;
	ConditionExtensiveUpdater3D eu(eu_sequence);

    ZeroForce3D force;
	CourantFriedrichsLewy tsf(0.25, 1, force);

    
	HDSim3D sim(tess, cells, eos, pm, tsf, fc, cu, eu, force, std::pair<std::vector<std::string>, std::vector<std::string>> (ComputationalCell3D::tracerNames, ComputationalCell3D::stickerNames), false
		, true);
	vector<DiagnosticAppendix3D *> appendices;
	WriteSnapshot3D(sim, "init.h5", appendices, true);

	return 0;
}
