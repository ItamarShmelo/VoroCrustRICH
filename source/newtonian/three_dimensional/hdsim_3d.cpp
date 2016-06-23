#include <cassert>
#include "hdsim_3d.hpp"

HDSim3D::ProgressTracker::ProgressTracker(void) :
	time_(0), cycle_(0) {}

void HDSim3D::ProgressTracker::update(double dt)
{
	++cycle_;
	time_ += dt;
}

double HDSim3D::ProgressTracker::getTime(void) const
{
	return time_;
}

double HDSim3D::ProgressTracker::getCycle(void) const
{
	return cycle_;
}

HDSim3D::HDSim3D(Tessellation3D& tess,
	const vector<ComputationalCell3D>& cells,
	const EquationOfState& eos,
	const PointMotion3D& pm,
	const TimeStepFunction3D& tsc,
	const FluxCalculator3D& fc,
	const CellUpdater3D& cu,
	const ExtensiveUpdater3D& eu,
	const TracerStickerNames tsn) :
	tess_(tess), eos_(eos), cells_(cells),extensive_(),pm_(pm), tsc_(tsc), fc_(fc), cu_(cu),eu_(eu),tsn_(tsn),pt_()
{
	assert(tess.GetPointNo() == cells.size());
	size_t N = tess.GetPointNo();
	extensive_.resize(N);
	for (size_t i = 0; i < N; ++i)
		PrimitiveToConserved(cells_[i], tess.GetVolume(i), extensive_[i], eos_, tsn_);
}

namespace
{
	void CalcFaceVelocities(Tessellation3D const& tess, vector<Vector3D> const& point_vel, vector<Vector3D> &res)
	{
		size_t N = tess.GetTotalFacesNumber();
		res.resize(N);
		for (size_t i = 0; i < N; ++i)
			if (tess.BoundaryFace(i))
				res[i] = Vector3D();
			else
				tess.CalcFaceVelocity(i, point_vel[tess.GetFaceNeighbors(i).first], point_vel[tess.GetFaceNeighbors(i).second]);
	}

	void UpdateTessellation(Tessellation3D &tess, vector<Vector3D> &point_vel,double dt)
	{
		vector<Vector3D> points = tess.GetMeshPoints();
		size_t N = tess.GetPointNo();
		points.resize(N);
		for (size_t i = 0; i < N; ++i)
			points[i] += point_vel[i] * dt;
		tess.Build(points);
	}
}

void HDSim3D::timeAdvance(void)
{
	vector<Vector3D> point_vel,face_vel;
	pm_(tess_, cells_, pt_.getTime(), tsn_, point_vel);
	CalcFaceVelocities(tess_, point_vel, face_vel);
	const double dt = tsc_(tess_, cells_, eos_,face_vel,pt_.getTime(),tsn_);
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel, tsn_);
	CalcFaceVelocities(tess_, point_vel, face_vel);
	vector<Conserved3D> fluxes;
	fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), dt, tsn_);
	eu_(fluxes, tess_, dt, cells_, extensive_, pt_.getTime(), tsn_);
	UpdateTessellation(tess_, point_vel, dt);
	cu_(cells_, eos_, tess_, extensive_, tsn_);
	pt_.update(dt);
}

const Tessellation3D& HDSim3D::getTesselation(void) const
{
	return tess_;
}

const vector<ComputationalCell3D>& HDSim3D::getCells(void) const
{
	return cells_;
}

double HDSim3D::GetTime(void)const
{
	return pt_.getTime();
}

TracerStickerNames HDSim3D::GetTracerStickerNames(void)const
{
	return tsn_;
}

size_t HDSim3D::GetCycle(void)const
{
	return pt_.getCycle();
}