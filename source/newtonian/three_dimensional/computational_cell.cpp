#include "computational_cell.hpp"

ComputationalCell3D::ComputationalCell3D(void):
  density(0), pressure(0),internal_energy(0),temperature(0),ID(0), velocity(), Erad(0), Erad_dt(0),
  	Erad_dt_dt(0), tracers(),stickers() {}

ComputationalCell3D::ComputationalCell3D(double density_i,
				     double pressure_i,double internal_energy_i,size_t ID_i,
				     const Vector3D& velocity_i):
  density(density_i), pressure(pressure_i),internal_energy(internal_energy_i),temperature(0),ID(ID_i),
  velocity(velocity_i), Erad(0), Erad_dt(0), Erad_dt_dt(0), tracers(),stickers() {}

ComputationalCell3D::ComputationalCell3D(double density_i,
				     double pressure_i, double internal_energy_i,size_t ID_i,
				     const Vector3D& velocity_i,
				     const std::array<double,MAX_TRACERS>& tracers_i,
					 const std::array<bool,MAX_STICKERS>& stickers_i):
  density(density_i), pressure(pressure_i),internal_energy(internal_energy_i),temperature(0),ID(ID_i),
  velocity(velocity_i), Erad(0), Erad_dt(0), Erad_dt_dt(0), tracers(tracers_i),stickers(stickers_i) {}

ComputationalCell3D::ComputationalCell3D(const ComputationalCell3D& other):
density(other.density),
pressure(other.pressure),
internal_energy(other.internal_energy),
temperature(other.temperature),
ID(other.ID),
velocity(other.velocity),
Erad(other.Erad),
Erad_dt(other.Erad_dt),
Erad_dt_dt(other.Erad_dt_dt),
tracers(other.tracers),
stickers(other.stickers) {}


ComputationalCell3D& ComputationalCell3D::operator=(ComputationalCell3D const& other)
{
	density = other.density;
	pressure = other.pressure;
	internal_energy = other.internal_energy;
	temperature = other.temperature;
	velocity = other.velocity;
	Erad = other.Erad;
	Erad_dt = other.Erad_dt;
	Erad_dt_dt = other.Erad_dt_dt;
	tracers = other.tracers;
	stickers = other.stickers;
	ID = other.ID;
	return *this;
}

ComputationalCell3D& ComputationalCell3D::operator+=(ComputationalCell3D const& other)
{
	this->density += other.density;
	this->pressure += other.pressure;
	this->internal_energy += other.internal_energy;
	this->velocity += other.velocity;
	this->temperature += other.temperature;
	this->Erad += other.Erad;
	this->Erad_dt += other.Erad_dt;
	this->Erad_dt_dt += other.Erad_dt_dt;
	//assert(this->tracers.size() == other.tracers.size());
	//size_t N = this->tracers.size();
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		this->tracers[j] += other.tracers[j];
	return *this;
}

ComputationalCell3D& ComputationalCell3D::operator-=(ComputationalCell3D const& other)
{
	this->density -= other.density;
	this->pressure -= other.pressure;
	this->internal_energy -= other.internal_energy;
	this->velocity -= other.velocity;
	this->temperature -= other.temperature;
	this->Erad -= other.Erad;
	this->Erad_dt += other.Erad_dt;
	this->Erad_dt_dt += other.Erad_dt_dt;	
	//assert(this->tracers.size() == other.tracers.size());
	//size_t N = this->tracers.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		this->tracers[j] -= other.tracers[j];
	return *this;
}

ComputationalCell3D& ComputationalCell3D::operator*=(double s)
{
	this->density *= s;
	this->pressure *= s;
	this->internal_energy *= s;
	this->velocity *= s;
	this->temperature *= s;
	this->Erad *= s;
	this->Erad_dt *= s;
	this->Erad_dt_dt *= s;
	//size_t N = this->tracers.size();
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		this->tracers[j] *= s;
	return *this;
}

vector<string> ComputationalCell3D::tracerNames;
vector<string> ComputationalCell3D::stickerNames;

#ifdef RICH_MPI
size_t ComputationalCell3D::getChunkSize(void) const
{
	return 11 + tracers.size() + stickers.size();
}

vector<double> ComputationalCell3D::serialize(void) const
{
	vector<double> res(getChunkSize());
	res.at(0) = density;
	res.at(1) = pressure;
	res.at(2) = velocity.x;
	res.at(3) = velocity.y;
	res.at(4) = velocity.z;
	res.at(5) = internal_energy;
	res.at(6) = temperature;
	res.at(7) = static_cast<double>(ID);
	res.at(8) = Erad;
	res.at(9) = Erad_dt;
	res.at(10) = Erad_dt_dt;
	size_t counter = 11;
	//size_t N = tracers.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		res[j + counter] = tracers[j];
	//size_t N2 = stickers.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t j = 0; j < MAX_STICKERS; ++j)
		res[j + counter + MAX_TRACERS] = stickers[j] ? 1 : 0;
	return res;
}

void ComputationalCell3D::unserialize
(const vector<double>& data)
{
	assert(data.size() == getChunkSize());
	density = data.at(0);
	pressure = data.at(1);
	velocity.x = data.at(2);
	velocity.y = data.at(3);
	velocity.z = data.at(4);
	internal_energy = data.at(5);
	temperature = data.at(6);
	ID = static_cast<size_t>(std::llround(data.at(7)));
	Erad = data.at(8);
	Erad_dt = data.at(9);
	Erad_dt_dt = data.at(10);
	size_t counter = 11;
	//size_t N = tracers.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		tracers[j] = data.at(counter + j);
	//size_t N2 = stickers.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t i = 0; i < MAX_STICKERS; ++i)
		stickers[i] = data.at(counter + MAX_TRACERS + i)>0.5;
}

size_t Slope3D::getChunkSize(void) const
{
	return xderivative.getChunkSize() * 3;
}

vector<double> Slope3D::serialize(void) const
{
	vector<double> res(getChunkSize());
	vector<double> temp(xderivative.serialize());
	res.reserve(temp.size() * 3);
	std::copy(temp.begin(), temp.end(), res.begin());
	temp = yderivative.serialize();
	std::copy(temp.begin(), temp.end(), res.begin() + xderivative.getChunkSize());
	temp = zderivative.serialize();
	std::copy(temp.begin(), temp.end(), res.begin() + xderivative.getChunkSize() + yderivative.getChunkSize());
	return res;
}

void Slope3D::unserialize(const vector<double>& data)
{
	size_t size = xderivative.getChunkSize();
	xderivative.unserialize(vector<double>(data.begin(), data.begin() + size));
	yderivative.unserialize(vector<double>(data.begin() + size, data.begin() + 2*size));
	zderivative.unserialize(vector<double>(data.begin() + 2*size, data.end()));
}

#endif // RICH_MPI


void ComputationalCellAddMult(ComputationalCell3D &res, ComputationalCell3D const& other, double scalar)
{
	res.density += other.density*scalar;
	res.pressure += other.pressure*scalar;
	res.internal_energy += other.internal_energy*scalar;
	res.velocity += other.velocity*scalar;
	res.temperature += other.temperature*scalar;
	res.Erad += other.Erad*scalar;
	res.Erad_dt += other.Erad_dt*scalar;
	res.Erad_dt_dt += other.Erad_dt_dt*scalar;
	//assert(res.tracers.size() == other.tracers.size());
	//size_t N = res.tracers.size();
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		res.tracers[j] += other.tracers[j] * scalar;
}

ComputationalCell3D operator+(ComputationalCell3D const& p1, ComputationalCell3D const& p2)
{
	ComputationalCell3D res(p1);
	res += p2;
	return res;
}

ComputationalCell3D operator-(ComputationalCell3D const& p1, ComputationalCell3D const& p2)
{
	ComputationalCell3D res(p1);
	res -= p2;
	return res;
}

ComputationalCell3D operator/(ComputationalCell3D const& p, double s)
{
	ComputationalCell3D res(p);
	double const s_1 = 1.0 / s;
	res.density *= s_1;
	res.pressure *= s_1;
	res.internal_energy *= s_1;
	res.temperature *= s_1;
	res.Erad *= s_1;
	res.Erad_dt *= s_1;
	res.Erad_dt_dt *= s_1;
	//size_t N = res.tracers.size();
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		res.tracers[j] *= s_1;
	res.velocity = res.velocity * s_1;
	return res;
}

ComputationalCell3D operator*(ComputationalCell3D const& p, double s)
{
	ComputationalCell3D res(p);
	res.density *= s;
	res.pressure *= s;
	res.internal_energy *= s;
	res.temperature *= s;
	res.Erad *= s;
	res.Erad_dt *= s;
	res.Erad_dt_dt *= s;
	//size_t N = res.tracers.size();
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		res.tracers[j] *= s;
	res.velocity = res.velocity * s;
	return res;
}

ComputationalCell3D operator*(double s, ComputationalCell3D const& p)
{
	return p*s;
}

void ReplaceComputationalCell(ComputationalCell3D & cell, ComputationalCell3D const& other)
{
	cell.density = other.density;
	cell.pressure = other.pressure;
	cell.internal_energy = other.internal_energy;
	cell.ID = other.ID;
	cell.velocity = other.velocity;
	cell.temperature = other.temperature;
	cell.Erad = other.Erad;
	cell.Erad_dt = other.Erad_dt;
	cell.Erad_dt_dt = other.Erad_dt_dt;
	//size_t N = other.tracers.size();
	//cell.tracers.resize(N);
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		cell.tracers[j] = other.tracers[j];
	//N = other.stickers.size();
	//cell.stickers.resize(N);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t i = 0; i < MAX_STICKERS; ++i)
		cell.stickers[i] = other.stickers[i];
}



Slope3D::Slope3D(void) : xderivative(ComputationalCell3D()), yderivative(ComputationalCell3D()), zderivative(ComputationalCell3D()) {}

Slope3D::Slope3D(ComputationalCell3D const & x, ComputationalCell3D const & y,ComputationalCell3D const & z) : xderivative(x), yderivative(y),zderivative(z)
{}
