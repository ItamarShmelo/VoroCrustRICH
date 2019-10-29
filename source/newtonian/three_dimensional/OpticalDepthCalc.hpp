#ifndef OPTDEPTHCALC_HPP
#define OPTDEPTHCALC_HPP 1

#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "computational_cell.hpp"

class OpticalDepthCalc
{
private:
	const double opening_;
	Tessellation3D const* tproc_;
	OpticalDepthCalc(const OpticalDepthCalc&/*other*/) :opening_(0), tproc_(0) {};
	OpticalDepthCalc& operator=(const OpticalDepthCalc /*other*/) { return *this; }
public:
	OpticalDepthCalc(double opening = 0.25, Tessellation3D const* tproc = 0);

	// returns dz * Sigma, so for time need to multiply by kappa and divide by c
	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		std::vector<double>& res) const;
	~OpticalDepthCalc() {}
};

#endif //OPTDEPTHCALC_HPP

