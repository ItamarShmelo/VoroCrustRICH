#include "Vector3D.hpp"

void Split(vector<Vector3D> const & vIn, vector<double> & vX, vector<double> & vY, vector<double> & vZ)
{
	vX.resize(vIn.size());
	vY.resize(vIn.size());
	vZ.resize(vIn.size());

	for (std::size_t ii = 0; ii < vIn.size(); ++ii)
	{
		vX[ii] = vIn[ii].x;
		vY[ii] = vIn[ii].y;
		vZ[ii] = vIn[ii].z;
	}
	return;
}


