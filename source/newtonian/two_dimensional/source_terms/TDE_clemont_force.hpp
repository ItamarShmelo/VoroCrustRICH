#include "ConservativeForce.hpp"
#include "TDE_clemont.hpp"

class TDE_clemont_force : public SourceTerm
{
public:
	TDE_clemont acc_;

	explicit TDE_clemont_force(double const Mbh, Vector2D const CM, Vector2D const Vcm, double const time, bool const full = false);

	vector<Extensive> operator()(const Tessellation& tess, const PhysicalGeometry& pg, 
        const CacheData& cd, const vector<ComputationalCell>& cells, const vector<Extensive>& fluxes,
        const vector<Vector2D>& point_velocities, const double t) const override;

private:
    bool const full_;
    ConservativeForce force_;
};