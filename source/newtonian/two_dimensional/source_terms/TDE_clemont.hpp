#include "ConservativeForce.hpp"
#include "../../..//misc/utils.hpp"

class TDE_clemont : public Acceleration
{
    public:
        mutable Vector2D CM_, Vcm_;
        mutable double time_, delta_, lambda_, omega_, V_;
        double const Mbh_;
        bool const full_;

    TDE_clemont(double const Mbh, double const time, Vector2D const CM, Vector2D const Vcm, bool const full = false);

    Vector2D operator()(const Tessellation &tess, const vector<ComputationalCell> &cells,
                        const vector<Extensive> &fluxes, const double time, const int point) const;
};