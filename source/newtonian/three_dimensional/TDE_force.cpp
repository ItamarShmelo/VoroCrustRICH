#include "TDE_force.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/numeric/odeint.hpp>
using state_type = std::array<double, 4>;

namespace
{
    class PaczynskiOrbit
    {
    private:
        double M_, Rg_;

    public:
        PaczynskiOrbit(double M) : M_(M), Rg_(0)
        {
            Rg_ = 4.21 * M / 1e6;
        }

        void operator()(const state_type &x, state_type &dxdt, const double /* t */)
        {
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = -x[0] * M_ / (r * (r - Rg_) * (r - Rg_));
            dxdt[3] = -x[1] * M_ / (r * (r - Rg_) * (r - Rg_));
        }
    };

    state_type GetTrueAnomaly(double t, double M, double Rp)
    {
        double Rg = 4.21 * M / 1e6;
        double vp = std::sqrt(2 * M / (Rp - Rg));
        typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;
        PaczynskiOrbit orbit(M);
        state_type x0;
        x0[0] = Rp;
        x0[1] = 0;
        x0[2] = 0;
        x0[3] = vp;
        boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled<error_stepper_type>(1.0e-11, 1.0e-8), orbit,
                                                   x0, 0.0, t, t * 1e-5);
        return x0;
    }
}
void TDEGravity::operator()(const Tessellation3D &tess, const vector<ComputationalCell3D> &cells,
                            const vector<Conserved3D> &fluxes, const double time, vector<Vector3D> &acc) const
{
    // Calc self gravity
    selfgravity_(tess, cells, fluxes, time, acc);

    // Calc the force on the CM
    Vector3D Acm, Rcm;
    double Rg = 4.21 * Mbh_ / 1e6;
    if (tide_on_)
    {
        double Rt = R_ * std::pow(Mbh_ / M_, 0.333333333);
        double Rp = Rt / beta_;
        state_type x0 = GetTrueAnomaly(-time, Mbh_, Rp);
        double r = std::sqrt(x0[0] * x0[0] + x0[1] * x0[1]);
        Acm = Mbh_ * Vector3D(x0[0] / (r * (r - Rg) * (r - Rg)), x0[1] / (r * (r - Rg) * (r - Rg)), 0);
        Rcm = Vector3D(x0[0], x0[1], 0);
    }
    std::pair<Vector3D, Vector3D> box = tess.GetBoxCoordinates();
    double mindensity = std::max(1e-20, 1e-5 * M_ / ((box.second.x - box.first.x) * (box.second.z - box.first.z) * (box.second.y - box.first.y)));
    // Calc the tidal force
    size_t N = acc.size();
    double smooth = 60;
    for (size_t i = 0; i < N; ++i)
    {
        Vector3D const &point = tess.GetCellCM(i);
        Vector3D full_point = point + Rcm;
        double r_i = std::max(abs(full_point), 10.0);
        if (r_i > smooth)
            acc[i] += -(Mbh_ / (r_i * (r_i - Rg) * (r_i - Rg))) * full_point + Acm;
        else
        {
            double h = smooth;
            acc[i] += -(Mbh_ / (h * (h - Rg) * (h - Rg))) * full_point + Acm;
        }
        if (cells[i].density < mindensity || cells[i].tracers[1] < 0.1)
            acc[i] = Vector3D(0, 0, 0);
    }
}