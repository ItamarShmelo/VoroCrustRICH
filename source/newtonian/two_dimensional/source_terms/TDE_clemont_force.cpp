#include "TDE_clemont_force.hpp"

TDE_clemont_force::TDE_clemont_force(double const Mbh, Vector2D const CM, Vector2D const Vcm, double const time, bool const full)
:acc_(Mbh, time, CM, Vcm, full), full_(full), force_(acc_) {}

std::vector<Extensive> TDE_clemont_force::operator()(const Tessellation& tess, const PhysicalGeometry& pg, 
        const CacheData& cd, const vector<ComputationalCell>& cells, const vector<Extensive>& fluxes,
        const vector<Vector2D>& point_velocities, const double t) const
{
    std::vector<Extensive> res = force_(tess, pg, cd, cells, fluxes, point_velocities, t);
    size_t v_index = binary_index_find(ComputationalCell::tracerNames, std::string("v_parallel"));
    if(v_index >= ComputationalCell::tracerNames.size())
        throw UniversalError("Bad v_index");
    size_t const N = tess.GetPointNo();
    double const R = abs(acc_.CM_);
    double const R_3 = 1.0 / (R * R * R);
    for(size_t i = 0; i < N; ++i)
    {
        double const V = tess.GetMeshPoint(i).x * acc_.omega_ + cells[i].tracers[v_index];
        if(full_)
        {
            Vector2D v_hat = acc_.Vcm_ * (1.0 / abs(acc_.Vcm_));
            Vector2D Acm = -acc_.Mbh_ * acc_.CM_ * R_3;
            Vector2D location = tess.GetMeshPoint(i);
            Vector2D full_location = acc_.CM_ + Vector2D(-location.x * v_hat.y, location.x * v_hat.x);
            double const full_r = std::max(30.0, std::sqrt(ScalarProd(full_location, full_location) + location.y * location.y));
            Vector2D tides = -acc_.Mbh_ * full_location / (full_r * full_r * full_r) - Acm;
            res[i].tracers[v_index] = ScalarProd(tides, v_hat);
        }
        else
            res[i].tracers[v_index] = 3 * acc_.Mbh_ * tess.GetMeshPoint(i).x * std::cos(acc_.delta_) * std::sin(acc_.delta_) * R_3;
        double const factor = std::min(1.0, std::pow(cells[i].density / 1e-13, 0.75));
        res[i].tracers[v_index] += cells[i].velocity.x * acc_.omega_ - V * acc_.lambda_;
        res[i].tracers[v_index] *= factor * cd.volumes[i] * cells[i].density;
    }
    return res;
}