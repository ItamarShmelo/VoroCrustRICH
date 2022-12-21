#include "TDE_clemont.hpp"
#include <limits>

TDE_clemont::TDE_clemont(double const Mbh, double const time, Vector2D const CM, Vector2D const Vcm, bool const full)
: Mbh_(Mbh), time_(time), CM_(CM), Vcm_(Vcm), full_(full)
{
    delta_ = std::acos(ScalarProd(CM_, Vcm_) / (abs(CM_) * abs(Vcm_)));
    double const R = abs(CM_);
    double const s_delta = std::sin(delta_);
    double const c_delta = std::cos(delta_);
    omega_ = -Mbh_ * s_delta / (R * R * abs(Vcm));
    lambda_ = -Mbh_ * c_delta / (R * R * abs(Vcm));
}

Vector2D TDE_clemont::operator()(const Tessellation &tess, const vector<ComputationalCell> &cells,
                                  const vector<Extensive> &/*fluxes*/, const double time, const int point) const
{
    // Advance in time if needed the orbital paramaters
    if(std::abs(time_ - time) > time * std::numeric_limits<double>::epsilon() * 10)
    {
        double r = abs(CM_);
        double const dt = time - time_;
        Vcm_.x -= 0.5 * Mbh_ * CM_.x * dt / (r * r * r);
        Vcm_.y -= 0.5 * Mbh_ * CM_.y * dt / (r * r * r);
        CM_.x += Vcm_.x * dt;
        CM_.y += Vcm_.y * dt;
        r = abs(CM_);
        Vcm_.x -= 0.5 * Mbh_ * CM_.x * dt / (r * r * r);
        Vcm_.y -= 0.5 * Mbh_ * CM_.y * dt / (r * r * r);
        //delta_ = std::asin((-CM_.x * Vcm_.y + CM_.y * Vcm_.x) / (r * abs(Vcm_)));;
        delta_ = -std::acos(-std::abs(ScalarProd(CM_, Vcm_)) / (r * abs(Vcm_)));
        double const delta2 = -M_PI + std::acos(-std::abs(ScalarProd(CM_, Vcm_)) / (r * abs(Vcm_)));
        omega_ = -Mbh_ * std::sin(delta_) / (r * r * abs(Vcm_));
        lambda_ = -Mbh_ * std::cos(delta_) / (r * r * abs(Vcm_));
        time_= time;
        std::cout<<"Omega "<<omega_<<" delta "<<delta_<<" delta2 "<<delta2<<std::endl;
    }
    Vector2D res;
    Vector2D location = tess.GetMeshPoint(point);
    double const R = abs(CM_);
    double const s_delta = std::sin(delta_);
    double const rg = 4.21;
    double const V = location.x * omega_ + cells[point].tracers[1];
    double const R_3 = 1.0 / (R * (R - rg) * (R - rg));
    Vector2D v_hat = Vcm_ * (1.0 / abs(Vcm_));
    if(full_)
    {
        Vector2D Acm = -Mbh_ * CM_ * R_3;
        Vector2D full_location = CM_ + Vector2D(-location.x * v_hat.y, location.x * v_hat.x);
        double const full_r = std::max(30.0, std::sqrt(ScalarProd(full_location, full_location) + location.y * location.y));
        res = -Mbh_ * full_location / (full_r * (full_r - rg) * (full_r  - rg)) - Acm;
        res.x = -res.x * v_hat.y + res.y * v_hat.x + location.x * omega_ * omega_
             -2 * V * omega_;
        res.y = -Mbh_ * location.y / (full_r * full_r * full_r);
        double const factor = std::min(1.0, std::pow(cells[point].density / 1e-13, 0.75));
        res *= factor;
        // res.x = -Mbh_ * location.x * (1 - 3 * s_delta * s_delta) * R_3 + location.x * omega_ * omega_
        //     -2 * V * omega_;
       // res.y = -Mbh_ * location.y * R_3;    
    }
    else
    {
        res.x = -Mbh_ * location.x * (1 - 3 * s_delta * s_delta) * R_3 + location.x * omega_ * omega_
            -2 * V * omega_;
        res.y = -Mbh_ * location.y * R_3;
    }
    return res;
}