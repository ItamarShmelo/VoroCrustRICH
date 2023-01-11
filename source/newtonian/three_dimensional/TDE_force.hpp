#include "ConservativeForce3D.hpp"
#ifndef TDE_FORCE_HPP
#define TDE_FORCE_HPP 1

class TDEGravity : public Acceleration3D
{
private:
    Acceleration3D const &selfgravity_;
    const double Mbh_, M_, R_, beta_;
    const bool tide_on_;

public:
    TDEGravity(double Mbh, double M, double R, double beta, Acceleration3D const &sg, bool tide) 
        : selfgravity_(sg), Mbh_(Mbh), M_(M), R_(R), beta_(beta), tide_on_(tide) {}

    void operator()(const Tessellation3D &tess, const vector<ComputationalCell3D> &cells,
                    const vector<Conserved3D> &fluxes, const double time, vector<Vector3D> &acc) const;
};
#endif