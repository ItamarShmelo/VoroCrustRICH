#include "Diffusion.hpp"
#include <boost/math/special_functions.hpp>

namespace CG
{
     // LP flux limiter taken from "EQUATIONS AND ALGORITHMS FOR MIXED-FRAME FLUX-LIMITED DIFFUSION RADIATION HYDRODYNAMICS"
    double CalcSingleFluxLimiter(Vector3D const& grad, double const D, double const cell_value)
    {
        double const R = std::max(3 * abs(grad) * D / (cell_value * CG::speed_of_light), 1e20 * std::numeric_limits<double>::min());
        if(R < 1e-2) //series expansion
            return 1 - R * R / 15 + 2 * R * R * R * R /315;
        return 3 * (1.0 / std::tanh(R) - 1.0 / R) / R;
    }

    double FleckFactor(double const dt, double const beta, double const sigma_a)
    {
        return 1.0 / (1 + beta * dt * sigma_a * CG::speed_of_light);
    }

    double FleckFactorCompton(double const dt, double const beta, double const sigma_a, double const sigma_s, double const Erad, double const Cv)
    {
        return 1.0 / (1 + beta * dt * sigma_a * CG::speed_of_light + 16 * sigma_s * CG::boltzmann_constant * Erad / (CG::electron_mass * CG::speed_of_light * Cv));
    }
}

void Diffusion::BuildMatrix(Tessellation3D const& tess, mat& A, size_t_mat& A_indeces, std::vector<ComputationalCell3D> const& cells,
    double const dt, std::vector<double>& b, std::vector<double>& x0, double const current_time) const
{
    size_t const Nlocal = tess.GetPointNo();
    std::vector<ComputationalCell3D> cells_cgs(cells);
    for(size_t i = 0; i < Nlocal; ++i)
    {
        cells_cgs[i].density *= mass_scale_ / (length_scale_ * length_scale_ * length_scale_);
        cells_cgs[i].Erad *= length_scale_ * length_scale_ / (time_scale_ * time_scale_);
        cells_cgs[i].velocity *= length_scale_ / time_scale_ ;
    }
#ifdef RICH_MPI
	ComputationalCell3D cdummy;
	MPI_exchange_data(tess, cells_cgs, true, &cdummy);	
#endif
    b.resize(Nlocal, 0);
    x0.resize(Nlocal, 0);
    D.resize(Nlocal);
    fleck_factor.resize(Nlocal);
    sigma_planck.resize(Nlocal);
    sigma_s.resize(Nlocal);
    std::vector<size_t> neighbors;
    face_vec faces;
    std::vector<size_t> zero_indeces;
    size_t const Nzero = zero_cells_.size();
    for(size_t i = 0; i < Nzero; ++i)
        zero_indeces.push_back(binary_index_find(ComputationalCell3D::stickerNames, zero_cells_[i]));
    double const zero_value = 1e-10;
    std::vector<double> new_Er(Nlocal, 0);
    for(size_t i = 0; i < Nlocal; ++i)
    {
        double const volume = tess.GetVolume(i) * length_scale_ * length_scale_ * length_scale_;
        bool set_to_zero = false;
        for(size_t j = 0; j < Nzero; ++j)
            if(cells_cgs[i].stickers[zero_indeces[j]])
                set_to_zero = true;
        double const Er = cells_cgs[i].Erad * cells_cgs[i].density * (set_to_zero ? zero_value : 1);
        new_Er[i] = cells_cgs[i].Erad * cells_cgs[i].density;

        D[i] = D_coefficient_calcualtor.CalcDiffusionCoefficient(cells_cgs[i]);
        double const T = cells_cgs[i].temperature;
        sigma_planck[i] = D_coefficient_calcualtor.CalcPlanckOpacity(cells_cgs[i]);
        sigma_s[i] = D_coefficient_calcualtor.CalcScatteringOpacity(cells_cgs[i]);
        double const Cv = mass_scale_ * eos_.dT2cv(cells[i].density, T, cells[i].tracers, ComputationalCell3D::tracerNames) / (time_scale_ * time_scale_ * length_scale_);
        double const beta = 4 * CG::radiation_constant * T * T * T / Cv;
        fleck_factor[i] = compton_on_ ? FleckFactorCompton(dt * time_scale_, beta, sigma_planck[i], sigma_s[i], Er, Cv) : FleckFactor(dt * time_scale_, beta, sigma_planck[i]);
        b[i] = volume * Er;
        x0[i] = (Er + 0.5 * std::min(fleck_factor[i] * dt * sigma_planck[i] * CG::speed_of_light * time_scale_, 1.0) * (CG::radiation_constant * T * T * T * T - Er));
        b[i] += volume * fleck_factor[i] * dt * CG::speed_of_light * sigma_planck[i] * T * T * T * T * CG::radiation_constant * time_scale_;
    }
#ifdef RICH_MPI
    MPI_exchange_data2(tess, D, true);
#endif
    size_t max_neigh = 0;
    // Find maximum number of neighbors and allocate data
    for(size_t i = 0; i < Nlocal; ++i)
        max_neigh = std::max(max_neigh, tess.GetNeighbors(i).size());
    ++max_neigh;
    A.clear();
    A.resize(Nlocal);
    A_indeces.clear();
    A_indeces.resize(Nlocal);
    R2.clear();
    R2.resize(Nlocal, 0);
    cell_flux_limiter.clear();
    cell_flux_limiter.resize(Nlocal, 0);

    // Build the matrix
    for(size_t i = 0; i < Nlocal; ++i)
    {
        A_indeces[i].push_back(i);
        double const volume = tess.GetVolume(i) * length_scale_ * length_scale_ * length_scale_;
        double const T = cells_cgs[i].temperature;
        A[i].push_back(volume * (1 + fleck_factor[i] * dt * CG::speed_of_light * sigma_planck[i] * time_scale_));
        if(compton_on_)
            A[i][0] += volume * sigma_s[i] *  fleck_factor[i] * dt * time_scale_ * 4 * CG::boltzmann_constant * (std::pow(new_Er[i] / CG::radiation_constant, 0.25) - T) / (CG::electron_mass * CG::speed_of_light);
    }
    Vector3D dummy_v;
    std::vector<Vector3D> gradE(Nlocal);
    for(size_t i = 0; i < Nlocal; ++i)
    {
        double const volume = tess.GetVolume(i) * length_scale_ * length_scale_ * length_scale_;
        faces = tess.GetCellFaces(i);
        tess.GetNeighbors(i, neighbors);
        size_t const Nneigh = neighbors.size();
        Vector3D const CM = tess.GetCellCM(i);
        Vector3D const point = tess.GetMeshPoint(i);
        gradE[i] = Vector3D(0, 0, 0) ;
        double const Dcell = D[i];
        double const Er = cells_cgs[i].Erad * cells_cgs[i].density;
        bool self_zero = false;
        for(size_t k = 0; k < Nzero; ++k)
            if(cells_cgs[i].stickers[zero_indeces[k]])
                self_zero = true;
        for(size_t j = 0; j < Nneigh; ++j)
        {
            // Here we assume no flux to outside cells, this needs to be changed to a general boundary condition
            size_t const neighbor_j = neighbors[j];
            Vector3D r_ij = point - tess.GetMeshPoint(neighbor_j);
            double const r_ij_size = abs(r_ij);
            r_ij *= 1.0 / r_ij_size;
            double Er_j = 0;
            if(!tess.IsPointOutsideBox(neighbor_j))
            {
                bool set_to_zero = false;
                for(size_t k = 0; k < Nzero; ++k)
                    if(cells_cgs[neighbor_j].stickers[zero_indeces[k]])
                        set_to_zero = true;
                Er_j = cells_cgs[neighbor_j].Erad * cells_cgs[neighbor_j].density * (set_to_zero ? zero_value : 1);
                if(i < neighbor_j)
                {
                    Vector3D const cm_ij = CM - tess.GetCellCM(neighbor_j);
                    Vector3D const grad_E = cm_ij * (1.0 / (length_scale_ * ScalarProd(cm_ij, cm_ij)));
                    
                    double mid_D = 0.5 * (D[neighbor_j] + Dcell);
                    double const flux_limiter = flux_limiter_ ? CalcSingleFluxLimiter(grad_E * (Er - Er_j), mid_D, 0.5 * (Er + Er_j)) : 1;
                    mid_D *= flux_limiter;
                    double const flux = ((self_zero || set_to_zero) ? tess.GetArea(faces[j]) * dt * CG::speed_of_light * 0.5 : ScalarProd(grad_E, r_ij) * tess.GetArea(faces[j]) * dt * mid_D) * length_scale_ * length_scale_ * time_scale_; 
                    if(neighbor_j < Nlocal)
                    {
                        A[i][0] += flux;
                        A[i].push_back(-flux);
                        A_indeces[i].push_back(neighbor_j);
                        A[neighbor_j].push_back(-flux);
                        A_indeces[neighbor_j].push_back(i);
                        A[neighbor_j][0] += flux;
                    }                  
                    else
                    {
                        A[i][0] += flux;
                        A[i].push_back(-flux);
                        A_indeces[i].push_back(neighbor_j);
                    }
                }
            }
            else
            {
                if(i < neighbor_j)
                    boundary_calc_.SetBoundaryValues(tess, i, neighbor_j, dt * time_scale_, cells_cgs, tess.GetArea(faces[j]) * length_scale_ * length_scale_, A[i][0], b[i], faces[j]);
                boundary_calc_.GetOutSideValues(tess, cells_cgs, i, neighbor_j, new_Er, Er_j, dummy_v);
            }
            gradE[i] += r_ij * (tess.GetArea(faces[j]) * 0.5 * (Er + Er_j) * length_scale_ * length_scale_);
        }
    }
    for(size_t i = 0; i < Nlocal; ++i)
    {
        double const volume = tess.GetVolume(i) * length_scale_ * length_scale_* length_scale_;
        gradE[i] *= -1.0 / volume;
        faces = tess.GetCellFaces(i);
        tess.GetNeighbors(i, neighbors);
        size_t const Nneigh = neighbors.size();
        Vector3D const point = tess.GetMeshPoint(i);
        double const Dcell = D[i];
        double const Er = cells_cgs[i].Erad * cells_cgs[i].density;     
        double const flux_limiter = flux_limiter_ ? CalcSingleFluxLimiter(gradE[i], Dcell, Er) : 1;
        cell_flux_limiter[i] = flux_limiter;
        for(size_t j = 0; j < Nneigh; ++j)
        {
            size_t const neighbor_j = neighbors[j];
            if(!tess.IsPointOutsideBox(neighbor_j))
            {
                Vector3D r_ij = point - tess.GetMeshPoint(neighbor_j);
                double const r_ij_size = abs(r_ij);
                r_ij *= 1.0 / r_ij_size;
                double const momentum_relativity_term = -0.5 * fleck_factor[i] * dt * flux_limiter * tess.GetArea(faces[j]) * (2 * 3 * sigma_planck[i] * Dcell / CG::speed_of_light - 1) * length_scale_ * length_scale_ * time_scale_
                    * ScalarProd(cells_cgs[i].velocity, r_ij) / 3;
                A[i][0] += momentum_relativity_term;
                auto it = std::find(A_indeces[i].begin(), A_indeces[i].end(), neighbor_j);
                if(it == A_indeces[i].end())
                    throw UniversalError("Key not equal in diffusion");
                size_t const neigh_counter = static_cast<size_t>(it - A_indeces[i].begin());
                if(A_indeces[i][neigh_counter] != neighbor_j)
                    throw UniversalError("Key not equal value in diffusion");
                A[i][neigh_counter] += momentum_relativity_term;
            }
            else
                boundary_calc_.SetMomentumTermBoundary(tess, i, neighbor_j, dt * time_scale_, cells_cgs[i],
                    tess.GetArea(faces[j]) * length_scale_ * length_scale_, A[i][0], b[i], faces[j], fleck_factor[i],
                    flux_limiter, Dcell, sigma_planck[i]);
        }
        R2[i] = flux_limiter_ ? flux_limiter / 3 + boost::math::pow<2>(flux_limiter * abs(gradE[i]) * Dcell / (CG::speed_of_light * Er)) : 1.0 / 3.0;
        A[i][0] -= volume * fleck_factor[i] * dt * 0.5 * (3 - R2[i]) * sigma_planck[i] * ScalarProd(cells_cgs[i].velocity, cells_cgs[i].velocity) * time_scale_ / CG::speed_of_light;
    }
    for(size_t i = 0; i < Nlocal; ++i)
    {
        A[i].resize(max_neigh, 0);
        A_indeces[i].resize(max_neigh, max_size_t);
    }
}

void Diffusion::PostCG(Tessellation3D const& tess, std::vector<Conserved3D>& extensives, double const dt, std::vector<ComputationalCell3D>& cells,
        std::vector<double>const& CG_result)const
{
    Vector3D dummy_v;
    std::vector<size_t> neighbors;
    face_vec faces;
    size_t const N = tess.GetPointNo();
    bool const entropy = !(std::find(ComputationalCell3D::tracerNames.begin(), ComputationalCell3D::tracerNames.end(), std::string("Entropy")) ==
		ComputationalCell3D::tracerNames.end());
    size_t const entropy_index = static_cast<size_t>(std::find(ComputationalCell3D::tracerNames.begin(),
        ComputationalCell3D::tracerNames.end(), std::string("Entropy")) - ComputationalCell3D::tracerNames.begin());
    std::vector<size_t> zero_indeces;
    size_t const Nzero = zero_cells_.size();
    for(size_t i = 0; i < Nzero; ++i)
        zero_indeces.push_back(binary_index_find(ComputationalCell3D::stickerNames, zero_cells_[i]));
  
    for(size_t i = 0; i < N; ++i)
    {
        double const volume = tess.GetVolume(i) * length_scale_ * length_scale_* length_scale_;
        extensives[i].Erad = CG_result[i] * volume * time_scale_ * time_scale_ / (length_scale_ * length_scale_ * mass_scale_);
        double const T = cells[i].temperature;
        double dE = fleck_factor[i] * CG::speed_of_light * dt * sigma_planck[i] * (CG_result[i] - T * T * T * T * CG::radiation_constant
            -0.5 * (3 - R2[i]) * ScalarProd(cells[i].velocity, cells[i].velocity) * CG_result[i] * length_scale_ * length_scale_ / (CG::speed_of_light * CG::speed_of_light * time_scale_ * time_scale_)) * volume * time_scale_;
        if(compton_on_)
        {
            double const old_Tr = std::pow(cells[i].Erad * cells[i].density * mass_scale_ / (CG::radiation_constant * time_scale_ * time_scale_ * length_scale_), 0.25);
            dE += time_scale_ * volume * sigma_s[i] * fleck_factor[i] * dt * 4 * CG::boltzmann_constant * (old_Tr - T) / (CG::electron_mass * CG::speed_of_light);
        }
        dE *= time_scale_ * time_scale_ / (length_scale_ * length_scale_ * mass_scale_);
        extensives[i].energy += dE;
        extensives[i].internal_energy += dE;
        cells[i].Erad = extensives[i].Erad / extensives[i].mass;
        if(extensives[i].internal_energy < 0 || !std::isfinite(extensives[i].internal_energy) || cells[i].Erad < 0)
        {
            UniversalError eo("Bad internal energy in Diffusion::PostCG");
            eo.addEntry("cell index", i);
            eo.addEntry("energy", extensives[i].internal_energy);
            eo.addEntry("CG_result", CG_result[i]);
            eo.addEntry("T", T);
            eo.addEntry("Density", cells[i].density);
            eo.addEntry("ID", cells[i].ID);
            throw eo;
        }

        tess.GetNeighbors(i, neighbors);
        size_t const Nneigh = neighbors.size();
        faces = tess.GetCellFaces(i);
        Vector3D const CM = tess.GetCellCM(i);
        Vector3D const point = tess.GetMeshPoint(i);
        Vector3D gradE(0, 0, 0);
        double const Dcell = D[i];
        for(size_t j = 0; j < Nneigh; ++j)
        {
            size_t const neighbor_j = neighbors[j];
            Vector3D r_ij = point - tess.GetMeshPoint(neighbor_j);
            double const r_ij_size = abs(r_ij);
            r_ij *= 1.0 / r_ij_size;
            double Er_j = 0;
            if(tess.IsPointOutsideBox(neighbor_j))
                boundary_calc_.GetOutSideValues(tess, cells, i, neighbor_j, CG_result, Er_j, dummy_v);
            else
                Er_j = CG_result[neighbor_j];
            gradE += (0.5 * tess.GetArea(faces[j]) * (Er_j + CG_result[i])) * r_ij * length_scale_ * length_scale_;
            double const momentum_term = (0.5 * fleck_factor[i] * dt * cell_flux_limiter[i] * tess.GetArea(faces[j]) * ScalarProd(cells[i].velocity, r_ij) * (Er_j + CG_result[i]) / 3) * (time_scale_ * time_scale_ * length_scale_ / mass_scale_);
            double const relativity_term = -0.5 * momentum_term * 2 * 3 * sigma_planck[i] * Dcell / CG::speed_of_light;
            extensives[i].energy += momentum_term + relativity_term;
            extensives[i].internal_energy += relativity_term;
        }
        if(hydro_on_)
        {
            extensives[i].momentum += (cell_flux_limiter[i] * dt * time_scale_ / 3) * gradE * (time_scale_ / (length_scale_ * mass_scale_));
            extensives[i].energy = extensives[i].internal_energy +  ScalarProd(extensives[i].momentum, extensives[i].momentum) / (2 * extensives[i].mass);
        }

        cells[i].internal_energy = extensives[i].internal_energy / extensives[i].mass;
        cells[i].temperature = eos_.de2T(cells[i].density, cells[i].internal_energy, cells[i].tracers, ComputationalCell3D::tracerNames);
        cells[i].pressure = eos_.de2p(cells[i].density, cells[i].internal_energy, cells[i].tracers, ComputationalCell3D::tracerNames);
        cells[i].velocity = extensives[i].momentum / extensives[i].mass;
        
        if(entropy)
        {
            cells[i].tracers[entropy_index] = eos_.dp2s(cells[i].density, cells[i].pressure, cells[i].tracers, ComputationalCell3D::tracerNames);
            extensives[i].tracers[entropy_index] = cells[i].tracers[entropy_index] * extensives[i].mass;
        }

    }
}

void DiffusionSideBoundary::SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt, 
        std::vector<ComputationalCell3D> const& /*cells*/, double const Area, double& A, double &b, size_t const /*face_index*/)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        A += 0.5 * CG::speed_of_light * dt * Area;
        b += 2 * Area * dt * CG::stefan_boltzman * T_ * T_ * T_ * T_;
    }
}

void DiffusionSideBoundary::SetMomentumTermBoundary(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt,
        ComputationalCell3D const& cell, double const Area, double& A, double &b, size_t const face_index, double const fleck_factor,
        double const flux_limiter, double const D, double const sigma_planck)const
{
    double const R = tess.GetWidth(index);
    Vector3D r_ij = tess.GetMeshPoint(index) - tess.GetMeshPoint(outside_point);
    double const r_ij_size = abs(r_ij);
    r_ij *= 1.0 / r_ij_size;
    double const momentum_relativity_term = -0.5 * fleck_factor * dt * flux_limiter * Area * 
        (2 * 3 * sigma_planck * D / CG::speed_of_light - 1) * ScalarProd(cell.velocity, r_ij) / 3;
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        A += momentum_relativity_term;
        b -= momentum_relativity_term * CG::radiation_constant * T_ * T_ * T_ * T_;
    }
    else
        A += 2 * momentum_relativity_term;
}

void DiffusionSideBoundary::GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
    std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
        E_outside = CG::radiation_constant * T_ * T_ * T_ * T_;
    else
        E_outside = new_E[index];
    v_outside = cells[index].velocity;
}

void DiffusionClosedBox::SetMomentumTermBoundary(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt,
        ComputationalCell3D const& cell, double const Area, double& A, 
        double& /*b*/, size_t const /*face_index*/, double const fleck_factor, double const flux_limiter, 
        double const D, double const sigma_planck)const
{
     Vector3D r_ij = tess.GetMeshPoint(index) - tess.GetMeshPoint(outside_point);
    double const r_ij_size = abs(r_ij);
    r_ij *= 1.0 / r_ij_size;
    double const momentum_relativity_term = -0.5 * fleck_factor * dt * flux_limiter * Area * 
        (2 * 3 * sigma_planck * D / CG::speed_of_light - 1) * ScalarProd(cell.velocity, r_ij) / 3;
    A += 2 * momentum_relativity_term;
}

void DiffusionClosedBox::SetBoundaryValues(Tessellation3D const& /*tess*/, size_t const /*index*/, size_t const /*outside_point*/, double const /*dt*/, 
        std::vector<ComputationalCell3D> const& /*cells*/, double const /*Area*/, double& /*A*/, double& /*b*/, size_t const /*face_index*/)const
{}

void DiffusionClosedBox::GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
    std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const
{
    E_outside = new_E[index];
    Vector3D normal = normalize(tess.GetMeshPoint(outside_point) - tess.GetMeshPoint(index));
    v_outside = cells[index].velocity;
    v_outside -= 2 * normal * ScalarProd(normal, v_outside);
}

double PowerLawOpacity::CalcDiffusionCoefficient(ComputationalCell3D const& cell) const
{
    return D0_ * std::pow(cell.density, alpha_) * std::pow(cell.temperature, beta_);
}

double PowerLawOpacity::CalcPlanckOpacity(ComputationalCell3D const& cell) const
{
    return CG::speed_of_light / (3 * CalcDiffusionCoefficient(cell));
}

void DiffusionXInflowBoundary::SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt, 
        std::vector<ComputationalCell3D> const& cells, double const Area, double& A, double &b, size_t const face_index)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        double const Er_j = left_state_.Erad * left_state_.density;
        double const Er = cells[index].Erad * cells[index].density;
        double const dx = tess.GetCellCM(index).x - tess.GetBoxCoordinates().first.x;
        Vector3D const grad_E = Vector3D(1, 0, 0) * (1.0 / (2 * dx));
        Vector3D const r_ij = tess.GetMeshPoint(index) - tess.GetMeshPoint(outside_point);
        double mid_D = 0.5 * (D_calc_.CalcDiffusionCoefficient(cells[index]) + D_calc_.CalcDiffusionCoefficient(left_state_));
        double const flux_limiter = CalcSingleFluxLimiter(grad_E * (Er - Er_j), mid_D, 0.5 * (Er + Er_j));
        mid_D *= flux_limiter;
        double const flux = ScalarProd(grad_E, r_ij * (tess.GetArea(face_index) / abs(r_ij))) * dt * mid_D; 
        A += flux;
        b += flux * Er_j;
    }
    else
    {
        if(tess.GetMeshPoint(index).x < (tess.GetMeshPoint(outside_point).x - R * 1e-4))
        {
            double const Er_j = right_state_.Erad * right_state_.density;
            double const Er = cells[index].Erad * cells[index].density;
            double const dx = tess.GetBoxCoordinates().second.x - tess.GetCellCM(index).x;
            Vector3D const grad_E = Vector3D(-1, 0, 0) * (1.0 / (2 * dx));
            Vector3D const r_ij = tess.GetMeshPoint(index) - tess.GetMeshPoint(outside_point);
            double mid_D = 0.5 * (D_calc_.CalcDiffusionCoefficient(cells[index]) + D_calc_.CalcDiffusionCoefficient(right_state_));
            double const flux_limiter = CalcSingleFluxLimiter(grad_E * (Er - Er_j), mid_D, 0.5 * (Er + Er_j));
            mid_D *= flux_limiter;
            double const flux = ScalarProd(grad_E, r_ij * (tess.GetArea(face_index) / abs(r_ij))) * dt * mid_D; 
            A += flux;
            b += flux * Er_j;
        }
    }
}

void DiffusionXInflowBoundary::SetMomentumTermBoundary(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt,
        ComputationalCell3D const& cell, double const Area, double& A, double &b, size_t const face_index, double const fleck_factor,
        double const flux_limiter, double const D, double const sigma_planck)const
{
    double const R = tess.GetWidth(index);
    Vector3D r_ij = tess.GetMeshPoint(index) - tess.GetMeshPoint(outside_point);
    double const r_ij_size = abs(r_ij);
    r_ij *= 1.0 / r_ij_size;
    double const momentum_relativity_term = -0.5 * fleck_factor * dt * flux_limiter * Area * 
        (2 * 3 * sigma_planck * D / CG::speed_of_light - 1) * ScalarProd(cell.velocity, r_ij) / 3;
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        A += momentum_relativity_term;
        b -= momentum_relativity_term * left_state_.Erad * left_state_.density;
    }
    else
    {
        if(tess.GetMeshPoint(index).x < (tess.GetMeshPoint(outside_point).x - R * 1e-4))
        {
            A += momentum_relativity_term;
            b -= momentum_relativity_term * right_state_.Erad * right_state_.density;
        }
        else
            A += 2 * momentum_relativity_term;
    }
}

void DiffusionXInflowBoundary::GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
    std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        E_outside = left_state_.Erad * left_state_.density;
        v_outside = left_state_.velocity;
    }
    else
    {
        if(tess.GetMeshPoint(index).x < (tess.GetMeshPoint(outside_point).x - R * 1e-4))
        {
            E_outside = right_state_.Erad * right_state_.density;
            v_outside = right_state_.velocity;
        }
        else
        {
            E_outside = new_E[index];
            Vector3D normal = normalize(tess.GetMeshPoint(outside_point) - tess.GetMeshPoint(index));
            v_outside = cells[index].velocity;
            v_outside -= 2 * normal * ScalarProd(normal, v_outside);
        }
    }
}

void DiffusionOpenBoundary::SetBoundaryValues(Tessellation3D const& /*tess*/, size_t const /*index*/, size_t const /*outside_point*/, double const dt,
    std::vector<ComputationalCell3D> const& /*cells*/, double const Area, double& A, double& /*b*/, size_t const /*face_index*/)const
{
    A += Area * dt * 0.5 * CG::speed_of_light;
}

void DiffusionOpenBoundary::GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
    std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const
{
    E_outside = new_E[index] * 1e-20;
    v_outside = cells[index].velocity;
}

void DiffusionOpenBoundary::SetMomentumTermBoundary(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt,
    ComputationalCell3D const& cell, double const Area, double& A, double &b, size_t const face_index, 
    double const fleck_factor, double const flux_limiter, double const D, double const sigma_planck)const
{
    Vector3D r_ij = tess.GetMeshPoint(index) - tess.GetMeshPoint(outside_point);
    double const r_ij_size = abs(r_ij);
    r_ij *= 1.0 / r_ij_size;
    double const momentum_relativity_term = -0.5 * fleck_factor * dt * flux_limiter * Area * (2 * 3 * 
        sigma_planck * D / CG::speed_of_light - 1) * ScalarProd(cell.velocity, r_ij) / 3;
    A += momentum_relativity_term;
}