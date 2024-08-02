/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(SOURCEFIELDS_HPP_)
#error "This file should only be included through SourceFields.hpp"
#endif

#ifndef SOURCEFIELDS_IMPL_HPP_
#define SOURCEFIELDS_IMPL_HPP_
#define PI (3.14159265358979323846264338327950288E0)

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
AMREX_GPU_DEVICE emtensor_t<data_t> SourceFields<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, d1, h_UU, chris_ULL);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi   = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    out.rho += V_of_phi;
    out.S   += -3.0 * V_of_phi;
    FOR (i, j)
    {
        out.Sij[i][j] += -vars.h[i][j] * V_of_phi / vars.chi;
    }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
SourceFields<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL)
{
    // Useful quantity Vt
    data_t Vt = -vars.Pi * vars.Pi;
    FOR (i, j)
    {
        Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR (i, j)
    {
        out.Sij[i][j] =
            -0.5 * vars.h[i][j] * Vt / vars.chi + d1.phi[i] * d1.phi[j];
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = - n^a T_ai
    FOR (i)
    {
        out.Si[i] = -d1.phi[i] * vars.Pi;
    }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi * vars.Pi + 0.5 * Vt;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
SourceFields<potential_t>::add_sources_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec,
    std::default_random_engine &random_generator) const
{
    // first get the non potential part of the rhs
    // this may seem a bit long winded, but it makes the function
    // work for more multiple fields

    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, d1, d2, advec);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi   = 0.0;
    my_potential.compute_potential(V_of_phi, dVdphi, vars);

    // adjust RHS for the potential term
    total_rhs.Pi += -vars.lapse * dVdphi;

    //  RHS for Rlin
    total_rhs.Rlin  = 0.0;

    // Let's see as a test what happens if:
    // std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
    // total_rhs.phi += vars.Rlin*distribution(generator);
    // double temp=distribution(random_generator);
    // amrex::Print() << "Random draw : " << temp << std::endl;
    // total_rhs.phi += pow(-vars.K/3,1.5)/(2*PI)*temp;
    total_rhs.phi += pow(-vars.K/3,1.5)/(2*PI)*distribution(random_generator);

}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
SourceFields<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec)
{
    using namespace TensorAlgebra;

    const auto h_UU  = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs.phi = vars.lapse * vars.Pi + advec.phi;
    rhs.Pi  = vars.lapse * vars.K * vars.Pi + advec.Pi;

    FOR (i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Pi += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi[i] +
                                vars.chi * vars.lapse * d2.phi[i][j] +
                                vars.chi * d1.lapse[i] * d1.phi[j]);
        FOR (k)
        {
            rhs.Pi += -vars.chi * vars.lapse * h_UU[i][j] * chris.ULL[k][i][j] *
                      d1.phi[k];
        }
    }
}

#endif /* SOURCEFIELDS_IMPL_HPP_ */
