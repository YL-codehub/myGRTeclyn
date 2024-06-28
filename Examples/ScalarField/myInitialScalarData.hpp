/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MYINITIALSCALARDATA_HPP_
#define MYINITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total no. components
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

template <class data_t> double Trilinear_Interpolation(Coordinates<data_t> coordinates, const double &dx1, const double &dx2, const int &n, const std::vector<double> &V) 
{
    // This assumes that dx2<dx1 and that the spacing of V is dx2 whereas coordinates come from dx1 spacing. Also assumes that n*dx2=N*dx1.
    double nx = std::fmod(coordinates.x/dx2-0.5,n);
    double ny = std::fmod(coordinates.y/dx2-0.5,n);
    double nz = std::fmod(coordinates.z/dx2-0.5,n);
    int nx_d = std::floor(nx); //d for down value, up is just +1
    int ny_d = std::floor(ny);
    int nz_d = std::floor(nz);
    int nx_u = nx_d+1; //d for down value, up is just +1
    int ny_u = ny_d+1;
    int nz_u = nz_d+1;
    double Dx = nx-nx_d; // decimal part
    double Dy = ny-ny_d;
    double Dz = nz-nz_d;

    return V[nx_d+ny_d*n+nz_d*pow(n,2)]*(1-Dx)*(1-Dy)*(1-Dz)
            +V[nx_d+ny_d*n+nz_u*pow(n,2)]*(1-Dx)*(1-Dy)*(Dz)
                +V[nx_d+ny_u*n+nz_d*pow(n,2)]*(1-Dx)*(Dy)*(1-Dz)
                    +V[nx_d+ny_u*n+nz_u*pow(n,2)]*(1-Dx)*(Dy)*(Dz)
                        +V[nx_u+ny_d*n+nz_d*pow(n,2)]*(Dx)*(1-Dy)*(1-Dz)
                            +V[nx_u+ny_d*n+nz_u*pow(n,2)]*(Dx)*(1-Dy)*(Dz)
                                +V[nx_u+ny_u*n+nz_d*pow(n,2)]*(Dx)*(Dy)*(1-Dz)
                                    +V[nx_u+ny_u*n+nz_u*pow(n,2)]*(Dx)*(Dy)*(Dz);
}

//! Class which sets the initial scalar field matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double background_phi; 
        double background_dphi;
        double hubble;
        double min_chi;
        double min_lapse;
        std::array<double, AMREX_SPACEDIM> center;   //! coordinates center
        double dx_input; // coordinates of the input grid. 
        std::vector<double> delta_phi; // flat array (vector) 
        std::vector<double> delta_Pi; 
        std::vector<double> delta_X; 
        std::vector<double> delta_K; 
        int n; // size of input box // may have to do interpolation later
        
    };

    //! The constructor
    InitialScalarData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t> &cell) const
    {
        MatterCCZ4RHS<ScalarField<>>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.); // Set only the non-zero components below

        
        // where am i?
        amrex::IntVect pos(i, j, k);
        Coordinates<data_t> coords(pos, m_dx, m_params.center);

        // set the field vars //background+ perturbation
        vars.phi = m_params.background_phi+Trilinear_Interpolation(coords, m_dx, m_params.dx_input, m_params.n, m_params.delta_phi);
        vars.Pi = m_params.background_dphi+Trilinear_Interpolation(coords, m_dx, m_params.dx_input, m_params.n, m_params.delta_Pi); 

        // start with unit lapse and flat metric (must be relaxed for chi)
        vars.lapse = 1;
        vars.chi = 1+Trilinear_Interpolation(coords, m_dx, m_params.dx_input, m_params.n, m_params.delta_X);
        vars.Pi = vars.Pi/vars.lapse; // no initial shift

        // conformal metric is flat
        FOR (index)
            vars.h[index][index] = 1.;

        // Inflationnary background
        vars.K = -3*m_params.hubble+Trilinear_Interpolation(coords, m_dx, m_params.dx_input, m_params.n, m_params.delta_K);

        // Store the initial values of the variables
        store_vars(cell.cellData(i, j, k), vars);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* MYINITIALSCALARDATA_HPP_ */
