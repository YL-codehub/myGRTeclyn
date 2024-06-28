/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

// auxiliary functions
// // resonance
template <class data_t> data_t bump(data_t x, const double &phi_f,const double &d) { return 1-pow((tanh((phi_f-x)/d)),2); };

template <class data_t> data_t dbump(data_t x, const double &phi_f,const double &d) { return -2*tanh((phi_f-x)/d)*bump(x,phi_f,d)*(-1./d);};

template <class data_t> data_t osc(data_t x, const double &phi_f, const double &c, const double &f) {return  c*sin((x-phi_f)/f);};

template <class data_t> data_t dosc(data_t x, const double &phi_f, const double &c, const double &f) {return c*cos((x-phi_f)/f)/f;};

// // axion monodromy

template <class data_t> data_t wind(data_t x, const double &phi_f,const double &d, const double &f) { return  0.25 * (1 + tanh((x - phi_f) / f)) * (1 + tanh((phi_f - x + d) / f)); };

template <class data_t> data_t dwind(data_t x, const double &phi_f,const double &d, const double &f) { return  0.25 / f * (1 - pow(tanh((x - phi_f) / f), 2)) * (1 + tanh((phi_f - x + d) / f)) - 0.25 / f * (1 + tanh((x - phi_f) / f)) * (1 - pow(tanh((phi_f - x + d) / f),2));};

template <class data_t> data_t osc2(data_t x, const double &phi_f, const double &f) {return  (cos((x - phi_f) / f) - 1);};

template <class data_t> data_t dosc2(data_t x, const double &phi_f, const double &f) {return -sin((x - phi_f) / f) / f;};

class Potential
{
  public:
    struct params_t
    {
        // double scalar_mass;
        double potential_param_1;
        double potential_param_2;
        double potential_param_3;
        double potential_param_4;
        double potential_param_5;
        int potential_type;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const vars_t<data_t> &vars) const
    {
        // // The potential value at phi
        // // 1/2 m^2 phi^2
        // V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0);

        // // The potential gradient at phi
        // // m^2 phi
        // dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;
        if (m_params.potential_type == 1){
        // Quadratic type 
        V_of_phi = 0.5 * pow(m_params.potential_param_1 * vars.phi, 2.0);
        dVdphi = pow(m_params.potential_param_1, 2.0) * vars.phi;
        }
        else if (m_params.potential_type == 2){
          // cubic type
        V_of_phi = m_params.potential_param_1 *(1-pow(vars.phi/m_params.potential_param_2,3));
        dVdphi = -3*m_params.potential_param_1*pow(vars.phi,2)/pow(m_params.potential_param_2,3);
        }
        else if (m_params.potential_type == 3){
        // Quadratic type + bumps
        // 3.141592653589793238462643383279502884197169399375105820974944
        V_of_phi = 0.5 * pow(m_params.potential_param_1 * vars.phi, 2.0)+m_params.potential_param_3*(1-cos(2*3.1415926535*(vars.phi-m_params.potential_param_2)/m_params.potential_param_4));
        dVdphi = pow(m_params.potential_param_1, 2.0) * vars.phi+2*3.1415926535*m_params.potential_param_3/m_params.potential_param_4*sin(2*3.1415926535*(vars.phi-m_params.potential_param_2)/m_params.potential_param_4);
        }
        // else if (m_params.potential_type == 4){
        // // Quadratic type + 1 bump
        // V_of_phi = 0.5 * pow(m_params.potential_param_1 * vars.phi, 2.0)+m_params.potential_param_3*pow(vars.phi-m_params.potential_param_2,2)*exp(-pow(vars.phi-m_params.potential_param_2,2)/m_params.potential_param_4);
        // dVdphi = pow(m_params.potential_param_1, 2.0) * vars.phi+2*m_params.potential_param_3*(vars.phi-m_params.potential_param_2)*(1-pow(vars.phi-m_params.potential_param_2,2)/m_params.potential_param_4)*exp(-pow(vars.phi-m_params.potential_param_2,2)/m_params.potential_param_4);
        // }
        else if (m_params.potential_type == 4){
        // Quadratic type + 1 bump
        V_of_phi = 0.5 * pow(m_params.potential_param_1 * vars.phi, 2.0)*(1+m_params.potential_param_3*exp(-0.5*pow(vars.phi-m_params.potential_param_2,2)/pow(m_params.potential_param_4,2)));
        dVdphi = pow(m_params.potential_param_1, 2.0) * vars.phi*(1+m_params.potential_param_3*exp(-0.5*pow(vars.phi-m_params.potential_param_2,2)/pow(m_params.potential_param_4,2)))
                        + 0.5 * pow(m_params.potential_param_1 * vars.phi, 2.0)*(-0.5*(vars.phi-m_params.potential_param_2)/pow(m_params.potential_param_4,2))*m_params.potential_param_3*exp(-0.5*pow(vars.phi-m_params.potential_param_2,2)/pow(m_params.potential_param_4,2));
        }
        else if (m_params.potential_type == 5){
        // Quadratic type 
        V_of_phi = m_params.potential_param_1 * pow(vars.phi, m_params.potential_param_2);
        dVdphi =  m_params.potential_param_2*m_params.potential_param_1 * pow(vars.phi, m_params.potential_param_2-1);
        }
        else if (m_params.potential_type == 6){
        // USR+slope type.
        // bool where = (vars.phi <= m_params.potential_param_2);
        auto condition = simd_compare_lt(vars.phi, m_params.potential_param_2);
        V_of_phi = simd_conditional(condition,m_params.potential_param_1, m_params.potential_param_1* (3*pow(vars.phi/m_params.potential_param_2,2)-2*pow(vars.phi/m_params.potential_param_2,3)));
        dVdphi =  simd_conditional(condition,0.0, m_params.potential_param_1* 6*(vars.phi/pow(m_params.potential_param_2,2)-pow(vars.phi/m_params.potential_param_2,2)/m_params.potential_param_2));
        }
         else if (m_params.potential_type == 7){
        // oscillatory quadratic potential
        data_t B = bump(vars.phi,m_params.potential_param_2,m_params.potential_param_3);
        data_t dB = dbump(vars.phi,m_params.potential_param_2,m_params.potential_param_3);
        data_t O = osc(vars.phi,m_params.potential_param_2,m_params.potential_param_4,m_params.potential_param_5);
        data_t dO = dosc(vars.phi,m_params.potential_param_2,m_params.potential_param_4,m_params.potential_param_5);
        V_of_phi = 0.5 * pow(m_params.potential_param_1 * vars.phi, 2.0)*(1.+B*O);
        dVdphi = pow(m_params.potential_param_1, 2.0) * vars.phi*(1.+B*O)+V_of_phi*(dB*O+B*dO);
        }
        else if (m_params.potential_type == 8){
        // Prokopec & Germani 2017 SR + USR
        data_t phi2 = pow(vars.phi,2);
        V_of_phi = m_params.potential_param_1*pow(m_params.potential_param_2,4)*phi2/3*(3*phi2+2*sqrt(2)*vars.phi*m_params.potential_param_2+6*pow(m_params.potential_param_2,2))/pow(3*phi2+2*pow(m_params.potential_param_2,2),2);
        dVdphi = 2*m_params.potential_param_1*pow(m_params.potential_param_2,5)*vars.phi*(-sqrt(2)*vars.phi-2*m_params.potential_param_2)*(phi2-2*pow(m_params.potential_param_2,2))/pow(3*phi2+2*pow(m_params.potential_param_2,2),3);
        }
        else if (m_params.potential_type == 9){
        // axion monodromy
        data_t W = wind(vars.phi,m_params.potential_param_2,m_params.potential_param_3,m_params.potential_param_5);
        data_t dW = dwind(vars.phi,m_params.potential_param_2,m_params.potential_param_3,m_params.potential_param_5);
        data_t O = osc2(vars.phi,m_params.potential_param_2,m_params.potential_param_5);
        data_t dO = dosc2(vars.phi,m_params.potential_param_2,m_params.potential_param_5);
        V_of_phi = 0.5 * pow(m_params.potential_param_1 * vars.phi, 2.0)+m_params.potential_param_4*W*O;
        dVdphi = pow(m_params.potential_param_1, 2.0) * vars.phi+m_params.potential_param_4*(dW*O+W*dO);
        }
        else{
        V_of_phi = 0.0;
        dVdphi = 0.0;
        }
    }
};

#endif /* POTENTIAL_HPP_ */
