/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "myInitialScalarData.hpp"
#include "Potential.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // // Initial scalar field data
        // initial_params.center =
        //     center; // already read in SimulationParametersBase
        // pp.load("G_Newton", G_Newton,
        //         0.0); // for now the example neglects backreaction
        // pp.load("scalar_amplitude", initial_params.amplitude, 0.1);
        // pp.load("scalar_width", initial_params.width, 1.0);
        // pp.load("scalar_mass", potential_params.scalar_mass, 0.1);
        // Initial scalar field data
        initial_params.center = center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton,1.0);
        pp.load("background_phi", initial_params.background_phi, 0.0);
        pp.load("background_dphi", initial_params.background_dphi, 0.0);
        pp.load("hubble", initial_params.hubble, 0.0);
        pp.load("min_chi", initial_params.min_chi, 1e-4);
        pp.load("min_lapse", initial_params.min_lapse, 1e-4);
        pp.load("potential_param_1", potential_params.potential_param_1, 0.0);
        pp.load("potential_param_2", potential_params.potential_param_2, 1.0);
        pp.load("potential_param_3", potential_params.potential_param_3, 1.0);
        pp.load("potential_param_4", potential_params.potential_param_4, 1.0);
        pp.load("potential_param_5", potential_params.potential_param_5, 1.0);
        pp.load("potential_type", potential_params.potential_type, 1);

        // Pertuabtion initial data
        // pp.load("dx_input",initial_params.dx_input,1.); //spacing input grid
        // pp.load("n",initial_params.n,1); //size per dimension
        // int ntot=pow(initial_params.n,3);
        // pp.load("delta_phi",initial_params.delta_phi, ntot); 
        // pp.load("delta_Pi",initial_params.delta_Pi,  ntot); 
        // pp.load("delta_X",initial_params.delta_X,  ntot);
        // pp.load("delta_K",initial_params.delta_K,  ntot); 
    }

    void check_params()
    {
        // warn_parameter("scalar_mass", potential_params.scalar_mass,
        //                potential_params.scalar_mass <
        //                    0.2 / coarsest_dx / dt_multiplier,
        //                "oscillations of scalar field do not appear to be "
        //                "resolved on coarsest level");
        // warn_parameter("scalar_width", initial_params.width,
        //                initial_params.width < 0.5 * L,
        //                "is greater than half the domain size");
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    InitialScalarData::params_t initial_params;
    Potential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
