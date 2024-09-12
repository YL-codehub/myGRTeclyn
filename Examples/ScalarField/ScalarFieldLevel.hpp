/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDLEVEL_HPP_
#define SCALARFIELDLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "DefaultPotential.hpp"
#include "Potential.hpp"
#include "SourceFields.hpp"
#include <random> // random numbers generation.
#include "SpectralModifier.hpp"  // applying a spectrum to a random array
#define PI (3.14159265358979323846264338327950288E0)

//!  A class for the evolution of a scalar field, minimally coupled to gravity
/*!
     The class takes some initial data for a scalar field (variables phi and Pi)
     and evolves it using the CCZ4 equations. It is possible to specify an
   initial period of relaxation for the conformal factor chi, for non analytic
   initial conditions (for example, a general field configuration at a moment of
   time symmetry assuming conformal flatness). \sa MatterCCZ4(),
   ConstraintsMatter(), ScalarField(), RelaxationChi()
*/
class ScalarFieldLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ScalarFieldLevel>;

  public:

    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    static void variableSetUp();

    // Typedef for scalar field
    typedef SourceFields<Potential> SourceFieldsWithPotential;

    using DefaultScalarField = SourceFields<DefaultPotential>;

    //! Things to do at the end of the advance step, after RK4 calculation
    void specificAdvance() override;

    //! Initialize data for the field and metric variables

    void initData() override;

    //! RHS routines used at each RK4 step
    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override;
                         
    void specificEvalStochasticRHS(amrex::MultiFab &a_soln,amrex::MultiFab &a_rhs);

    //! Things to do in UpdateODE step, after soln + rhs update
    void specificUpdateODE(amrex::MultiFab &a_soln) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    void preTagCells();

    //! Tell GRTeclyn how to tag cells for regridding
    void errorEst(amrex::TagBoxArray &tag_box_array, int clearval, int tagval,
                  amrex::Real time, int n_error_buf = 0, int ngrow = 0) final;

    void derive(const std::string &name, amrex::Real time,
                amrex::MultiFab &multifab, int dcomp) override;
  
  private:
    // Random draw tools
    std::mt19937 random_engine;
    unsigned int seed;
    void initializeRandomEngine(); // Method to initialize the random seed
    // Stochastic grids (Stored on main node)
    amrex::MultiFab gaussian_grid;
    amrex::MultiFab stochastic_rhs_R; // see eq. (47-48) of 10.1103/PhysRevD.109.123523 
    amrex::MultiFab stochastic_rhs_Pi;
    amrex::MultiFab stochastic_rhs_K;
    void initializeStochasticMultiFabs(const MultiFab& existing_mf);
    // spectral method class
    SpectralModifier spectral_modifier;
    // utils
    int comp_Pi = 26; //This should be the momentum's var num
    int comp_K = 7; // same for extrinsic curvature trace
    int comp_chi = 0; //This should be the conformal factor's var num
    double num_cells = 0.0;
    double Mpl = 1.0; //Reduced Planck mass
    double sigma = 1.0;

};

#endif /* SCALARFIELDLEVEL_HPP_ */
