/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"
// //#include "SixthOrderDerivatives.hpp"

// // For RHS update
#include "MultiSourceCCZ4RHS.hpp"

// // For constraints calculation
#include "Constraints.hpp"
#include "EMTensor.hpp"
#include "NewMatterConstraints.hpp"

// For tagging cells
#include "ChiExtractionTaggingCriterion.hpp"
#include "FixedGridsTaggingCriterion.hpp"

// // Problem specific includes
#include "myInitialScalarData.hpp"
#include "Potential.hpp"
#include "SourceFields.hpp"

#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif


// Power function (not used yet)
void CustomElementWiseOperation(MultiFab& mf, const int& component, double (*func)(double)){
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();  // Get the valid box for the current grid
        auto const& arr = mf.array(mfi); // Access the grid data

        // Apply lambda function element-wise within the box using ParallelFor
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Apply a lambda function to the 6th component of each grid cell
            arr(i, j, k, component) = func(arr(i, j, k, component));  //although not so sure that.
        });
    }
}

double PowNegativeHalf(double x) {
    return std::pow(x, -0.5);
}


// Example of function.
double myFourierAmplitude(double k) {
    double hubble = 1.0e-3;
    if (k<= hubble) {
        return hubble*pow(k,-1.5)/sqrt(2.0);
    } else{
        return 0.0;
        }
}

double Heaviside(double x) {
    return (x >= 0.0) ? 1.0 : 0.0;
}

double sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

// Windowing functions examples:
double dWindow(double dk_t_up, double dk_t_down, std::string type) {
    if (type == "heaviside") {
        return Heaviside(dk_t_up)-Heaviside(dk_t_down);
    } 
    else if (type == "erf") {
        return std::erf(dk_t_up)-std::erf(dk_t_down);
        }
    else if (type == "sigmoid"){
        return sigmoid(dk_t_up)-sigmoid(dk_t_down);
    }
    else if (type == "gaussian"){
        return std::exp(-pow(dk_t_up,2)/2.0)-std::exp(-pow(dk_t_down,2)/2.0);
    }
    }

// Windowing functions examples:
double DWindowDt(double k_ratio, double dlog_ksig_dt, std::string type) {
    if (type == "gaussian"){
        return pow(k_ratio,2.0)*dlog_ksig_dt*std::exp(-pow(k_ratio,2)/2.0);
    }
    }


// for later auto g = [](double x) { return f(x, 3.0); };

void ScalarFieldLevel::variableSetUp()
{
    BL_PROFILE("ScalarFieldLevel::variableSetUp()");

    // Set up the state variables
    stateVariableSetUp();

    const int nghost = simParams().num_ghosts;

    // // Add the constraints to the derive list
    static inline const amrex::Vector<std::string> var_names = {"Ham","Mom", "Ham_abs_terms", "Mom_abs_terms"};
    derive_lst.add(
        "constraints", amrex::IndexType::TheCellType(),
        static_cast<int>(var_names.size()), var_names, //static_cast<int>(Constraints::var_names_norm.size()),Constraints::var_names,
        amrex::DeriveFuncFab(), // null function because we won't use
                                // it.
        [=](const amrex::Box &box) { return amrex::grow(box, nghost); },
        &amrex::cell_quartic_interp);

    // Determine the number of components in State_Type
    const int num_state_components = desc_lst[State_Type].nComp();

    // We only need the non-gauge CCZ4 variables to calculate the constraints
    // derive_lst.addComponent("constraints", desc_lst, State_Type, 0, c_lapse);
    derive_lst.addComponent("constraints", desc_lst, State_Type, 0, num_state_components);
    
    int my_n = static_cast<int>(1);  
    derive_lst.add(
        "density", amrex::IndexType::TheCellType(),
        my_n, {"rho"},
        amrex::DeriveFuncFab(), // null function because we won't use
                                // it.
        [=](const amrex::Box &box) { return amrex::grow(box, nghost); },
        &amrex::cell_quartic_interp);
    
    derive_lst.addComponent("density", desc_lst, State_Type, 0, num_state_components);
    
}

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    BL_PROFILE("ScalarFieldLevel::specificAdvance");
    // Enforce trace free A_ij and positive chi and alpha
    amrex::MultiFab &S_new = get_new_data(State_Type);
    const auto &arrs       = S_new.arrays();

    // Enforce the trace free A_ij condition and positive chi and alpha
    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndAlpha(simParams().initial_params.min_chi, simParams().initial_params.min_lapse)(cell);
                       });

    // Check for nan's
    if (simParams().nan_check)
    {
        if (S_new.contains_nan(0, S_new.nComp(), amrex::IntVect(0), true))
        {
            amrex::Abort("NaN in specificAdvance");
        }
    }
}

// Initial data for field and metric variables
void ScalarFieldLevel::initData()
{
    BL_PROFILE("ScalarFieldLevel::initData");
    if (m_verbosity)
        amrex::Print() << "ScalarFieldLevel::initialData " << Level()
                       << std::endl;

    const auto dx = geom.CellSizeArray();
    // InitialScalarData gaussian_pulse(simParams().initial_params, dx[0]);
    InitialScalarData my_initial_data(simParams().initial_params, dx[0]);

    amrex::MultiFab &state  = get_new_data(State_Type);
    auto const &state_array = state.arrays();

    amrex::ParallelFor(
        state, state.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_ind, int i, int j, int k) noexcept
        {
            amrex::CellData<amrex::Real> cell =
                state_array[box_ind].cellData(i, j, k);
            for (int n = 0; n < cell.nComp(); ++n)
            {
                cell[n] = 0.;
            }

            my_initial_data.compute(i, j, k, state_array[box_ind]);
        });

    if (simParams().nan_check)
    {
        if (state.contains_nan(0, state.nComp(), amrex::IntVect(0), true))
        {
            amrex::Abort("NaN in initData");
        }
    }
    //Initialize stochastic engine/seed
    initializeRandomEngine();
    // Setting up the spectral modifier and its random grid.
    initializeStochasticMultiFabs(state);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                       amrex::MultiFab &a_rhs,
                                       const double a_time)
{
    BL_PROFILE("ScalarFieldLevel::specificEvalRHS()");

    // std::string plot_file_root =
    //     "/home/dc-kwan1/rds/rds-dirac-dp002/dc-kwan1/GRTeclyn/"
    //     "ScalarField/test_file";

    amrex::Vector<std::string> var_names;
    for (int i = 0; i < NUM_VARS; i++)
    {
        if (a_soln.contains_nan(i, 1, amrex::IntVect(0), true))
            amrex::Print() << "Nan found in component " << i << std::endl;
        var_names.push_back(StateVariables::names[i]);
    }

    // amrex::WriteSingleLevelPlotfileHDF5(plot_file_root, a_rhs,
    // 				      var_names, Geom(), 0.0, 0);

    const auto &soln_arrs   = a_soln.arrays();
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();

    // Enforce positive chi and alpha and trace free A
    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndAlpha(simParams().initial_params.min_chi, simParams().initial_params.min_lapse)(cell);
                       });

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(simParams().potential_params);
    SourceFieldsWithPotential source_fields(potential);

    // Calculate CCZ4 right hand side
    if (simParams().max_spatial_derivative_order == 4)
    {
        MultiSourceCCZ4RHS<SourceFieldsWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            matter_ccz4_rhs(source_fields, simParams().ccz4_params,
                            Geom().CellSize(0), simParams().sigma,
                            simParams().formulation, simParams().G_Newton);
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) {
                matter_ccz4_rhs.compute(i, j, k, rhs_arrs[box_no],
                                        soln_c_arrs[box_no]);
            });
    }
    else if (simParams().max_spatial_derivative_order == 6)
    {
        amrex::Abort("xxxxx max_spatial_derivative_order == 6 todo");
#if 0
        MultiSourceCCZ4RHS<SourceFieldsWithPotential, MovingPunctureGauge, SixthOrderDerivatives>
	  matter_ccz4_rhs(source_fields, simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
			  simParams().formulation, simParams().G_Newton);
        amrex::ParallelFor(a_rhs,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
        {
            amrex::CellData<amrex::Real const> state = soln_c_arrs[box_no].cellData(i,j,k);
            amrex::CellData<amrex::Real> rhs = rhs_arrs[box_no].cellData(i,j,k);
            matter_ccz4_rhs.compute(i,j,k,rhs_arrs[box_no], soln_c_arrs[box_no]);
        });
#endif
    }

    /////////////////////// STOCHASTIC SOURCES ////////////////////////////////////
    specificEvalStochasticRHS(a_soln,a_rhs);


    if (simParams().nan_check)
    {   
        if (a_soln.contains_nan(0, a_soln.nComp(), amrex::IntVect(0), true))
        {

            // const std::string& pltfile = amrex::Concatenate(plot_file_root,
            //                                           file_name_digits);

            amrex::Vector<std::string> var_names;
            for (int i = 0; i < NUM_VARS; i++)
            {
                if (a_soln.contains_nan(i, 1, amrex::IntVect(0), true))
                    amrex::Print()
                        << "Nan found in component " << i << std::endl;
                //	      var_names.push_back(StateVariables::names[i]);
            }

            // amrex::WriteSingleLevelPlotfileHDF5(plot_file_root, a_rhs,
            // 				      var_names, Geom(), 0.0, 0);

            //	  amrex::WriteSingleLevelPlotfile(pltfile);
            amrex::Abort("NaN in specificUpdateRHS");
        }
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    BL_PROFILE("ScalarFieldLevel::specificUpdateODE()");
    // Enforce the trace free A_ij condition
    const auto &soln_arrs = a_soln.arrays();
    amrex::ParallelFor(a_soln, amrex::IntVect(0), // zero ghost cells
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                       });
}

void ScalarFieldLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void ScalarFieldLevel::errorEst(amrex::TagBoxArray &tagging_criterion,
                                int /*clearval*/, int /*tagval*/,
                                amrex::Real /*time*/, int /*n_error_buf*/,
                                int /*ngrow*/)

{
    BL_PROFILE("ScalarFieldLevel::errorEst()");

    amrex::MultiFab &state_new = get_new_data(State_Type);
    const auto curr_time       = get_state_data(State_Type).curTime();

    const int nghost =
        state_new.nGrow(); // Need ghost cells to compute gradient
    const int ncomp = state_new.nComp();

    // I filled all the ghost cells in case but could also just fill the ones
    // used for tagging We only use chi in the tagging criterion so only fill
    // the ghosts for chi
    FillPatch(*this, state_new, nghost, curr_time, State_Type, 0, ncomp);

    const auto &simpar = simParams();

    const auto &tags           = tagging_criterion.arrays();
    const auto &state_new_arrs = state_new.const_arrays();
    const auto tagval          = amrex::TagBox::SET;

    amrex::Real dx     = Geom().CellSize(0);
    int curr_level     = Level();
    const auto probhi  = Geom().ProbHiArray();
    const auto problo  = Geom().ProbLoArray();
    amrex::Real length = probhi[0] - problo[0];

    const auto test_dx = Geom().CellSizeArray();

    FixedGridsTaggingCriterion tagger(test_dx[0], curr_level, length,
                                      simParams().initial_params.center);
    // ChiExtractionTaggingCriterion tagger(Geom().CellSize(0), Level(),
    //                                      simpar.extraction_params,
    //                                      simpar.activate_extraction);

    amrex::Real threshold = simpar.regrid_thresholds[Level()];
    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::Real criterion =
                               tagger.compute(i, j, k, state_new_arrs[box_no]);

                           // amrex::Real criterion =
                           //     tagger(i, j, k, state_new_arrs[box_no]);

                           if (criterion >= threshold)
                           {
                               tags[box_no](i, j, k) = tagval;
                           }
                       });

    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::derive(const std::string &name, amrex::Real time,
                              amrex::MultiFab &multifab, int dcomp)
{
    BL_PROFILE("ScalarFieldLevel::derive()");

    BL_ASSERT(dcomp < multifab.nComp());

    const int num_ghosts = multifab.nGrow();

    const amrex::DeriveRec *rec = derive_lst.get(name);
    if (rec != nullptr)
    {
        int state_idx, derive_scomp, derive_ncomp;

        // we only have one state so state_idx will be State_Type = 0
        rec->getRange(0, state_idx, derive_scomp, derive_ncomp);

        // work out how many extra ghost cells we need
        const amrex::BoxArray &src_ba = state[state_idx].boxArray();

        int num_extra_ghosts = num_ghosts;
        {
            amrex::Box box0   = src_ba[0];
            amrex::Box box1   = rec->boxMap()(box0);
            num_extra_ghosts += box0.smallEnd(0) - box1.smallEnd(0);
        }

        // Make a Multifab with enough extra ghosts to calculated derived
        // quantity. For now use NUM_VARS in case the enum mapping loads more
        // vars than is actually needed
        amrex::MultiFab src_mf(src_ba, dmap, NUM_VARS, num_extra_ghosts,
                               amrex::MFInfo(), *m_factory);

        // Fill the multifab with the needed state data including the ghost
        // cells
        FillPatch(*this, src_mf, num_extra_ghosts, time, state_idx,
                  derive_scomp, derive_ncomp);

        const auto &src_arrays = src_mf.const_arrays();

        Potential potential(simParams().potential_params);
        SourceFieldsWithPotential source_fields(potential);

        if (name == "constraints")
        {
            const auto &out_arrays = multifab.arrays();
            // Interval imom = Interval(dcomp + 1, dcomp + AMREX_SPACEDIM); // This neds to be modified to same begin and end to get sqrt Mom^2
            // MatterConstraints<SourceFieldsWithPotential> constraints(
                // source_fields, Geom().CellSize(0), simParams().G_Newton, iham,
                // imom); // This needs to also have an input of iham_abs and imom_abs
            int iham = dcomp;
            Interval imom = Interval(dcomp + 1, dcomp + 1);
            int iham_abs = dcomp+2;
            Interval imom_abs = Interval(dcomp + 3, dcomp + 3);

            MatterConstraints<SourceFieldsWithPotential> constraints(
                source_fields, Geom().CellSize(0), simParams().G_Newton, iham,
                imom, iham_abs, imom_abs);
            amrex::ParallelFor(
                multifab, multifab.nGrowVect(),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
                    constraints.compute(i, j, k, out_arrays[box_no],
                                        src_arrays[box_no]);
                });
        }
        else if (name == "density")
        {
            const auto &out_arrays = multifab.arrays();

            EMTensor<SourceFieldsWithPotential> emtensor(
                source_fields, Geom().CellSize(0), dcomp);

            amrex::ParallelFor(
                multifab, multifab.nGrowVect(),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
                    emtensor.compute(i, j, k, out_arrays[box_no],
                                  src_arrays[box_no]);
                });
        }
        else
        {
            amrex::Abort("Unknown derived variable");
        }
    }
    else
    {
        amrex::Abort("Unknown derived variable");
    }
    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::initializeRandomEngine()
{
    std::random_device rd;
    seed = rd(); // Generate a random seed
    random_engine.seed(seed); // Initialize the random engine with the seed

    amrex::Print() << "Initialized random engine with seed: " << seed << std::endl;
}

void ScalarFieldLevel::initializeStochasticMultiFabs(const MultiFab& existing_mf) {

    int ncomp = 1;               // Only one component
    int ngrow = 0; //existing_mf.nGrow(); // Number of ghost cells same as existing_mf

    gaussian_grid.define(existing_mf.boxArray(), existing_mf.DistributionMap(), ncomp, ngrow);
    stochastic_rhs_R.define(existing_mf.boxArray(), existing_mf.DistributionMap(), ncomp, ngrow);
    // stochastic_rhs_Pi.define(existing_mf.boxArray(), existing_mf.DistributionMap(), ncomp, ngrow);
    // stochastic_rhs_K.define(existing_mf.boxArray(), existing_mf.DistributionMap(), ncomp, ngrow);
    buffer.define(existing_mf.boxArray(), existing_mf.DistributionMap(), ncomp, existing_mf.nGrow());

    gaussian_grid.setVal(0.0);
    stochastic_rhs_R.setVal(0.0);
    // stochastic_rhs_Pi.setVal(0.0);
    // stochastic_rhs_K.setVal(0.0);
    buffer.setVal(0.0);

    spectral_modifier = SpectralModifier(gaussian_grid, geom, 0);

    if (StateVariables::names[comp_Pi] != "Pi") 
        {// {amrex::Print() << "Current comp  of Pi is " << StateVariables::names[comp_Pi] << std::endl;
        amrex::Abort("Wrong component for Pi.");
        }
    if (StateVariables::names[comp_K] != "K" )
        {// {amrex::Print() << "Current comp of K is  " << StateVariables::names[comp_K] << std::endl;

        amrex::Abort("Wrong component for K.");
        }
    if (StateVariables::names[comp_chi] != "chi" )
        {// {amrex::Print() << "Current comp of K is  " << StateVariables::names[comp_K] << std::endl;

        amrex::Abort("Wrong component for chi.");
        }
    
    for (MFIter mfi(existing_mf); mfi.isValid(); ++mfi) {
        num_cells += double(mfi.validbox().numPts());  // Counting the valid cells in the MultiFab, used later for averaging
    }

    Mpl = pow(8*PI*simParams().G_Newton,-0.5);
    amrex::Print() << "Mpl = " << Mpl << std::endl;
    sigma = simParams().initial_params.sigma;
    amrex::Print() << "sigma = " << sigma << std::endl;
    spectral_modifier.FillInputWithRandomNoise(random_engine); // NEW /!!

}

void ScalarFieldLevel::specificEvalStochasticRHS(amrex::MultiFab &a_soln,amrex::MultiFab &a_rhs){

    // Background quantities for factors //ignore non-commutativity of <> and operations for now.
    amrex::Real current_dt = parent->dtLevel(level);

    double av_a = pow(a_soln.sum(comp_chi,true)/num_cells,-0.5);
    double av_H = a_soln.sum(comp_K, true)/num_cells/-3.0; // average H without ghost cells
    double av_inv_rad = av_a*av_H; //inverse Hubble radius
    double av_eps1 = a_rhs.sum(comp_K, true)/num_cells/pow(av_H,2.0)/3.0;  //see 14/08/24 draft 6, lapse = 1 for now.
    // amrex::Print() << "Average computations: " << std::endl;
    // spectral_modifier.FillInputWithRandomNoise(random_engine);

    // Build the window function's boundaries.
    double k_cross = sigma*av_inv_rad;
    double delta_k = k_cross*av_H*current_dt*(1.0-av_eps1)/2.0; // see draft 6, 9/9/24
    // double ksup = k_cross+delta_k;
    // double kinf = k_cross-delta_k;
    double ksup = k_cross+2.0*delta_k;
    double kinf = k_cross;
    amrex::Print() << "kinf = " << kinf << std::endl;
    amrex::Print() << "ksup = " << ksup << std::endl;

    auto CurrentFourierAmplitude = [=](double k) -> double {
        // if (k < ksup && k > kinf) {
        //     // double window_factor = std::sqrt(av_H*k_cross)/current_dt;
        //     double window_factor = 1/current_dt;
        //     // return av_H/Mpl * std::pow(k, -1.5) / std::sqrt(4.0*av_eps1)*window_factor;
        //     return 3.0*std::pow(av_H, 2.0)/Mpl * std::pow(k, -1.5) / std::sqrt(4.0*av_eps1)*window_factor;
        // } else {
        //     return 0.0;
        // }

        // return 3.0*std::pow(av_H, 2.0)/Mpl * std::pow(k, -1.5) / std::sqrt(4.0*av_eps1)*dWindow(ksup-k,kinf-k,"heaviside")/current_dt;
        // return 3.0*std::pow(av_H, 2.0)/Mpl * std::pow(k, -1.5) / std::sqrt(4.0*av_eps1)*dWindow(k/ksup,k/kinf,"gaussian")/sqrt(current_dt);

        // return 3.0*std::pow(av_H, 2.0)/Mpl * std::pow(k, -1.5) / std::sqrt(4.0*av_eps1)*dWindow(k/ksup,k/kinf,"gaussian")* std::pow(current_dt, -1.5);
        // return 3.0*std::pow(av_H, 2.0)/Mpl * std::pow(k, -1.5) / std::sqrt(4.0*av_eps1)*DWindowDt(k/kinf,av_H*(1.0-av_eps1),"gaussian")* std::pow(current_dt, -0.5);
        // return 3.0*std::pow(av_H, 2.0)/Mpl * std::pow(k, -1.5) / std::sqrt(4.0*av_eps1)*dWindow(ksup-k,kinf-k,"heaviside")* std::pow(current_dt, -1.5);

        return 3.0*std::pow(av_H, 2.0)/Mpl * std::pow(k, -1.5) / std::sqrt(4.0*av_eps1)*DWindowDt(k/kinf,av_H*(1.0-av_eps1),"gaussian");
    };

    spectral_modifier.apply_func(CurrentFourierAmplitude, stochastic_rhs_R);

    // Wiener process scaling of the noise, before discrete *dt. See Specs.
    // stochastic_rhs_R.mult(pow(current_dt,-0.5), 0, 1, 0);
    // stochastic_rhs_R.mult(pow(current_dt,-1), 0, 1, 0);

    // Make stoch. rhs of Pi
    stochastic_rhs_R.mult(-1.0*pow(2.0*av_eps1,0.5)*Mpl, 0, 1, 0);
    amrex::MultiFab::Add(a_rhs, stochastic_rhs_R, 0, comp_Pi, 1, 0);

    // Make stoch. rhs of K from the previous one
    stochastic_rhs_R.mult(-1.0*pow(av_eps1/2.0,0.5)/Mpl, 0, 1, 0);

    amrex::MultiFab::Add(a_rhs, stochastic_rhs_R, 0, comp_K, 1, 0);
    // amrex::Print() << "Stochastic source  added" << std::endl;

    // const int nghost = a_rhs.nGrow(); // Need ghost cells to compute gradient
    // const auto curr_time = get_state_data(State_Type).curTime();

    // FillPatch(*this, a_rhs, nghost, curr_time, State_Type, comp_Pi, 1);
    // FillPatch(*this, a_rhs, nghost, curr_time, State_Type, comp_K, 1);
    // a_rhs.FillBoundary(geom.periodicity());

}