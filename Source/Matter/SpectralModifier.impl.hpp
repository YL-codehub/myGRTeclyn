#if !defined(SPECTRALMODIFIER_HPP)
#error "This file should only be included through SpectralModifier.hpp"
#endif

#ifndef SPECTRALMODIFIER_IMPL_HPP
#define SPECTRALMODIFIER_IMPL_HPP
#define PI (3.14159265358979323846264338327950288E0)


#include "SpectralModifier.hpp"


// Default Constructor
inline SpectralModifier::SpectralModifier (): 
                  input(nullptr), geom(), verbose(), ba(nullptr), dmap(),
                  n_box_dim(), domain(), numb_boxes_dim(),
                  numb_boxes(), rank_mapping(), h(), a(), b()
    {
    }

// Constructor
inline SpectralModifier::SpectralModifier (MultiFab& input, Geometry& geom, int verbose): 
                  input(&input), geom(geom), verbose(verbose), ba(&(input.boxArray())), dmap(input.DistributionMap()),
                  n_box_dim((*ba)[0].size()), domain(geom.Domain()), numb_boxes_dim(domain.length() / n_box_dim),
                  numb_boxes(numb_boxes_dim[0]*numb_boxes_dim[1]*numb_boxes_dim[2]), rank_mapping(), h(geom.CellSize(0)), a(), b()
    {

    if (input.nGrow() != 0)
       amrex::Error("Current implementation requires that both input and output have no ghost cells");
      // Ericka said fillGhosts could be useful to use to fill periodic Ghost cells once I have filled in with my noise.

    if (numb_boxes != (*ba).size()) { // size of Box object = nx*ny*nz
       amrex::Print() <<  ba <<std::endl;
       amrex::Print() << numb_boxes << std::endl;
       amrex::Error("numb_boxes not right");}

    rank_mapping.resize(numb_boxes); 
    remap();
  
    a.resize(n_box_dim[0]*n_box_dim[1]*n_box_dim[2]);
    b.resize(n_box_dim[0]*n_box_dim[1]*n_box_dim[2]); 
    // amrex::Print() << "dx?"<< h << std::endl;

    }

inline void SpectralModifier::remap()
{
    //Building the FFT's new C-ordered rank mapping.
    // AMREX has Fortran/ column-major ordering (first index, ie x, changing fastest in memory: arr[i_x][i_y][i_z] --> arr[i_x+ n_x*i_y+n_x*n_y*i_z])
    // while dfft has C/row-major ordering (last index, ie x, changing fastest in memory: arr[i_x][i_y][i_z] --> arr[i_z+ n_z*i_y+n_z*n_y*i_x])   

    for (int ib = 0; ib < numb_boxes; ++ib) // For each subbox
    {
        // The smallEnd method provides the indices of the lower corner of the Box in each dimension.
        int i = (*ba)[ib].smallEnd(0) / n_box_dim[0]; 
        int j = (*ba)[ib].smallEnd(1) / n_box_dim[1];
        int k = (*ba)[ib].smallEnd(2) / n_box_dim[2];

        // This would be the "correct" local index if the data wasn't being transformed
        // int local_index = k*nbx*nby + j*nbx + i;
        int local_index = k*numb_boxes_dim[0]*numb_boxes_dim[1] + j*numb_boxes_dim[1] + i;
        // This is what we pass to dfft to compensate for the Fortran ordering
        //      of amrex data in MultiFabs.
        // int local_index = i*numb_boxes_dim[1]*numb_boxes_dim[2] + j*numb_boxes_dim[2] + k;

        rank_mapping[local_index] = dmap[ib]; // copy each rank in the F-ordered Distribution Mapping to its C-order version.
        if (verbose)
          amrex::Print() << "LOADING RANK NUMBER " << dmap[ib] << " FOR GRID NUMBER " << ib
                         << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
    }
}

inline MultiFab SpectralModifier::apply_func(double (*amp_func)(double))
    {
    MultiFab output((*input).boxArray(), (*input).DistributionMap(), (*input).nComp(), (*input).nGrow());
    output.ParallelCopy((*input));

    // Assume for now that nx = ny = nz
    int Ndims[3] = {numb_boxes_dim[2], numb_boxes_dim[1], numb_boxes_dim[0]};
    int     n[3] = {domain.length(2), domain.length(1), domain.length(0)}; // C-order

   // Hardware/Hybrid Cosmology Code (HACC) provides a functionality for performing Discrete Fast Fourier Transforms in AMReX (tailored for cosmological simulations.)
    hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
    hacc::Dfft dfft(d); // creates the dfftiser


    for (MFIter mfi((*input),false); mfi.isValid(); ++mfi) // iterating over the FArrayBox objects (subboxes = Box) that are stored locally on the current MPI rank.
    {
       int gid = mfi.index();

       size_t local_size  = dfft.local_size();


       dfft.makePlans(&a[0],&b[0],&a[0],&b[0]); // designs the fft to be optimal given inputs and outputs dims

       // *******************************************
       // Copy real data from input into real part of a -- no ghost cells and
       // put into C++ ordering (not Fortran)
       // *******************************************
       complex_t zero(0.0, 0.0);
       size_t local_indx = 0;
       for(size_t k=0; k<(size_t)n_box_dim[2]; k++) {
        for(size_t j=0; j<(size_t)n_box_dim[1]; j++) {
         for(size_t i=0; i<(size_t)n_box_dim[0]; i++) {

           complex_t temp((*input)[mfi].dataPtr()[local_indx],0.);
           a[local_indx] = temp;
           local_indx++;

         }
       }
      }

       dfft.forward(&a[0]); // surely this is done in parallel
    //  *******************************************
//  Now apply amp_func to the coefficients of the transform
//  *******************************************
// Warning: DFT such that kx,ky,kz range from 0 to N-1
    local_indx = 0;
    const int *self = dfft.self_kspace(); // rank location in k-space?
    const int *local_ng = dfft.local_ng_kspace();  // local grid dimensions in k-space?
    const int *global_ng = dfft.global_ng();

   // global_i/j/k are the global indices for the k modes.
    for(size_t i=0; i<(size_t)local_ng[0]; i++) {
     size_t global_i = local_ng[0]*self[0] + i;

     for(size_t j=0; j<(size_t)local_ng[1]; j++) {
      size_t global_j = local_ng[1]*self[1] + j;

      for(size_t k=0; k<(size_t)local_ng[2]; k++) {
        size_t global_k = local_ng[2]*self[2] + k;

        if (global_i == 0 && global_j == 0 & global_k == 0) {
           a[local_indx] = 0;
        } else {
         // This is k but remember that gradients in discrete space are not that easy. It's the effective k (see Caravano21), not a problem if there is a cutoff to keep global_i,j,k<<global_ng
           double fac = sqrt(-2. * (
                        (cos(2*PI*double(global_i)/double(global_ng[0])) - 1.) +
                        (cos(2*PI*double(global_j)/double(global_ng[1])) - 1.) +
                        (cos(2*PI*double(global_k)/double(global_ng[2])) - 1.) ));
         // dimensionful effective k is fac/dx
           a[local_indx] = a[local_indx] * std::abs(amp_func(fac/h)); // only keep the square root of the power spectrum.
        }
        local_indx++;

      }
     }
    }

    dfft.backward(&a[0]);

    size_t global_size  = dfft.global_size(); // is this Nx*Ny*Nz? yes
  //  double fac = pow(h,2) / global_size; // why is that how it was orignally normalised? because laplacian?
    double fac = pow(sqrt(2*PI)/h,1.5)/sqrt(global_size);

    local_indx = 0;
    for(size_t k=0; k<(size_t)n_box_dim[2]; k++) {
    for(size_t j=0; j<(size_t)n_box_dim[1]; j++) {
      for(size_t i=0; i<(size_t)n_box_dim[0]; i++) {

        output[mfi].dataPtr()[local_indx] = fac * std::real(a[local_indx]);
        local_indx++;

      }
    }
    }
   // END of parallel
    }
    return output;
    }


// inline void SpectralModifier::apply_func(double (*amp_func)(double), MultiFab& output)
inline void SpectralModifier::apply_func(std::function<double(double)> amp_func, MultiFab& output)
    {

    // Assume for now that nx = ny = nz
    int Ndims[3] = {numb_boxes_dim[2], numb_boxes_dim[1], numb_boxes_dim[0]};
    int     n[3] = {domain.length(2), domain.length(1), domain.length(0)}; // C-order

   // Hardware/Hybrid Cosmology Code (HACC) provides a functionality for performing Discrete Fast Fourier Transforms in AMReX (tailored for cosmological simulations.)
    hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
    hacc::Dfft dfft(d); // creates the dfftiser


    for (MFIter mfi((*input),false); mfi.isValid(); ++mfi) // iterating over the FArrayBox objects (subboxes = Box) that are stored locally on the current MPI rank.
    {
       int gid = mfi.index();

       size_t local_size  = dfft.local_size();


       dfft.makePlans(&a[0],&b[0],&a[0],&b[0]); // designs the fft to be optimal given inputs and outputs dims

       // *******************************************
       // Copy real data from input into real part of a -- no ghost cells and
       // put into C++ ordering (not Fortran)
       // *******************************************
       complex_t zero(0.0, 0.0);
       size_t local_indx = 0;
       for(size_t k=0; k<(size_t)n_box_dim[2]; k++) {
        for(size_t j=0; j<(size_t)n_box_dim[1]; j++) {
         for(size_t i=0; i<(size_t)n_box_dim[0]; i++) {

           complex_t temp((*input)[mfi].dataPtr()[local_indx],0.);
           a[local_indx] = temp;
           local_indx++;

         }
       }
      }

       dfft.forward(&a[0]); // surely this is done in parallel


    //  *******************************************
//  Now apply amp_func to the coefficients of the transform
//  *******************************************
// Warning: DFT such that kx,ky,kz range from 0 to N-1
    local_indx = 0;
    const int *self = dfft.self_kspace(); // rank location in k-space?
    const int *local_ng = dfft.local_ng_kspace();  // local grid dimensions in k-space?
    const int *global_ng = dfft.global_ng();

   // global_i/j/k are the global indices for the k modes.
    for(size_t i=0; i<(size_t)local_ng[0]; i++) {
     size_t global_i = local_ng[0]*self[0] + i;

     for(size_t j=0; j<(size_t)local_ng[1]; j++) {
      size_t global_j = local_ng[1]*self[1] + j;

      for(size_t k=0; k<(size_t)local_ng[2]; k++) {
        size_t global_k = local_ng[2]*self[2] + k;

        if (global_i == 0 && global_j == 0 & global_k == 0) {
           a[local_indx] = 0.0;
        } else {
         // This is k but remember that gradients in discrete space are not that easy. It's the effective k (see Caravano21), not a problem if there is a cutoff to keep global_i,j,k<<global_ng
           double fac = sqrt(-2. * (
                        (cos(2*PI*double(global_i)/double(global_ng[0])) - 1.) +
                        (cos(2*PI*double(global_j)/double(global_ng[1])) - 1.) +
                        (cos(2*PI*double(global_k)/double(global_ng[2])) - 1.) ));
         // dimensionful effective k is fac/dx
           a[local_indx] = a[local_indx] * std::abs(amp_func(fac/h)); // only keep the square root of the power spectrum.
        }
        local_indx++;

      }
     }
    }

    dfft.backward(&a[0]);

    size_t global_size  = dfft.global_size(); // is this Nx*Ny*Nz? yes
  //  double fac = pow(h,2) / global_size; // why is that how it was orignally normalised? because laplacian?
    double fac = pow(sqrt(2*PI)/h,1.5)/sqrt(global_size);

    local_indx = 0;

    for(size_t k=0; k<(size_t)n_box_dim[2]; k++) {
    for(size_t j=0; j<(size_t)n_box_dim[1]; j++) {
      for(size_t i=0; i<(size_t)n_box_dim[0]; i++) {

        output[mfi].dataPtr()[local_indx] = fac * std::real(a[local_indx]);
        // if (output[mfi].dataPtr()[local_indx]!=0.0) amrex::Print() << "filled output:" << a[local_indx] << std::endl;
        local_indx++;

      }
    }
    }
   // END of parallel
    }
}

// Function to fill MultiFab with random noise
inline void SpectralModifier::FillInputWithRandomNoise(std::mt19937& gen) // could use a parallel for there no? no
{
    double stddev = pow(double(n_box_dim[0]*n_box_dim[1]*n_box_dim[2]),-0.5); // Normalisation to get a rayleigh draw once in Fourier space.
    std::normal_distribution<Real>  dis(0.0, stddev); 
    for (MFIter mfi((*input)); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        auto arr = ((*input)[mfi]).array();

        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    arr(i,j,k) = dis(gen);
                }
            }
        }
    }
}


#endif /* SPECTRALMODIFIER_IMPL_HPP_ */

///////////////////////////////
   //  How do you do fft on small domains and then gather for bigger wavelengths?? 
   //  Well The F to C ordering conversion was not the end!!
   //  "SWFFT takes three-dimensional arrays of data distributed across block-structured grids,
   //  and redistributes the data into “pencil” grids in and then, belonging to different MPI processes. 
   // After each pencil conversion, a 1D FFT is performed on the data along the pencil direction using calls to the FFTW 3 library."
   // See https://amrex-codes.github.io/amrex/docs_html/SWFFT.html 
//////////////////////////////

