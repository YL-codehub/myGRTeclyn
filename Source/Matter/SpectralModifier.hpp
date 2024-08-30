// Credits to Yoann L. Launay, adapted from HACC/SWFFT or amrex-tutorials

#ifndef SPECTRALMODIFIER_HPP
#define SPECTRALMODIFIER_HPP

#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

using namespace amrex;

class SpectralModifier
{
public:
    inline SpectralModifier (MultiFab& input,
                Geometry& geom,
                int verbose);
                
    inline SpectralModifier ();
    inline MultiFab apply_func (double (*amp_func)(double));
    inline void apply_func(double (*amp_func)(double), MultiFab& output);
    inline void apply_array ();
    inline void FillInputWithRandomNoise(std::mt19937 gen);

private:
    // pointers are sometimes useful because they can be initialized as nullptr
    MultiFab* input;
    Geometry geom;
    int verbose = 2;
    const BoxArray* ba;
    DistributionMapping dmap;
    IntVect n_box_dim; // n for each subbox (assumed the same everywhere) for each dim  // SIZE OF EACH SUBGRID (one per rank). ba[0] gives the first subbox. Same nx,ny,nz on all ranks // size of BoxArray = number of subboxes
    Box domain;
    IntVect numb_boxes_dim; // number of subboxes per direction //length of Box object = associated nx/ny/nz
    int numb_boxes; // total number of subboxes
    Vector<int> rank_mapping; // to be passed from amrex rank mapping to dfft's
    Real h;
};

#include "SpectralModifier.impl.hpp"
#endif /* SPECTRALMODIFIER_HPP */