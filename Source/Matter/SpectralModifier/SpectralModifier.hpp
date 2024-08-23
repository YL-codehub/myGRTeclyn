// Credits to Yoann L. Launay, adapted from HACC/SWFFT or amrex-tutorials

#ifndef SPECTRALMODIFIER_HPP_
#define SPECTRALMODIFIER_HPP_

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
    SpectralModifier (const MultiFab& input,
                Geometry& geom,
                int verbose);

    MultiFab apply_func (double (*amp_func)(double));
    void apply_array ();

private:

    const MultiFab& input;
    Geometry geom;
    int verbose = 2;
    const BoxArray& ba;
    DistributionMapping dmap;
    IntVect n_box_dim; // n for each subbox (assumed the same everywhere) for each dim  // SIZE OF EACH SUBGRID (one per rank). ba[0] gives the first subbox. Same nx,ny,nz on all ranks // size of BoxArray = number of subboxes
    Box domain;
    IntVect numb_boxes_dim; // number of subboxes per direction //length of Box object = associated nx/ny/nz
    int numb_boxes; // total number of subboxes
    Vector<int> rank_mapping;
    Real h;

    void remap ();
};


#include "SpectralModifier.impl.hpp"
#endif /* SPECTRALMODIFIER_HPP_ */