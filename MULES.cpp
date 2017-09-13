//#####################################################################
// Copyright 2017, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include "MULES.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MULES<TV>::
MULES()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MULES<TV>::
~MULES()
{
}
//#####################################################################
// Limiter
//#####################################################################
template<class TV> bool MULES<TV>::
Limiter(const GRID<TV>& grid,
        const MPI_UNIFORM_GRID<T_GRID> *mpi_grid,
        const T rDeltaT,
        const T_ARRAYS_SCALAR& alpha_field_ghost,
        const T_ARRAYS_SCALAR& alpha_previous_field_ghost,
        const T_FACE_ARRAYS_SCALAR& face_velocity_flux_field,
        const T_FACE_ARRAYS_SCALAR& face_alpha_upwind_flux_field,
        //const T_FACE_ARRAYS_SCALAR& face_alpha_flux_field,
        const T_FACE_ARRAYS_SCALAR& face_alpha_flux_corr_field,
        const T_ARRAYS_SCALAR& Sp,
        const T_ARRAYS_SCALAR& Su,
        const T alpha_max,
        const T alpha_min,
        T_BOUNDARY_SCALAR *boundary,
        T_FACE_ARRAYS_SCALAR& lambda
       )
{
    T result = true;

    int n_limiter_iter = 3;
    T extrema_coeff = (T)0;
    T smooth_limiter = (T)0;

    // Initialize lambda to 1
    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        //TV_INT first_cell = iterator.First_Cell_Index();
        //TV_INT second_cell = iterator.Second_Cell_Index();

        lambda(axis, face_index) = (T)1;
    }

    int num_ghost_cells = 1;
    T_ARRAYS_SCALAR psiMaxn_ghost(grid.Domain_Indices(num_ghost_cells));
    T_ARRAYS_SCALAR psiMinn_ghost(grid.Domain_Indices(num_ghost_cells));
    T_ARRAYS_SCALAR sumPhiBD_ghost(grid.Domain_Indices(num_ghost_cells));
    T_ARRAYS_SCALAR sumPhip_ghost(grid.Domain_Indices(num_ghost_cells));
    T_ARRAYS_SCALAR mSumPhim_ghost(grid.Domain_Indices(num_ghost_cells));

    T_ARRAYS_SCALAR sumlPhip_ghost(grid.Domain_Indices(num_ghost_cells));
    T_ARRAYS_SCALAR mSumlPhim_ghost(grid.Domain_Indices(num_ghost_cells));

    for (CELL_ITERATOR iterator(grid, num_ghost_cells); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();
        psiMaxn_ghost(cell) = alpha_min;
        psiMinn_ghost(cell) = alpha_max;
        sumPhiBD_ghost(cell) = (T)0;
        sumPhip_ghost(cell) = (T)0;
        mSumPhim_ghost(cell) = (T)0;

        sumlPhip_ghost(cell) = (T)0;
        mSumlPhim_ghost(cell) = (T)0;
    }

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        psiMaxn_ghost(first_cell) = max(psiMaxn_ghost(first_cell), alpha_field_ghost(second_cell));
        psiMinn_ghost(first_cell) = min(psiMinn_ghost(first_cell), alpha_field_ghost(second_cell));

        psiMaxn_ghost(second_cell) = max(psiMaxn_ghost(second_cell), alpha_field_ghost(first_cell));
        psiMinn_ghost(second_cell) = min(psiMinn_ghost(second_cell), alpha_field_ghost(first_cell));

        sumPhiBD_ghost(first_cell) += face_alpha_upwind_flux_field(axis, face_index);
        sumPhiBD_ghost(second_cell) -= face_alpha_upwind_flux_field(axis, face_index);

        if (face_alpha_flux_corr_field(axis, face_index) > (T)0) {
            sumPhip_ghost(first_cell) += face_alpha_flux_corr_field(axis, face_index);
            mSumPhim_ghost(second_cell) += face_alpha_flux_corr_field(axis, face_index);
        } else {
            mSumPhim_ghost(first_cell) -= face_alpha_flux_corr_field(axis, face_index);
            sumPhip_ghost(second_cell) -= face_alpha_flux_corr_field(axis, face_index);
        }
    }

    for (CELL_ITERATOR iterator(grid, num_ghost_cells); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();

        psiMaxn_ghost(cell) = min(psiMaxn_ghost(cell) + extrema_coeff * (alpha_max - alpha_min), alpha_max);
        psiMinn_ghost(cell) = max(psiMinn_ghost(cell) - extrema_coeff * (alpha_max - alpha_min), alpha_min);
    }

    if (smooth_limiter > SMALL_NUMBER) {
        for (CELL_ITERATOR iterator(grid, num_ghost_cells); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();

            psiMaxn_ghost(cell) = min(smooth_limiter * alpha_field_ghost(cell) +
                    ((T)1 - smooth_limiter) * psiMaxn_ghost(cell), alpha_max);
            psiMinn_ghost(cell) = max(smooth_limiter * alpha_field_ghost(cell) +
                    ((T)1 - smooth_limiter) * psiMinn_ghost(cell), alpha_min);
        }
    }

    T vol = grid.DX().Product();
    for (CELL_ITERATOR iterator(grid, num_ghost_cells); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();

        psiMaxn_ghost(cell) = vol * ( (rDeltaT - Sp(cell)) * psiMaxn_ghost(cell) -
                Su(cell) - (rDeltaT * alpha_previous_field_ghost(cell)) ) +
            sumPhiBD_ghost(cell);

        psiMinn_ghost(cell) = vol * ( Su(cell) - (rDeltaT - Sp(cell)) * psiMinn_ghost(cell) +
                (rDeltaT * alpha_previous_field_ghost(cell)) ) -
            sumPhiBD_ghost(cell);
    }

    for (int i_limiter_count = 0; i_limiter_count < n_limiter_iter; ++i_limiter_count) {

        for (CELL_ITERATOR iterator(grid, num_ghost_cells); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();

            sumlPhip_ghost(cell) = (T)0;
            mSumlPhim_ghost(cell) = (T)0;
        }

        T lambda_phi_corr = (T)0;

        for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            FACE_INDEX<TV::dimension> face = iterator.Full_Index();
            int axis = face.axis;
            TV_INT face_index = face.index;
            TV_INT first_cell = iterator.First_Cell_Index();
            TV_INT second_cell = iterator.Second_Cell_Index();

            lambda_phi_corr = lambda(axis, face_index) * face_alpha_flux_corr_field(axis, face_index);

            if (lambda_phi_corr > (T)0) {
                sumlPhip_ghost(first_cell) += lambda_phi_corr;
                mSumlPhim_ghost(second_cell) += lambda_phi_corr;
            } else {
                mSumlPhim_ghost(first_cell) -= lambda_phi_corr;
                sumlPhip_ghost(second_cell) -= lambda_phi_corr;
            }
        }

        for (CELL_ITERATOR iterator(grid, num_ghost_cells); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();

            sumlPhip_ghost(cell) = max(min(
                        (sumlPhip_ghost(cell) + psiMaxn_ghost(cell)) /
                        (mSumPhim_ghost(cell) + (T)ROOT_V_SMALL), (T)1.), (T)0.);

            mSumlPhim_ghost(cell) = max(min(
                        (mSumlPhim_ghost(cell) + psiMinn_ghost(cell)) /
                        (sumPhip_ghost(cell) + (T)ROOT_V_SMALL), (T)1.), (T)0.);
        }

        for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            FACE_INDEX<TV::dimension> face = iterator.Full_Index();
            int axis = face.axis;
            TV_INT face_index = face.index;
            TV_INT first_cell = iterator.First_Cell_Index();
            TV_INT second_cell = iterator.Second_Cell_Index();

            if (face_alpha_flux_corr_field(axis, face_index) > (T)0) {
                lambda(axis, face_index) = min(lambda(axis, face_index),
                        min(mSumlPhim_ghost(first_cell), sumlPhip_ghost(second_cell)));
            } else {
                lambda(axis, face_index) = min(lambda(axis, face_index),
                        min(sumlPhip_ghost(first_cell), mSumlPhim_ghost(second_cell)));
            }
        }

    }

    //for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        //FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        //int axis = face.axis;
        //TV_INT face_index = face.index;
        //TV_INT first_cell = iterator.First_Cell_Index();
        //TV_INT second_cell = iterator.Second_Cell_Index();

        //TV surface_area_vector = TV();
        //surface_area_vector(axis) = surface_area_field(axis, face_index);

        //TV face_velocity_vector = TV();
        //face_velocity_vector(axis) = face_velocities(axis, face_index);

        //face_velocity_flux(axis, face_index) = Dot_Product<TV>(face_velocity_vector, surface_area_vector);
    //}

    return result;
}
//#####################################################################
template class MULES<VECTOR<float,1> >;
template class MULES<VECTOR<float,2> >;
template class MULES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MULES<VECTOR<double,1> >;
template class MULES<VECTOR<double,2> >;
template class MULES<VECTOR<double,3> >;
#endif
