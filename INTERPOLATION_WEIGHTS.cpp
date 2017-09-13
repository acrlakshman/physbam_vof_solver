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
#include "INTERPOLATION_WEIGHTS.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERPOLATION_WEIGHTS<TV>::
INTERPOLATION_WEIGHTS()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERPOLATION_WEIGHTS<TV>::
~INTERPOLATION_WEIGHTS()
{
}
//#####################################################################
// Compute_Geometric_Weights
//#####################################################################
template<class TV> bool INTERPOLATION_WEIGHTS<TV>::
Compute_Geometric_Weights(const GRID<TV>& grid, const T_FACE_ARRAYS_SCALAR& surface_area_field,
        T_FACE_ARRAYS_SCALAR& geometric_weights_field)
{
    bool status = true;

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {

        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        // NOTE: sign of s_f and other vectors do not matter here, since
        // we need magnitudes only
        TV x_face = (T)0.5 * (grid.X(first_cell) + grid.X(second_cell));
        TV x_2f = grid.X(second_cell) - x_face;
        TV x_1f = x_face - grid.X(first_cell);
        TV s_f = TV();
        s_f(axis) = surface_area_field(axis, face_index);

        geometric_weights_field(axis, face_index) = abs(Dot_Product<TV>(s_f, x_2f));
        geometric_weights_field(axis, face_index) /= (abs(Dot_Product<TV>(s_f, x_2f)) + abs(Dot_Product<TV>(s_f, x_1f)));

    }

    return status;
}
//#####################################################################
// Compute_Interpolation_Weight
//#####################################################################
template<class TV> typename TV::SCALAR INTERPOLATION_WEIGHTS<TV>::
Compute_Interpolation_Weight(const GRID<TV>& grid,
        const T geometric_weight,
        const T face_velocity_flux,
        const T alpha1_first_cell,
        const T alpha1_second_cell,
        const TV grad_alpha1_first_cell,
        const TV grad_alpha1_second_cell,
        const TV cell_distance_second_first,
        const SurfaceInterpolationScheme surface_interpolation_scheme,
        const T surface_interpolation_scheme_limiter_factor)
{
    T interpolation_weight = (T)0;

    if (face_velocity_flux >= (T)0) {
        interpolation_weight = (T)1;
    } else {
        interpolation_weight = (T)0;
    }

    switch (surface_interpolation_scheme) {
        case UPWIND:
            break;

        case LINEAR:
            interpolation_weight = geometric_weight;
            break;

        default: {
                     T limiter = LIMITER_FUNCTIONS<TV>::Limiter(geometric_weight,
                             face_velocity_flux,
                             alpha1_first_cell,
                             alpha1_second_cell,
                             grad_alpha1_first_cell,
                             grad_alpha1_second_cell,
                             cell_distance_second_first,
                             surface_interpolation_scheme,
                             surface_interpolation_scheme_limiter_factor);

                     interpolation_weight = (limiter * geometric_weight) +
                         (((T)1 - limiter) * pos0(face_velocity_flux));
                 }
                 break;
    }

    return interpolation_weight;
}
//#####################################################################
template class INTERPOLATION_WEIGHTS<VECTOR<float,1> >;
template class INTERPOLATION_WEIGHTS<VECTOR<float,2> >;
template class INTERPOLATION_WEIGHTS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INTERPOLATION_WEIGHTS<VECTOR<double,1> >;
template class INTERPOLATION_WEIGHTS<VECTOR<double,2> >;
template class INTERPOLATION_WEIGHTS<VECTOR<double,3> >;
#endif
