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
#include "SURFACE_INTERPOLATION.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_INTERPOLATION<TV>::
SURFACE_INTERPOLATION()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SURFACE_INTERPOLATION<TV>::
~SURFACE_INTERPOLATION()
{
}
//#####################################################################
// Interpolate_Cell_To_Face_Scalar
//#####################################################################
template<class TV> bool SURFACE_INTERPOLATION<TV>::
Interpolate_Cell_To_Face_Scalar(const GRID<TV>& grid,
        const T_FACE_ARRAYS_SCALAR& face_velocity_flux,
        const T_FACE_ARRAYS_SCALAR& surface_area_field,
        const T_ARRAYS_SCALAR& cell_centered_field_ghost,
        const T_ARRAYS_VECTOR& grad_cell_centered_field_ghost,
        const SurfaceInterpolationScheme surface_interpolation_scheme,
        const T surface_interpolation_scheme_limiter_factor,
        T_FACE_ARRAYS_SCALAR& face_centered_field)
{
    T result = true;

    T lambda_face = (T)0;
    T_FACE_ARRAYS_SCALAR geometric_weights_field(grid, 0, false);
    INTERPOLATION_WEIGHTS<TV> interpolation_weights;

    interpolation_weights.Compute_Geometric_Weights(grid, surface_area_field, geometric_weights_field);

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        //if (surface_interpolation_scheme == UPWIND) {

            //if (face_velocity_flux(axis, face_index) >= 0) {
                //face_centered_field(axis, face_index) = cell_centered_field_ghost(first_cell);
            //} else {
                //face_centered_field(axis, face_index) = cell_centered_field_ghost(second_cell);
            //}

        //} else {

            lambda_face = interpolation_weights.Compute_Interpolation_Weight(grid,
                    geometric_weights_field(axis, face_index),
                    face_velocity_flux(axis, face_index),
                    cell_centered_field_ghost(first_cell),
                    cell_centered_field_ghost(second_cell),
                    grad_cell_centered_field_ghost(first_cell),
                    grad_cell_centered_field_ghost(second_cell),
                    grid.X(second_cell) - grid.X(first_cell),
                    surface_interpolation_scheme,
                    surface_interpolation_scheme_limiter_factor);

            face_centered_field(axis, face_index) =
                lambda_face * (cell_centered_field_ghost(first_cell) - cell_centered_field_ghost(second_cell)) +
                cell_centered_field_ghost(second_cell);

        //}
    }

    return result;
}
//#####################################################################
template class SURFACE_INTERPOLATION<VECTOR<float,1> >;
template class SURFACE_INTERPOLATION<VECTOR<float,2> >;
template class SURFACE_INTERPOLATION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SURFACE_INTERPOLATION<VECTOR<double,1> >;
template class SURFACE_INTERPOLATION<VECTOR<double,2> >;
template class SURFACE_INTERPOLATION<VECTOR<double,3> >;
#endif
