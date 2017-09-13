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
#include "VELOCITY_FIELD.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VELOCITY_FIELD<TV>::
VELOCITY_FIELD()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> VELOCITY_FIELD<TV>::
~VELOCITY_FIELD()
{
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> bool VELOCITY_FIELD<TV>::
Compute(const GRID<TV>& grid,
        const CUSTOM_PARSE_ARGS<TV>& custom_parse_args,
        const T dt,
        const T time,
        T_ARRAYS_VECTOR& cell_centered_velocity_field)
{
    T result = true;

    if (custom_parse_args.velocity == 0) {
        T coeff = pi / (T)3.15;
        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();

            for (int axis = 1; axis <= TV::dimension; ++axis) {
                if (axis == 1) {
                    cell_centered_velocity_field(cell)(axis) =
                        coeff * (custom_parse_args.velocity_center(axis) - grid.X(cell)(axis+1));
                } else if (axis == 2) {
                    cell_centered_velocity_field(cell)(axis) =
                        coeff * (grid.X(cell)(axis-1) - custom_parse_args.velocity_center(axis));
                }
            }
        }
    } else if (custom_parse_args.velocity == 1) {
        // Single vortex
        T time_period = custom_parse_args.time_period;

        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();

            for (int axis = 1; axis <= TV::dimension; ++axis) {
                if (axis == 1) {
                    cell_centered_velocity_field(cell)(axis) =
                        sqr(sin(pi*grid.X(cell)(axis))) * sin((T)2*pi*grid.X(cell)(axis+1)) * cos(pi*time/time_period);
                } else if (axis == 2) {
                    cell_centered_velocity_field(cell)(axis) =
                        -sqr(sin(pi*grid.X(cell)(axis))) * sin((T)2*pi*grid.X(cell)(axis-1)) * cos(pi*time/time_period);
                }
            }
        }

    } else if (custom_parse_args.velocity == 3) {
        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();

            cell_centered_velocity_field(cell) = custom_parse_args.velocity_vector;
        }
    }

    return result;
}
//#####################################################################
template class VELOCITY_FIELD<VECTOR<float,1> >;
template class VELOCITY_FIELD<VECTOR<float,2> >;
template class VELOCITY_FIELD<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VELOCITY_FIELD<VECTOR<double,1> >;
template class VELOCITY_FIELD<VECTOR<double,2> >;
template class VELOCITY_FIELD<VECTOR<double,3> >;
#endif
