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
#include "INITIALIZE_ALPHA1.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INITIALIZE_ALPHA1<TV>::
INITIALIZE_ALPHA1()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INITIALIZE_ALPHA1<TV>::
~INITIALIZE_ALPHA1()
{
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> bool INITIALIZE_ALPHA1<TV>::
Initialize(const GRID<TV>& grid,
        const CUSTOM_PARSE_ARGS<TV>& custom_parse_args,
        T_ARRAYS_SCALAR& alpha1_field)
{
    T result = true;

    if (custom_parse_args.object == 1) {
        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();

            T distance = (T)0;
            for (int axis = 1; axis <= TV::dimension; ++axis) {
                distance += sqr(grid.X(cell)(axis) - custom_parse_args.object_center(axis));
            }
            distance = sqrt(distance);
            distance -= custom_parse_args.object_radius;

            alpha1_field(cell) = (distance <= (T)0) ? (T)1 : (T)0;
        }
    }

    return result;
}
//#####################################################################
template class INITIALIZE_ALPHA1<VECTOR<float,1> >;
template class INITIALIZE_ALPHA1<VECTOR<float,2> >;
template class INITIALIZE_ALPHA1<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INITIALIZE_ALPHA1<VECTOR<double,1> >;
template class INITIALIZE_ALPHA1<VECTOR<double,2> >;
template class INITIALIZE_ALPHA1<VECTOR<double,3> >;
#endif
