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
#include "SURFACE_NORMAL_GRADIENT.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_NORMAL_GRADIENT<TV>::
SURFACE_NORMAL_GRADIENT()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SURFACE_NORMAL_GRADIENT<TV>::
~SURFACE_NORMAL_GRADIENT()
{
}
//#####################################################################
// Compute_For_Scalar_Field (one cell centered scalar field)
//#####################################################################
template<class TV> bool SURFACE_NORMAL_GRADIENT<TV>::
Compute_For_Scalar_Field(const GRID<TV>& grid,
            const T_ARRAYS_SCALAR& cell_centered_scalar_field_ghost,
            T_FACE_ARRAYS_SCALAR& face_centered_scalar_field)
{
    bool result = true;

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        T one_by_distance = grid.X(second_cell)(axis) - grid.X(first_cell)(axis);
        one_by_distance = (T)1. / one_by_distance;

        face_centered_scalar_field(axis, face_index) = (cell_centered_scalar_field_ghost(second_cell) -
            cell_centered_scalar_field_ghost(first_cell)) * one_by_distance;
    }

    return result;
}
//#####################################################################
template class SURFACE_NORMAL_GRADIENT<VECTOR<float,1> >;
template class SURFACE_NORMAL_GRADIENT<VECTOR<float,2> >;
template class SURFACE_NORMAL_GRADIENT<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SURFACE_NORMAL_GRADIENT<VECTOR<double,1> >;
template class SURFACE_NORMAL_GRADIENT<VECTOR<double,2> >;
template class SURFACE_NORMAL_GRADIENT<VECTOR<double,3> >;
#endif
