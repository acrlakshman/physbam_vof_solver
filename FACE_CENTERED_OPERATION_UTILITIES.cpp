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
#include "FACE_CENTERED_OPERATION_UTILITIES.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FACE_CENTERED_OPERATION_UTILITIES<TV>::
FACE_CENTERED_OPERATION_UTILITIES()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FACE_CENTERED_OPERATION_UTILITIES<TV>::
~FACE_CENTERED_OPERATION_UTILITIES()
{
}
//#####################################################################
// Vector_Interpolate_Cell_To_Face_Linear
//#####################################################################
template<class TV> bool FACE_CENTERED_OPERATION_UTILITIES<TV>::
Vector_Interpolate_Cell_To_Face_Linear(const GRID<TV>& grid,
        const T_ARRAYS_VECTOR& vector_field_ghost,
        T_FACE_ARRAYS_SCALAR& face_field_for_vector)
{
    bool result = true;

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        face_field_for_vector(axis, face_index) =
            (T)0.5 * (vector_field_ghost(first_cell)(axis) + vector_field_ghost(second_cell)(axis));
    }

    return result;
}
//#####################################################################
// Dot_Product_Face_Fields
//#####################################################################
template<class TV> bool FACE_CENTERED_OPERATION_UTILITIES<TV>::
Dot_Product_Face_Fields(const GRID<TV>& grid, const T_FACE_ARRAYS_SCALAR& face_field_one,
        const T_FACE_ARRAYS_SCALAR& face_field_two, T_FACE_ARRAYS_SCALAR& face_field_dot_product)
{
    bool result = true;

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        TV face_vector_one = TV();
        TV face_vector_two = TV();

        face_vector_one(axis) = face_field_one(axis, face_index);
        face_vector_two(axis) = face_field_two(axis, face_index);

        face_field_dot_product(axis, face_index) =
            Dot_Product<TV>(face_vector_one, face_vector_two);
    }

    return result;
}
//#####################################################################
template class FACE_CENTERED_OPERATION_UTILITIES<VECTOR<float,1> >;
template class FACE_CENTERED_OPERATION_UTILITIES<VECTOR<float,2> >;
template class FACE_CENTERED_OPERATION_UTILITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FACE_CENTERED_OPERATION_UTILITIES<VECTOR<double,1> >;
template class FACE_CENTERED_OPERATION_UTILITIES<VECTOR<double,2> >;
template class FACE_CENTERED_OPERATION_UTILITIES<VECTOR<double,3> >;
#endif
