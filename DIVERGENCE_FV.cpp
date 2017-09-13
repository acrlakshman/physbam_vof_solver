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
#include "DIVERGENCE_FV.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DIVERGENCE_FV<TV>::
DIVERGENCE_FV()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DIVERGENCE_FV<TV>::
~DIVERGENCE_FV()
{
}
//#####################################################################
// Compute_For_Scalar_Field (single surface field)
//#####################################################################
template<class TV> bool DIVERGENCE_FV<TV>::
Compute_For_Scalar_Field(const GRID<TV>& grid,
        const T_FACE_ARRAYS_SCALAR& face_scalar_field,
        T_ARRAYS_SCALAR& divergence_value_field
       )
{
    bool result = true;

    T one_by_volume = (T)1. / grid.DX().Product();

    //for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        //FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        //int axis = face.axis;
        //TV_INT face_index = face.index;
        //TV_INT first_cell = iterator.First_Cell_Index();
        //TV_INT second_cell = iterator.Second_Cell_Index();

        //divergence_value_field(first_cell) += face_scalar_field(axis, face_index);
        //divergence_value_field(second_cell) -= face_scalar_field(axis, face_index);
    //}
    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();
        divergence_value_field(cell) = (T)0;

        for (int axis = 1; axis <= TV::dimension; ++axis) {
            divergence_value_field(cell) += face_scalar_field(axis, iterator.Second_Face_Index(axis));
            divergence_value_field(cell) -= face_scalar_field(axis, iterator.First_Face_Index(axis));
        }
        divergence_value_field(cell) *= one_by_volume;
    }

    return result;
}
//#####################################################################
// Compute_For_Scalar_Field (two surface fields)
//#####################################################################
template<class TV> bool DIVERGENCE_FV<TV>::
Compute_For_Scalar_Field(const GRID<TV>& grid,
        const T_FACE_ARRAYS_SCALAR& face_scalar_field_one,
        const T_FACE_ARRAYS_SCALAR& face_scalar_field_two,
        T_ARRAYS_SCALAR& divergence_value_field
       )
{
    bool result = true;

    T one_by_volume = (T)1. / grid.DX().Product();

    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();
        divergence_value_field(cell) = (T)0;

        for (int axis = 1; axis <= TV::dimension; ++axis) {
            divergence_value_field(cell) += (face_scalar_field_one(axis, iterator.Second_Face_Index(axis)) *
                    face_scalar_field_two(axis, iterator.Second_Face_Index(axis)));
            divergence_value_field(cell) -= (face_scalar_field_one(axis, iterator.First_Face_Index(axis)) *
                    face_scalar_field_two(axis, iterator.First_Face_Index(axis)));
        }
        divergence_value_field(cell) *= one_by_volume;
    }

    return result;
}
//#####################################################################
// Compute_For_Vector_Field (one face vector field)
//#####################################################################
template<class TV> bool DIVERGENCE_FV<TV>::
Compute_For_Vector_Field(const GRID<TV>& grid,
        const T_FACE_ARRAYS_SCALAR& face_vector_field,
        const T_FACE_ARRAYS_SCALAR& surface_area_field,
        T_ARRAYS_SCALAR& divergence_value_field
       )
{
    bool result = true;

    T_FACE_ARRAYS_SCALAR face_scalar_field_product(grid, 0, false);

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        TV surface_area_vector = TV();
        surface_area_vector(axis) = surface_area_field(axis, face_index);

        TV face_vector = TV();
        face_vector(axis) = face_vector_field(axis, face_index);

        face_scalar_field_product(axis, face_index) = Dot_Product<TV>(face_vector, surface_area_vector);
    }

    result = Compute_For_Scalar_Field(grid, face_scalar_field_product, divergence_value_field);

    return result;
}
//#####################################################################
// Compute_For_Vector_Field (one surface field, one face vector field)
//#####################################################################
template<class TV> bool DIVERGENCE_FV<TV>::
Compute_For_Vector_Field(const GRID<TV>& grid,
        const T_FACE_ARRAYS_SCALAR& face_scalar_field,
        const T_FACE_ARRAYS_SCALAR& face_vector_field,
        const T_FACE_ARRAYS_SCALAR& surface_area_field,
        T_ARRAYS_SCALAR& divergence_value_field
       )
{
    bool result = true;

    T_FACE_ARRAYS_SCALAR face_scalar_field_product(grid, 0, false);

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        TV surface_area_vector = TV();
        surface_area_vector(axis) = surface_area_field(axis, face_index);

        TV face_vector = TV();
        face_vector(axis) = face_vector_field(axis, face_index);

        face_scalar_field_product(axis, face_index) = (face_scalar_field(axis, face_index) *
            Dot_Product<TV>(face_vector, surface_area_vector));
    }

    result = Compute_For_Scalar_Field(grid, face_scalar_field_product, divergence_value_field);

    return result;
}
//#####################################################################
template class DIVERGENCE_FV<VECTOR<float,1> >;
template class DIVERGENCE_FV<VECTOR<float,2> >;
template class DIVERGENCE_FV<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DIVERGENCE_FV<VECTOR<double,1> >;
template class DIVERGENCE_FV<VECTOR<double,2> >;
template class DIVERGENCE_FV<VECTOR<double,3> >;
#endif
