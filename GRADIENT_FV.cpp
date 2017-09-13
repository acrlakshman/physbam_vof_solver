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
#include "GRADIENT_FV.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRADIENT_FV<TV>::
GRADIENT_FV()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GRADIENT_FV<TV>::
~GRADIENT_FV()
{
}
//#####################################################################
// Compute_For_Scalar_Field (one face centered scalar field)
//#####################################################################
template<class TV> bool GRADIENT_FV<TV>::
Compute_For_Scalar_Field(const GRID<TV>& grid,
        const T_FACE_ARRAYS_SCALAR& face_centered_scalar_field,
        const T_FACE_ARRAYS_SCALAR& surface_area_field,
        T_ARRAYS_VECTOR& gradient_value_field)
{
    bool result = true;

    T one_by_volume = (T)1. / grid.DX().Product();

    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();
        gradient_value_field(cell) = TV();

        for (int axis = 1; axis <= TV::dimension; ++axis) {
            TV face_centered_vector = TV();
            face_centered_vector(axis) = (face_centered_scalar_field(axis, iterator.Second_Face_Index(axis)) *
                surface_area_field(axis, iterator.Second_Face_Index(axis)));
            gradient_value_field(cell) += face_centered_vector;

            face_centered_vector(axis) = (face_centered_scalar_field(axis, iterator.First_Face_Index(axis)) *
                    surface_area_field(axis, iterator.First_Face_Index(axis)));
            gradient_value_field(cell) -= face_centered_vector;
        }
        gradient_value_field(cell) *= one_by_volume;
    }

    return result;
}
//#####################################################################
template class GRADIENT_FV<VECTOR<float,1> >;
template class GRADIENT_FV<VECTOR<float,2> >;
template class GRADIENT_FV<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRADIENT_FV<VECTOR<double,1> >;
template class GRADIENT_FV<VECTOR<double,2> >;
template class GRADIENT_FV<VECTOR<double,3> >;
#endif
