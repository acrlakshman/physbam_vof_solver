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
#include "ALPHA1_UTILITIES.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ALPHA1_UTILITIES<TV>::
ALPHA1_UTILITIES()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ALPHA1_UTILITIES<TV>::
~ALPHA1_UTILITIES()
{
}
//#####################################################################
// Compute_Curvature
//#####################################################################
template<class TV> bool ALPHA1_UTILITIES<TV>::
Compute_Curvature(const GRID<TV>& grid,
        const T_ARRAYS_VECTOR& grad_alpha1_field_ghost,
        const T_FACE_ARRAYS_SCALAR& surface_area_field,
        T_ARRAYS_SCALAR& curvature_field
        )
{
    bool result = true;

    // TODO: HERE
    // Procedure: kappa = -divergence(normal)
    // = -(1/dV) sum(n_f \cdot s_f)
    // = -(1/dV) sum(\nabla \alpha1_f \cdot s_f / (|\nabla \alpha1_f| + SMALL))
    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
       const TV_INT cell = iterator.Cell_Index();

       curvature_field(cell) = (T)0;
       T deltaN_ = (T)SMALL_NUMBER;

       TV face_normal = TV();
       TV face_surface_area = TV();
       T face_normal_dot_face_surface_area = (T)0;
       T one_by_volume = (T)1. / grid.DX().Product();

       for (int axis = 1; axis <= TV::dimension; ++axis) {
          // NOTE: Linear interpolation
          // Second face
          face_normal = (T)0.5 * (grad_alpha1_field_ghost(cell) +
                grad_alpha1_field_ghost(cell + TV_INT::Axis_Vector(axis)));
          face_normal /= (face_normal.Magnitude() + deltaN_);

          face_surface_area = TV();
          face_surface_area(axis) = surface_area_field(axis, iterator.Second_Face_Index(axis));

          face_normal_dot_face_surface_area =
             Dot_Product<TV>(face_normal, face_surface_area);

          curvature_field(cell) += face_normal_dot_face_surface_area;

          // First face
          face_normal = (T)0.5 * (grad_alpha1_field_ghost(cell) +
                grad_alpha1_field_ghost(cell - TV_INT::Axis_Vector(axis)));
          face_normal /= (face_normal.Magnitude() + deltaN_);

          face_surface_area = TV();
          face_surface_area(axis) = surface_area_field(axis, iterator.First_Face_Index(axis));

          face_normal_dot_face_surface_area =
             Dot_Product<TV>(face_normal, face_surface_area);

          curvature_field(cell) -= face_normal_dot_face_surface_area;
       }
       curvature_field(cell) *= (-one_by_volume);
    }

    return result;
}
//#####################################################################
template class ALPHA1_UTILITIES<VECTOR<float,1> >;
template class ALPHA1_UTILITIES<VECTOR<float,2> >;
template class ALPHA1_UTILITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ALPHA1_UTILITIES<VECTOR<double,1> >;
template class ALPHA1_UTILITIES<VECTOR<double,2> >;
template class ALPHA1_UTILITIES<VECTOR<double,3> >;
#endif
