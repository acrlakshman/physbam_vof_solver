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
#include "CELL_CENTERED_OPERATION_UTILITIES.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CELL_CENTERED_OPERATION_UTILITIES<TV>::
CELL_CENTERED_OPERATION_UTILITIES()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> CELL_CENTERED_OPERATION_UTILITIES<TV>::
~CELL_CENTERED_OPERATION_UTILITIES()
{
}
//#####################################################################
// Normalize_Cell_Centered_Vector_Field
//#####################################################################
template<class TV> bool CELL_CENTERED_OPERATION_UTILITIES<TV>::
Normalize_Cell_Centered_Vector_Field(const GRID<TV>& grid,
        const T_ARRAYS_VECTOR& cell_centered_vector_field,
        T_ARRAYS_VECTOR& cell_centered_normal_field)
{
    bool result = true;

    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const TV_INT cell = iterator.Cell_Index();
        cell_centered_normal_field(cell) = cell_centered_vector_field(cell).Normalized();
    }

    return result;
}
//#####################################################################
// Compute_Scalar_Gradient
//#####################################################################
template<class TV> bool CELL_CENTERED_OPERATION_UTILITIES<TV>::
Compute_Scalar_Gradient(const GRID<TV>& grid,
        const MPI_UNIFORM_GRID<T_GRID> *mpi_grid,
        const T_ARRAYS_SCALAR& cell_centered_scalar_field,
        const T_ARRAYS_SCALAR& cell_centered_scalar_field_ghost,
        const short gradient_approximation_order,
        T_BOUNDARY_SCALAR *boundary,
        T_ARRAYS_VECTOR& grad_cell_centered_scalar_field)
{
    // NOTE: Correct below ASSERT if LHS has ghost values
    PHYSBAM_ASSERT(grad_cell_centered_scalar_field.Size()==cell_centered_scalar_field.Size());

    bool result = true;

    TV_INT cells = grid.Counts();

    //T_ARRAYS cell_centered_scalar_field_ghost(grid.Domain_Indices(3));
    //T_ARRAYS::Put(cell_centered_scalar_field,cell_centered_scalar_field_ghost);
    //if(mpi_grid) boundary->Fill_Ghost_Cells(grid,cell_centered_scalar_field,cell_centered_scalar_field_ghost,(T)0.,(T)0.,3);

    if (gradient_approximation_order == 2) {

        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT index = iterator.Cell_Index();

            for (int axis = 1; axis <= TV::dimension; ++axis) {

                if (!cell_centered_scalar_field.Valid_Index(index-TV_INT::Axis_Vector(axis)) && (!mpi_grid || (mpi_grid && !mpi_grid->Neighbor(axis,1)))) {
                    grad_cell_centered_scalar_field(index)(axis) =
                        (cell_centered_scalar_field_ghost(index+TV_INT::Axis_Vector(axis)) -
                         cell_centered_scalar_field_ghost(index))*grid.One_Over_DX()(axis);
                } else if (!cell_centered_scalar_field.Valid_Index(index+TV_INT::Axis_Vector(axis)) && (!mpi_grid || (mpi_grid && !mpi_grid->Neighbor(axis,2)))) {
                    grad_cell_centered_scalar_field(index)(axis) =
                        (cell_centered_scalar_field_ghost(index) -
                         cell_centered_scalar_field_ghost(index-TV_INT::Axis_Vector(axis)))*grid.One_Over_DX()(axis);
                } else {
                    grad_cell_centered_scalar_field(index)(axis) =
                        (cell_centered_scalar_field_ghost(index+TV_INT::Axis_Vector(axis)) -
                         cell_centered_scalar_field_ghost(index-TV_INT::Axis_Vector(axis)))*(T)0.5*grid.One_Over_DX()(axis);
                }
            }
        }

    } else if (gradient_approximation_order == 3) {

        T division_factor = (T)one_sixth;

        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT index = iterator.Cell_Index();

            for (int axis = 1; axis <= TV::dimension; ++axis) {

                if ((index(axis)==1) && (!mpi_grid || (mpi_grid && !mpi_grid->Neighbor(axis,1)))) {
                    grad_cell_centered_scalar_field(index)(axis) =
                        ((-(T)11*cell_centered_scalar_field_ghost(index)) +
                         ((T)18*cell_centered_scalar_field_ghost(index+TV_INT::Axis_Vector(axis))) -
                         ((T)9*cell_centered_scalar_field_ghost(index+TV_INT::Axis_Vector(axis)*(T)2)) +
                         ((T)2*cell_centered_scalar_field_ghost(index+TV_INT::Axis_Vector(axis)*(T)3))) *
                        division_factor*grid.One_Over_DX()(axis);
                } else if ((index(axis) == 2) && (!mpi_grid || (mpi_grid && !mpi_grid->Neighbor(axis,1)))) {
                    grad_cell_centered_scalar_field(index)(axis) =
                        ((-(T)3*cell_centered_scalar_field_ghost(index)) +
                         ((T)6*cell_centered_scalar_field_ghost(index+TV_INT::Axis_Vector(axis))) -
                         (cell_centered_scalar_field_ghost(index+TV_INT::Axis_Vector(axis)*(T)2)) -
                         ((T)2*cell_centered_scalar_field_ghost(index-TV_INT::Axis_Vector(axis)))) *
                        division_factor*grid.One_Over_DX()(axis);
                } else if ((index(axis) == cells(axis)) && (!mpi_grid || (mpi_grid && !mpi_grid->Neighbor(axis,2)))) {
                    grad_cell_centered_scalar_field(index)(axis) =
                        (((T)11*cell_centered_scalar_field_ghost(index)) -
                         ((T)18*cell_centered_scalar_field_ghost(index-TV_INT::Axis_Vector(axis))) +
                         ((T)9*cell_centered_scalar_field_ghost(index-TV_INT::Axis_Vector(axis)*(T)2)) -
                         ((T)2*cell_centered_scalar_field_ghost(index-TV_INT::Axis_Vector(axis)*(T)3))) *
                        division_factor*grid.One_Over_DX()(axis);
                } else {
                    grad_cell_centered_scalar_field(index)(axis) =
                        (((T)3*cell_centered_scalar_field_ghost(index)) +
                         ((T)2*cell_centered_scalar_field_ghost(index+TV_INT::Axis_Vector(axis))) -
                         ((T)6*cell_centered_scalar_field_ghost(index-TV_INT::Axis_Vector(axis))) +
                         (cell_centered_scalar_field_ghost(index-TV_INT::Axis_Vector(axis)*(T)2))) *
                        division_factor*grid.One_Over_DX()(axis);
                }
            }
        }
    }

    return result;
}
//#####################################################################
// Compute_Vector_Gradient
//#####################################################################
template<class TV> bool CELL_CENTERED_OPERATION_UTILITIES<TV>::
Compute_Vector_Gradient(const GRID<TV>& grid,
        const MPI_UNIFORM_GRID<T_GRID> *mpi_grid,
        const T_ARRAYS_VECTOR& cell_centered_vector_field,
        const T_ARRAYS_VECTOR& cell_centered_vector_field_ghost,
        const short gradient_approximation_order,
        T_ARRAYS_TENSOR& grad_cell_centered_vector_field)
{
    bool result = true;

    //TV_ARRAYS velocity_ghost(grid.Domain_Indices(1));
    //TV_ARRAYS::Put(velocity,velocity_ghost);
    //if(mpi_grid) vector_boundary->Fill_Ghost_Cells(grid,velocity,velocity_ghost,(T)0.,(T)0.,1);

    if (gradient_approximation_order == 2) {
        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();
            for (int axis = 1; axis <= TV::dimension; ++axis) {
                for (int cmpt = 1; cmpt <= TV::dimension; ++cmpt) {
                    if (!cell_centered_vector_field.Valid_Index(cell - TV_INT::Axis_Vector(cmpt)) &&
                            (!mpi_grid || (mpi_grid && !mpi_grid->Neighbor(cmpt, 1)))) {
                        grad_cell_centered_vector_field(cell)(axis)(cmpt) =
                            (cell_centered_vector_field(cell + TV_INT::Axis_Vector(cmpt))(axis) -
                             cell_centered_vector_field(cell)(axis)) * grid.One_Over_DX()(cmpt);
                    } else if (!cell_centered_vector_field.Valid_Index(cell + TV_INT::Axis_Vector(cmpt)) &&
                            (!mpi_grid || (mpi_grid && !mpi_grid->Neighbor(cmpt, 2)))) {
                        grad_cell_centered_vector_field(cell)(axis)(cmpt) =
                            (cell_centered_vector_field(cell)(axis) -
                             cell_centered_vector_field(cell - TV_INT::Axis_Vector(cmpt))(axis)) * grid.One_Over_DX()(cmpt);
                    } else {
                        grad_cell_centered_vector_field(cell)(axis)(cmpt) = (T).5 * (
                                cell_centered_vector_field_ghost(cell + TV_INT::Axis_Vector(cmpt))(axis) -
                                cell_centered_vector_field_ghost(cell - TV_INT::Axis_Vector(cmpt))(axis)) *
                            grid.One_Over_DX()(cmpt);
                    }
                }
            }
        }
    }

    return result;
}
//#####################################################################
template class CELL_CENTERED_OPERATION_UTILITIES<VECTOR<float,1> >;
template class CELL_CENTERED_OPERATION_UTILITIES<VECTOR<float,2> >;
template class CELL_CENTERED_OPERATION_UTILITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CELL_CENTERED_OPERATION_UTILITIES<VECTOR<double,1> >;
template class CELL_CENTERED_OPERATION_UTILITIES<VECTOR<double,2> >;
template class CELL_CENTERED_OPERATION_UTILITIES<VECTOR<double,3> >;
#endif
