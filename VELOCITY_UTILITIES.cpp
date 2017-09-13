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
#include "VELOCITY_UTILITIES.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VELOCITY_UTILITIES<TV>::
VELOCITY_UTILITIES()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> VELOCITY_UTILITIES<TV>::
~VELOCITY_UTILITIES()
{
}
//#####################################################################
// Interpolate_Cell_To_Face_Linear
//#####################################################################
template<class TV> bool VELOCITY_UTILITIES<TV>::
Interpolate_Cell_To_Face_Linear(const GRID<TV>& grid, const T_ARRAYS_VECTOR& velocity_ghost,
        T_FACE_ARRAYS_SCALAR& face_velocities)
{
    bool result = true;

    FACE_CENTERED_OPERATION_UTILITIES<TV> face_centered_operation_utilities;
    result = face_centered_operation_utilities.Vector_Interpolate_Cell_To_Face_Linear(grid, velocity_ghost, face_velocities);

    return result;
}
//#####################################################################
// Populate_Ghost_Values_Using_Boundary_Conditions
//#####################################################################
template<class TV> bool VELOCITY_UTILITIES<TV>::
Populate_Ghost_Values_Using_Boundary_Conditions(const GRID<TV>& grid,
        const int number_of_ghost_cells,
        const MPI_UNIFORM_GRID<T_GRID> *mpi_grid,
        const VECTOR<VECTOR<VECTOR<int,TV::dimension>,2>,TV::dimension> face_velocity_boundary_conditions,
        const VECTOR<VECTOR<VECTOR<T,TV::dimension>,2>,TV::dimension> face_velocity_boundary_values,
        T_ARRAYS_VECTOR& cell_centered_velocity_field_ghost
        )
{
    bool result = true;

    // Populate ghost velocity values using boundary conditions
    for (int cmpt = 1; cmpt <= TV::dimension; ++cmpt) {
        for (int axis = 1; axis <= TV::dimension; ++axis) {
            for (int axis_side = 1; axis_side <= 2; ++axis_side) {
                int side = 2*(axis - 1) + axis_side;
                TV_INT interior_cell_offset = (axis_side == 1) ? TV_INT() : -TV_INT::Axis_Vector(axis);
                TV_INT exterior_cell_offset = (axis_side == 1) ? -TV_INT::Axis_Vector(axis) : TV_INT();
                TV_INT boundary_face_offset = (axis_side == 1) ? TV_INT::Axis_Vector(axis) : -TV_INT::Axis_Vector(axis);

                if (!mpi_grid || (mpi_grid && !mpi_grid->Neighbor(axis, axis_side))) {
                    for (FACE_ITERATOR iterator(grid, 1, GRID<TV>::BOUNDARY_REGION, side); iterator.Valid(); iterator.Next()) {
                        if (face_velocity_boundary_conditions[axis][axis_side][cmpt] == 1) {
                            // Dirichlet boundary condition
                            TV_INT face = iterator.Face_Index() + boundary_face_offset;
                            TV_INT cell = face + exterior_cell_offset;
                            cell_centered_velocity_field_ghost(cell)(cmpt) =
                                ((T)2 * face_velocity_boundary_values[axis][axis_side][cmpt] -
                                 cell_centered_velocity_field_ghost(face + interior_cell_offset)(cmpt));

                            if (number_of_ghost_cells > 1) {
                                cell_centered_velocity_field_ghost(cell + exterior_cell_offset)(cmpt) =
                                    cell_centered_velocity_field_ghost(cell)(cmpt);
                            }
                            if (number_of_ghost_cells > 2) {
                                cell_centered_velocity_field_ghost(cell + 2*exterior_cell_offset)(cmpt) =
                                    cell_centered_velocity_field_ghost(cell + exterior_cell_offset)(cmpt);
                            }
                        } else if (face_velocity_boundary_conditions[axis][axis_side][cmpt] == 2) {
                            // Neumann boundary condition
                            TV_INT face = iterator.Face_Index() + boundary_face_offset;
                            TV_INT cell = face + exterior_cell_offset;
                            cell_centered_velocity_field_ghost(cell)(cmpt) =
                                cell_centered_velocity_field_ghost(face + interior_cell_offset)(cmpt);

                            if (number_of_ghost_cells > 1) {
                                cell_centered_velocity_field_ghost(cell + exterior_cell_offset)(cmpt) =
                                    cell_centered_velocity_field_ghost(cell)(cmpt);
                            }
                            if (number_of_ghost_cells > 2) {
                                cell_centered_velocity_field_ghost(cell + 2*exterior_cell_offset)(cmpt) =
                                    cell_centered_velocity_field_ghost(cell + exterior_cell_offset)(cmpt);
                            }
                        }
                    }
                }
            }
        }
    }

    return result;
}
//#####################################################################
template class VELOCITY_UTILITIES<VECTOR<float,1> >;
template class VELOCITY_UTILITIES<VECTOR<float,2> >;
template class VELOCITY_UTILITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VELOCITY_UTILITIES<VECTOR<double,1> >;
template class VELOCITY_UTILITIES<VECTOR<double,2> >;
template class VELOCITY_UTILITIES<VECTOR<double,3> >;
#endif
