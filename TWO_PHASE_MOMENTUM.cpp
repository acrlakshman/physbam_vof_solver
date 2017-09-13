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
#include "TWO_PHASE_MOMENTUM.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TWO_PHASE_MOMENTUM<TV>::
TWO_PHASE_MOMENTUM()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TWO_PHASE_MOMENTUM<TV>::
~TWO_PHASE_MOMENTUM()
{
}
//#####################################################################
// Compute_Phi
//#####################################################################
template<class TV> bool TWO_PHASE_MOMENTUM<TV>::
Compute_Au(const GRID<TV>& grid, const T_FACE_ARRAYS_SCALAR& face_velocities,
        const T_FACE_ARRAYS_SCALAR& surface_area_field,
        T_FACE_ARRAYS_SCALAR& face_velocity_flux)
{
    bool result = true;

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        TV surface_area_vector = TV();
        surface_area_vector(axis) = surface_area_field(axis, face_index);

        TV face_velocity_vector = TV();
        face_velocity_vector(axis) = face_velocities(axis, face_index);

        face_velocity_flux(axis, face_index) = Dot_Product<TV>(face_velocity_vector, surface_area_vector);
    }

    return result;
}
//#####################################################################
template class TWO_PHASE_MOMENTUM<VECTOR<float,1> >;
template class TWO_PHASE_MOMENTUM<VECTOR<float,2> >;
template class TWO_PHASE_MOMENTUM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TWO_PHASE_MOMENTUM<VECTOR<double,1> >;
template class TWO_PHASE_MOMENTUM<VECTOR<double,2> >;
template class TWO_PHASE_MOMENTUM<VECTOR<double,3> >;
#endif
