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
#include "SURFACE_AREA.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_AREA<TV>::
SURFACE_AREA()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SURFACE_AREA<TV>::
~SURFACE_AREA()
{
}
//#####################################################################
// Compute
//#####################################################################
template<class TV> bool SURFACE_AREA<TV>::
Compute(const GRID<TV>& grid, T_FACE_ARRAYS_SCALAR& surface_area_field)
{
    bool status = true;

    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        TV_INT first_cell = iterator.First_Cell_Index();
        TV_INT second_cell = iterator.Second_Cell_Index();

        surface_area_field(axis, face_index) = (T)1;

        for (int d = 1; d <= TV::dimension; ++d) {
            if (d != axis) {
                surface_area_field(axis, face_index) *= grid.DX()(d);
            }
        }
    }

    return status;
}
//#####################################################################
// Get_Area_Vector
//#####################################################################
template<class TV> bool SURFACE_AREA<TV>::
Get_Area_Vector(const GRID<TV>& grid, const T_FACE_ARRAYS_SCALAR& surface_area_field, const T_FACE_LOOKUP& face_lookup, const TV_INT& cell,
                const int axis, VECTOR<TV,2>& area_vector)
{
    bool status = false;

    /* // TODO
    const typename T_FACE_LOOKUP::LOOKUP& lookup = face_lookup.Starting_Point_Cell(cell);
    TV area_vector_first_face = TV();
    area_vector_first_face(axis) = lookup.First_Face_Index(axis);
    area_vector(1) = lookup.First_Face_Index(axis);
    area_vector(2) = lookup.Second_Face_Index(axis);
    */

    return status;
}
//#####################################################################
template class SURFACE_AREA<VECTOR<float,1> >;
template class SURFACE_AREA<VECTOR<float,2> >;
template class SURFACE_AREA<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SURFACE_AREA<VECTOR<double,1> >;
template class SURFACE_AREA<VECTOR<double,2> >;
template class SURFACE_AREA<VECTOR<double,3> >;
#endif
