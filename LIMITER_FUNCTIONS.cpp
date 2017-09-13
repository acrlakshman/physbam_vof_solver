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
#include "LIMITER_FUNCTIONS.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LIMITER_FUNCTIONS<TV>::
LIMITER_FUNCTIONS()
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0)*/
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LIMITER_FUNCTIONS<TV>::
~LIMITER_FUNCTIONS()
{
}
//#####################################################################
template class LIMITER_FUNCTIONS<VECTOR<float,1> >;
template class LIMITER_FUNCTIONS<VECTOR<float,2> >;
template class LIMITER_FUNCTIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LIMITER_FUNCTIONS<VECTOR<double,1> >;
template class LIMITER_FUNCTIONS<VECTOR<double,2> >;
template class LIMITER_FUNCTIONS<VECTOR<double,3> >;
#endif
