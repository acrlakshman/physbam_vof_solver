//#####################################################################
// Copyright 2017, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_INTERPOLATION
//#####################################################################
#ifndef __SURFACE_INTERPOLATION__
#define __SURFACE_INTERPOLATION__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include "INTERPOLATION_WEIGHTS.h"
#include "Tools/utilities.h"
#include "Tools/constants_enums.h"

namespace PhysBAM {

template<class TV>
class SURFACE_INTERPOLATION
{
    typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef VECTOR<T,TV::mixed_partial_derivatives> T_MIXED_DERIVATIVES;
    typedef ARRAY<T_MIXED_DERIVATIVES,TV_INT> T_MIXED_DERIVATIVES_ARRAYS;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef ARRAY<int,TV_INT> INT_ARRAYS;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef ARRAY<TV,TV_INT> T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    //typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_VECTOR_ARRAYS T_FACE_ARRAYS_VECTOR; // TODO
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename REBIND_LENGTH<T_BOUNDARY_SCALAR,TV::dimension>::TYPE T_BOUNDARY_VECTOR;
    typedef typename REBIND_LENGTH<T_BOUNDARY_VECTOR,TV::dimension>::TYPE T_BOUNDARY_TENSOR;
    typedef typename REBIND_LENGTH<T_BOUNDARY_SCALAR,TV::mixed_partial_derivatives>::TYPE T_BOUNDARY_MIXED_DERIVATIVES;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    enum {d=TV::dimension};

  public:

    /*
    T_GRID grid;
    MPI_UNIFORM_GRID<T_GRID> *mpi_grid;

    T_BOUNDARY_SCALAR boundary_scalar;
    T_BOUNDARY_SCALAR *boundary;
    T_BOUNDARY_VECTOR vector_boundary_scalar;
    T_BOUNDARY_VECTOR *vector_boundary;
    */

    SURFACE_INTERPOLATION();
    ~SURFACE_INTERPOLATION();

//#####################################################################
    // Interpolate_Cell_To_Face_Scalar
    bool Interpolate_Cell_To_Face_Scalar(const GRID<TV>& grid,
            const T_FACE_ARRAYS_SCALAR& face_velocity_flux,
            const T_FACE_ARRAYS_SCALAR& surface_area_field,
            const T_ARRAYS_SCALAR& cell_centered_field_ghost,
            const T_ARRAYS_VECTOR& grad_cell_centered_field_ghost,
            const SurfaceInterpolationScheme surface_interpolation_scheme,
            const T surface_interpolation_scheme_limiter_factor, // similar to limitedLinear 1 (<= 1 is the factor here)
            T_FACE_ARRAYS_SCALAR& face_centered_field);
//#####################################################################
};
}
#endif // __SURFACE_INTERPOLATION__
