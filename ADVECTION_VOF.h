//#####################################################################
// Copyright 2017, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_VOF
//#####################################################################
#ifndef __ADVECTION_VOF__
#define __ADVECTION_VOF__

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
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include "SURFACE_INTERPOLATION.h"
#include "CELL_CENTERED_OPERATION_UTILITIES.h"
#include "FACE_CENTERED_OPERATION_UTILITIES.h"
#include "MULES.h"
#include "Parallel_Computation/MPI_UNIFORM_COMMUNICATION_HELPER.h"
#include "Tools/utilities.h"
#include "Tools/constants_enums.h"

namespace PhysBAM {

template<class TV>
class ADVECTION_VOF
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

    // Parameters for solution
    TemporalDiscretization temporal_discretization;
    SpatialDiscretization alpha1_convection_term_discretization;
    SpatialDiscretization alpha1_compression_term_discretization;
    SurfaceInterpolationScheme alpha1_convection_term_surface_interpolation;
    SurfaceInterpolationScheme alpha1_compression_term_surface_interpolation;
    SurfaceInterpolationScheme alpha1_velocity_convection_term_surface_interpolation;
    T alpha1_convection_term_surface_interpolation_scheme_limiter_factor;
    T alpha1_compression_term_surface_interpolation_scheme_limiter_factor;
    bool solve_alpha1_convection_term;
    bool solve_alpha1_compression_term;
    bool alpha1_use_mules;

    ADVECTION_VOF();
    ~ADVECTION_VOF();

//#####################################################################
    void Set_Options(CUSTOM_PARSE_ARGS<TV>& custom_parse_args);

    bool Explicit_Solve(const GRID<TV>& grid,
            const MPI_UNIFORM_GRID<T_GRID> *mpi_grid,
            const short num_cells_transfer,
            const T dt,
            const T_ARRAYS_VECTOR& velocity_field_ghost,
            const T_FACE_ARRAYS_SCALAR& face_velocities_field,
            const T_FACE_ARRAYS_SCALAR& face_velocity_flux_field,
            const T_FACE_LOOKUP& face_velocity_flux_field_lookup,
            const T_ARRAYS_SCALAR& alpha1_previous_field_ghost,
            const T_ARRAYS_VECTOR& grad_alpha1_previous_field_ghost,
            const T_FACE_ARRAYS_SCALAR& surface_area_field,
            T_BOUNDARY_SCALAR *boundary,
            T_BOUNDARY_VECTOR *vector_boundary,
            T_FACE_ARRAYS_SCALAR& face_alpha1_flux_field,
            T_ARRAYS_SCALAR& alpha1_field);

    bool Forward_Euler(const GRID<TV>& grid,
            const MPI_UNIFORM_GRID<T_GRID> *mpi_grid,
            const short num_cells_transfer,
            const T dt, const T_ARRAYS_VECTOR& velocity_field_ghost,
            const T_FACE_ARRAYS_SCALAR& face_velocities_field,
            const T_FACE_ARRAYS_SCALAR& face_velocity_flux_field,
            const T_FACE_LOOKUP& face_velocity_flux_field_lookup,
            const T_ARRAYS_SCALAR& alpha1_previous_field_ghost,
            const T_ARRAYS_VECTOR& grad_alpha1_previous_field_ghost,
            const T_FACE_ARRAYS_SCALAR& surface_area_field,
            T_BOUNDARY_SCALAR *boundary,
            T_BOUNDARY_VECTOR *vector_boundary,
            T_FACE_ARRAYS_SCALAR& face_alpha1_flux_field,
            T_ARRAYS_SCALAR& alpha1_field);

    bool Get_MULES_Limiter(const GRID<TV>& grid,
            const MPI_UNIFORM_GRID<T_GRID> *mpi_grid,
            const T rDeltaT,
            const T_ARRAYS_SCALAR& alpha1_field_ghost,
            const T_ARRAYS_SCALAR& alpha1_previous_field_ghost,
            const T_FACE_ARRAYS_SCALAR& face_velocity_flux_field,
            const T_FACE_ARRAYS_SCALAR& face_alpha1_upwind_flux_field,
            //const T_FACE_ARRAYS_SCALAR& face_alpha1_flux_field,
            const T_FACE_ARRAYS_SCALAR& face_alpha1_flux_corr_field,
            const T_ARRAYS_SCALAR& Sp,
            const T_ARRAYS_SCALAR& Su,
            const T alpha_max,
            const T alpha_min,
            T_BOUNDARY_SCALAR *boundary,
            T_FACE_ARRAYS_SCALAR& lambda
            );
//#####################################################################
};
}
#endif // __ADVECTION_VOF__
