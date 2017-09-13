//#####################################################################
// Copyright 2017, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LIMITER_FUNCTIONS
//#####################################################################
#ifndef __LIMITER_FUNCTIONS__
#define __LIMITER_FUNCTIONS__

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
#include "Tools/utilities.h"
#include "Tools/constants_enums.h"

namespace PhysBAM {

template<class TV>
class LIMITER_FUNCTIONS
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

    LIMITER_FUNCTIONS();
    ~LIMITER_FUNCTIONS();

    // Compute r(/R)
    static T R(const T face_velocity_flux, const T alpha1_first_cell,
            const T alpha1_second_cell, const TV grad_alpha1_first_cell,
            const TV grad_alpha1_second_cell, const TV cell_distance_second_first)
    {
        T r = (T)0;

        T diff_alpha1_f = alpha1_second_cell - alpha1_first_cell;
        T diff_alpha1_cf = (T)0;

        if (face_velocity_flux > (T)0) {
            diff_alpha1_cf = Dot_Product<TV>(cell_distance_second_first, grad_alpha1_first_cell);
        } else {
            diff_alpha1_cf = Dot_Product<TV>(cell_distance_second_first, grad_alpha1_second_cell);
        }

        if (abs(diff_alpha1_cf) >= (T)1000*abs(diff_alpha1_f)) {
            r = (T)2 * (T)1000 * sign(diff_alpha1_cf) * sign(diff_alpha1_f) - (T)1;
        } else {
            r = (T)2 * (diff_alpha1_cf / diff_alpha1_f) - (T)1;
        }

        return r;
    }

    // Compute limiter
    static T Limiter(const T geometric_weight, const T face_velocity_flux,
            const T alpha1_first_cell, const T alpha1_second_cell,
            const TV grad_alpha1_first_cell, const TV grad_alpha1_second_cell,
            const TV cell_distance_second_first,
            const SurfaceInterpolationScheme surface_interpolation_scheme,
            const T factor)
    {
        T limiter = (T)0;

        T r = R(face_velocity_flux, alpha1_first_cell,
                alpha1_second_cell, grad_alpha1_first_cell,
                grad_alpha1_second_cell, cell_distance_second_first);

        switch (surface_interpolation_scheme) {
            case UPWIND:
                if (face_velocity_flux >= (T)0)
                    limiter = (T)0;
                else
                    // NOTE: assuming uniform grid
                    //limiter = (T)2;
                    limiter = (T)1 / geometric_weight;
                break;

            case LINEAR:
                limiter = (T)1;
                break;

            case LIMITED_LINEAR: {
                    T twoByFactor = (T)2 / max(factor, (T)SMALL_NUMBER);
                    limiter = max(min(twoByFactor*r, (T)1), (T)0);
                }
                break;

            case SUPERBEE:
                limiter = max(max(min((T)2*r, (T)1), min(r, (T)2)), (T)0);
                break;

            case MINMOD:
                limiter = max(min(r, (T)1), (T)0);
                break;

            case VANALBADA:
                limiter = r*(r + (T)1) / (sqr(r) + (T)1);
                break;

            case VANLEER:
                limiter = (r + abs(r)) / ((T)1 + abs(r));
                break;

            case MUSCL:
                limiter = max(min(min((T)2*r, (T)0.5*r + (T)0.5), (T)2), (T)0);
                break;

            case INTERFACE_COMPRESSION:
                // TODO: Check
                limiter = min(max((T)1 - max(sqr((T)1 - (T)4*alpha1_first_cell*((T)1 - alpha1_first_cell)),
                                sqr((T)1 - (T)4*alpha1_second_cell*((T)1 - alpha1_second_cell))), (T)0), (T)1);
                break;

            default:
                PHYSBAM_NOT_IMPLEMENTED();
        }

        return limiter;
    }

//#####################################################################
//#####################################################################
};

}
#endif // __LIMITER_FUNCTIONS__
