//#####################################################################
// Copyright 2016, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_REINITIALIZATION
//#####################################################################
#ifndef __LEVELSET_REINITIALIZATION__
#define __LEVELSET_REINITIALIZATION__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include "Tools/utilities.h"

namespace PhysBAM{

template<class TV>
class LEVELSET_REINITIALIZATION{
  public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<TV,TV::dimension> T_TENSOR;
    typedef VECTOR<int,TV::dimension> T_INDEX;
    typedef VECTOR<T,TV::mixed_partial_derivatives> T_MIXED_DERIVATIVES;
    typedef ARRAY<T,T_INDEX> T_ARRAYS;
    typedef ARRAY<TV,T_INDEX> TV_ARRAYS;
    typedef ARRAY<T_TENSOR,T_INDEX> T_TENSOR_ARRAYS;
    typedef ARRAY<T_MIXED_DERIVATIVES,T_INDEX> T_MIXED_DERIVATIVES_ARRAYS;
    typedef ARRAY<int,T_INDEX> INT_ARRAYS;
    typedef GRID<TV> T_GRID;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;

    T_GRID grid;
    T_BOUNDARY_SCALAR *boundary;
    MPI_UNIFORM_GRID<T_GRID> *mpi_grid;
    CUSTOM_PARSE_ARGS<TV> custom_parse_args;
    VECTOR<VECTOR<int,2>,TV::dimension> domain_boundary_conditions;
    VECTOR<VECTOR<T,2>,TV::dimension> domain_boundary_values;

    LEVELSET_REINITIALIZATION(
            T_GRID &grid_input,
            MPI_UNIFORM_GRID<T_GRID> *mpi_grid_input,
            T_BOUNDARY_SCALAR *boundary_input,
            CUSTOM_PARSE_ARGS<TV> &custom_parse_args_input
        );
    ~LEVELSET_REINITIALIZATION();

    T Smooth_Sign_Function_One(T phi, T dx) {
        return phi / (sqrt((phi*phi) + (dx*dx)));
    }

    T Smooth_Sign_Function_Two(T phi, T mag_grad_phi, T dx) {
        return phi / (sqrt((phi*phi) + (mag_grad_phi*mag_grad_phi*dx*dx)));
    }

//#####################################################################
    bool Solve_Levelset_Reinitialization(
            T_ARRAYS &phi,
            const T dtau,
            const int temporal_scheme,
            const int hj_discrete_approximation,
            const int convection_discretization_method,
            const int max_iterations,
            const int number_of_ghost_cells = 3
        );
    bool Solve_Levelset_Reinitialization_Forward_Euler(
            T_ARRAYS &phi,
            const T dtau,
            const int hj_discrete_approximation,
            const int convection_discretization_method,
            const int max_iterations,
            const int number_of_ghost_cells = 3
        );
    bool Solve_Levelset_Reinitialization_Forward_TVD_RK_Second_Order(
            T_ARRAYS &phi,
            const T dtau,
            const int hj_discrete_approximation,
            const int convection_discretization_method,
            const int max_iterations,
            const int number_of_ghost_cells = 3
        );
    bool Compute_Convection_Term(
            T_ARRAYS &convection_term,
            const T_ARRAYS &phi_initial,
            const T_ARRAYS &phi_ghost,
            const int hj_discrete_approximation,
            const int discretization_method
        );
    bool Compute_Convection_Term_Godunov_WENO5(
            T_ARRAYS &convection_term,
            const T_ARRAYS &phi_initial,
            const T_ARRAYS &phi_ghost
        );
//#####################################################################
};
}
#endif
