//#####################################################################
// Copyright 2016, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_REINITIALIZATION
//#####################################################################

#include "LEVELSET_REINITIALIZATION.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_REINITIALIZATION<TV>::
LEVELSET_REINITIALIZATION(
        GRID<TV> &grid_input,
       MPI_UNIFORM_GRID<T_GRID> *mpi_grid_input,
       T_BOUNDARY_SCALAR *boundary_input,
       CUSTOM_PARSE_ARGS<TV> &custom_parse_args_input
    ) : grid(grid_input),
        boundary(boundary_input),
        mpi_grid(mpi_grid_input),
        custom_parse_args(custom_parse_args_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_REINITIALIZATION<TV>::
~LEVELSET_REINITIALIZATION()
{
}
//#####################################################################
// Solve_Levelset_Advection
//#####################################################################
template<class TV> bool LEVELSET_REINITIALIZATION<TV>::
Solve_Levelset_Reinitialization(
        T_ARRAYS &phi,
        const T dtau,
        const int temporal_scheme,
        const int hj_discrete_approximation,
        const int convection_discretization_method,
        const int max_iterations,
        const int number_of_ghost_cells
    )
{
    bool return_value = false;

    switch(temporal_scheme) {
    case 0:
        // No advection
        return_value = true;
        break;

    case 1:
        // Forward Euler
        return_value = Solve_Levelset_Reinitialization_Forward_Euler(phi, dtau,
                                                              hj_discrete_approximation,
                                                              convection_discretization_method,
                                                              max_iterations,
                                                              number_of_ghost_cells);
        break;

    case 2:
        // Second order TVD RK
        return_value = Solve_Levelset_Reinitialization_Forward_TVD_RK_Second_Order(phi, dtau,
                                                              hj_discrete_approximation,
                                                              convection_discretization_method,
                                                              max_iterations,
                                                              number_of_ghost_cells);
        break;

    default:
        PHYSBAM_NOT_IMPLEMENTED();
    }

    return return_value;
}
//#####################################################################
// Solve_Levelset_Reinitialization_Forward_Euler
//#####################################################################
template<class TV> bool LEVELSET_REINITIALIZATION<TV>::
Solve_Levelset_Reinitialization_Forward_Euler(
        T_ARRAYS &phi,
        const T dtau,
        const int hj_discrete_approximation,
        const int convection_discretization_method,
        const int max_iterations,
        const int number_of_ghost_cells
    )
{
    bool return_value = false;

    T_ARRAYS phi_initial = phi;
    // phi_ghost
    T_ARRAYS phi_ghost(grid.Domain_Indices(number_of_ghost_cells));

    // convection_term
    T_ARRAYS convection_term(grid.Domain_Indices());

    for (int iteration = 1; iteration <= max_iterations; ++iteration) {

        // Populate phi_ghost
        T_ARRAYS::Put(phi, phi_ghost);
        boundary->Fill_Ghost_Cells(grid, phi, phi_ghost, (T)0, (T)0, number_of_ghost_cells);

        return_value = Compute_Convection_Term(convection_term, phi_initial, phi_ghost, hj_discrete_approximation, convection_discretization_method);
        if (return_value) {
            for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
                const T_INDEX &cell = iterator.Cell_Index();
                phi(cell) -= (dtau * convection_term(cell));
            }
        }

    }

    return return_value;
}
//#####################################################################
// Solve_Levelset_Reinitialization_Forward_TVD_RK_Second_Order
//#####################################################################
template<class TV> bool LEVELSET_REINITIALIZATION<TV>::
Solve_Levelset_Reinitialization_Forward_TVD_RK_Second_Order(
        T_ARRAYS &phi,
        const T dtau,
        const int hj_discrete_approximation,
        const int convection_discretization_method,
        const int max_iterations,
        const int number_of_ghost_cells
    )
{
    bool return_value = false;

    T_ARRAYS phi_initial = phi;
    // phi_ghost
    T_ARRAYS phi_ghost(grid.Domain_Indices(number_of_ghost_cells));

    // convection_term
    T_ARRAYS convection_term(grid.Domain_Indices());

    for (int iteration = 1; iteration <= max_iterations; ++iteration) {

        T_ARRAYS phi_current = phi;

        for (int rk_step = 1; rk_step <= 2; ++rk_step) {
            // Populate phi_ghost
            T_ARRAYS::Put(phi, phi_ghost);
            boundary->Fill_Ghost_Cells(grid, phi, phi_ghost, (T)0, (T)0, number_of_ghost_cells);

            return_value = Compute_Convection_Term(convection_term, phi_initial, phi_ghost, hj_discrete_approximation, convection_discretization_method);
            if (return_value) {
                for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
                    const T_INDEX &cell = iterator.Cell_Index();
                    phi(cell) -= (dtau * convection_term(cell));
                }
            }
        }

        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const T_INDEX &cell = iterator.Cell_Index();
            phi(cell) += phi_current(cell);
            phi(cell) *= (T)0.5;
        }

    }

    return return_value;
}

//#####################################################################
// Compute_Convection_Term
//#####################################################################
template<class TV> bool LEVELSET_REINITIALIZATION<TV>::
Compute_Convection_Term(
        T_ARRAYS &convection_term,
        const T_ARRAYS &phi_initial,
        const T_ARRAYS &phi_ghost,
        const int hj_discrete_approximation,
        const int discretization_method
    )
{
    bool return_value = false;

    switch(hj_discrete_approximation) {
    case 0:
        // Use this to skip convection discretization
        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const T_INDEX &cell = iterator.Cell_Index();
            convection_term(cell) = (T)0;
        }
        return_value = true;
        break;

    case 1:
        // Godunov's method
        return_value = Compute_Convection_Term_Godunov_WENO5(convection_term, phi_initial, phi_ghost);
        return_value = true;
        break;

    default:
        PHYSBAM_NOT_IMPLEMENTED();
    }

    return return_value;
}
//#####################################################################
// Compute_Convection_Term_Godunov_WENO5
//#####################################################################
template<class TV> bool LEVELSET_REINITIALIZATION<TV>::
Compute_Convection_Term_Godunov_WENO5(
        T_ARRAYS &convection_term,
        const T_ARRAYS &phi_initial,
        const T_ARRAYS &phi_ghost
    )
{
    bool return_value = false;

    for (int axis = 1; axis <= TV::dimension; ++axis) {
        PHYSBAM_ASSERT(phi_ghost.Size()(axis) == grid.Counts()(axis)+6);
    }

    T dx_min = grid.DX().Min();
    TV one_over_dx = grid.One_Over_DX();

    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const T_INDEX &cell = iterator.Cell_Index();
        convection_term(cell) = (T)0;

        TV variable_grad_plus, variable_grad_minus;

        for (int cmpt = 1; cmpt <= TV::dimension; ++cmpt) {
            T variable_WENO5_stencil[7];
            for (int i = 0, offset = -3; i < 7; ++i, ++offset)
                variable_WENO5_stencil[i] = phi_ghost(cell+offset*T_INDEX::Axis_Vector(cmpt));

            T small_num_=(T)1e-6;
            // dvariable_by_dx_minus
            // WENO "v's"
            T v1 = one_over_dx(cmpt)*(variable_WENO5_stencil[1] - variable_WENO5_stencil[0]);
            T v2 = one_over_dx(cmpt)*(variable_WENO5_stencil[2] - variable_WENO5_stencil[1]);
            T v3 = one_over_dx(cmpt)*(variable_WENO5_stencil[3] - variable_WENO5_stencil[2]);
            T v4 = one_over_dx(cmpt)*(variable_WENO5_stencil[4] - variable_WENO5_stencil[3]);
            T v5 = one_over_dx(cmpt)*(variable_WENO5_stencil[5] - variable_WENO5_stencil[4]);

            // Smoothness
            T S1 = thirteen_over_twelve*sqr(v1-(T)2.*v2+v3) + (T).25*sqr(v1-(T)4.*v2+(T)3.*v3);
            T S2 = thirteen_over_twelve*sqr(v2-(T)2.*v3+v4) + (T).25*sqr(v2-v4);
            T S3 = thirteen_over_twelve*sqr(v3-(T)2.*v4+v5) + (T).25*sqr((T)3.*v3-(T)4.*v4+v5);

            // Weights
            T a1 = (T).1/sqr(small_num_+S1);
            T a2 = (T).6/sqr(small_num_+S2);
            T a3 = (T).3/sqr(small_num_+S3);
            T one_over_a123 = (T)1./(a1+a2+a3);
            T w1 = a1*one_over_a123;
            T w2 = a2*one_over_a123;
            T w3 = a3*one_over_a123;

            T dvariable_by_dx_minus = (w1*(one_third*v1-(T)7.*one_sixth*v2+(T)11.*one_sixth*v3)) +
                (w2*(-one_sixth*v2+(T)5.*one_sixth*v3+one_third*v4)) +
                (w3*(one_third*v3+(T)5.*one_sixth*v4-one_sixth*v5));

            // dvariable_by_dx_plus
            // WENO "v's"
            v1 = one_over_dx(cmpt)*(variable_WENO5_stencil[6] - variable_WENO5_stencil[5]);
            v2 = one_over_dx(cmpt)*(variable_WENO5_stencil[5] - variable_WENO5_stencil[4]);
            v3 = one_over_dx(cmpt)*(variable_WENO5_stencil[4] - variable_WENO5_stencil[3]);
            v4 = one_over_dx(cmpt)*(variable_WENO5_stencil[3] - variable_WENO5_stencil[2]);
            v5 = one_over_dx(cmpt)*(variable_WENO5_stencil[2] - variable_WENO5_stencil[1]);

            // Smoothness
            S1 = thirteen_over_twelve*sqr(v1-(T)2.*v2+v3)+(T).25*sqr(v1-(T)4.*v2+(T)3.*v3);
            S2 = thirteen_over_twelve*sqr(v2-(T)2.*v3+v4)+(T).25*sqr(v2-v4);
            S3 = thirteen_over_twelve*sqr(v3-(T)2.*v4+v5)+(T).25*sqr((T)3.*v3-(T)4.*v4+v5);

            // Weights
            a1 = (T).1/sqr(small_num_+S1);
            a2 = (T).6/sqr(small_num_+S2);
            a3 = (T).3/sqr(small_num_+S3);
            one_over_a123 = (T)1./(a1+a2+a3);
            w1 = a1*one_over_a123;
            w2 = a2*one_over_a123;
            w3 = a3*one_over_a123;

            T dvariable_by_dx_plus = (w1*(one_third*v1-(T)7.*one_sixth*v2+(T)11.*one_sixth*v3)) +
                (w2*(-one_sixth*v2+(T)5.*one_sixth*v3+one_third*v4)) +
                (w3*(one_third*v3+(T)5.*one_sixth*v4-one_sixth*v5));

            variable_grad_plus[cmpt] = dvariable_by_dx_plus;
            variable_grad_minus[cmpt] = dvariable_by_dx_minus;
        }

        // Compute convection term
        if (phi_initial(cell) > (T)0) {
            for (int cmpt = 1; cmpt <= TV::dimension; ++cmpt) {
                T max_term = max( sqr(max(variable_grad_minus[cmpt], (T)0)),
                                  sqr(min(variable_grad_plus[cmpt], (T)0)) );
                convection_term(cell) += max_term;
            }
            convection_term(cell) = sqrt(convection_term(cell))-(T)1;
        } else if (phi_initial(cell) < (T)0) {
            for (int cmpt = 1; cmpt <= TV::dimension; ++cmpt) {
                T max_term = max( sqr(min(variable_grad_minus[cmpt], (T)0)),
                                  sqr(max(variable_grad_plus[cmpt], (T)0)) );
                convection_term(cell) += max_term;
            }
            convection_term(cell) = sqrt(convection_term(cell))-(T)1;
        }
        convection_term(cell) *= Smooth_Sign_Function_One(phi_initial(cell), dx_min);

    }
    return_value = true;

    return return_value;
}
//#####################################################################
template class LEVELSET_REINITIALIZATION<VECTOR<float,1> >;
template class LEVELSET_REINITIALIZATION<VECTOR<float,2> >;
template class LEVELSET_REINITIALIZATION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_REINITIALIZATION<VECTOR<double,1> >;
template class LEVELSET_REINITIALIZATION<VECTOR<double,2> >;
template class LEVELSET_REINITIALIZATION<VECTOR<double,3> >;
#endif
