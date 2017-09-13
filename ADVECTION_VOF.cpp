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
#include "ADVECTION_VOF.h"

using namespace PhysBAM;

//#####################################################################
// Constructor
//#####################################################################
template<class TV> ADVECTION_VOF<TV>::
ADVECTION_VOF() :
    /*grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    mpi_grid(0),
    boundary(0),
    vector_boundary(0),*/
    temporal_discretization(FORWARD_EULER),
    alpha1_convection_term_discretization(GAUSS),
    alpha1_compression_term_discretization(GAUSS),
    alpha1_convection_term_surface_interpolation(SUPERBEE),
    alpha1_compression_term_surface_interpolation(INTERFACE_COMPRESSION),
    alpha1_velocity_convection_term_surface_interpolation(LINEAR),
    alpha1_convection_term_surface_interpolation_scheme_limiter_factor((T)1),
    alpha1_compression_term_surface_interpolation_scheme_limiter_factor((T)1),
    solve_alpha1_convection_term(true),
    solve_alpha1_compression_term(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ADVECTION_VOF<TV>::
~ADVECTION_VOF()
{
}
//#####################################################################
// Set_Options
//#####################################################################
template<class TV> void ADVECTION_VOF<TV>::
Set_Options(CUSTOM_PARSE_ARGS<TV>& custom_parse_args)
{
    this->temporal_discretization = static_cast<TemporalDiscretization>(custom_parse_args.alpha1_temporal_discretization);
    this->solve_alpha1_convection_term = (bool)custom_parse_args.solve_alpha1_convection_term;
    this->solve_alpha1_compression_term = (bool)custom_parse_args.solve_alpha1_compression_term;
    this->alpha1_use_mules = (bool)custom_parse_args.alpha1_use_mules;
    this->alpha1_convection_term_surface_interpolation_scheme_limiter_factor =
        custom_parse_args.alpha1_convection_term_surface_interpolation_scheme_limiter_factor;
    this->alpha1_compression_term_surface_interpolation_scheme_limiter_factor =
        custom_parse_args.alpha1_compression_term_surface_interpolation_scheme_limiter_factor;
    this->alpha1_convection_term_discretization =
        static_cast<SpatialDiscretization>(custom_parse_args.alpha1_convection_term_discretization);
    this->alpha1_compression_term_discretization =
        static_cast<SpatialDiscretization>(custom_parse_args.alpha1_compression_term_discretization);
    this->alpha1_convection_term_surface_interpolation =
        static_cast<SurfaceInterpolationScheme>(custom_parse_args.alpha1_convection_term_surface_interpolation);
    this->alpha1_compression_term_surface_interpolation =
        static_cast<SurfaceInterpolationScheme>(custom_parse_args.alpha1_compression_term_surface_interpolation);
    this->alpha1_velocity_convection_term_surface_interpolation =
        static_cast<SurfaceInterpolationScheme>(custom_parse_args.alpha1_velocity_convection_term_surface_interpolation);
}
//#####################################################################
// Explicit_Solve
//#####################################################################
template<class TV> bool ADVECTION_VOF<TV>::
Explicit_Solve(const GRID<TV>& grid,
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
        T_ARRAYS_SCALAR& alpha1_field)
{
    bool result = false;

    switch (temporal_discretization) {
    case FORWARD_EULER:
        result = Forward_Euler(grid,
                mpi_grid,
                num_cells_transfer,
                dt,
                velocity_field_ghost,
                face_velocities_field,
                face_velocity_flux_field,
                face_velocity_flux_field_lookup,
                alpha1_previous_field_ghost,
                grad_alpha1_previous_field_ghost,
                surface_area_field,
                boundary,
                vector_boundary,
                face_alpha1_flux_field,
                alpha1_field);
        break;

    default:
        break;
    }

    return result;
}
//#####################################################################
// Forward_Euler
//#####################################################################
template<class TV> bool ADVECTION_VOF<TV>::
Forward_Euler(const GRID<TV>& grid,
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
        T_ARRAYS_SCALAR& alpha1_field)
{
    bool result = true;

    T convection_term = (T)0;
    T compression_term = (T)0;
    T alpha1_flux = (T)0;
    T alpha1_face = (T)0;
    //T phi_face = (T)0;
    T phi_r_face = (T)0;
    T alpha1_r_face = (T)0;
    T one_minus_alpha1_r_face = (T)0;
    T c_alpha1 = (T)1;

    SURFACE_INTERPOLATION<TV> surface_interpolation;
    CELL_CENTERED_OPERATION_UTILITIES<TV> cell_centered_operation_utilities;
    FACE_CENTERED_OPERATION_UTILITIES<TV> face_centered_operation_utilities;

    T_ARRAYS_SCALAR zero_field(grid.Domain_Indices()); // TODO: Make this static
    T_ARRAYS_VECTOR alpha1_normalized_field_ghost;
    T_FACE_ARRAYS_SCALAR face_alpha1_r_flux_field;
    T_FACE_ARRAYS_SCALAR face_alpha1_upwind_flux_field;
    T_FACE_ARRAYS_SCALAR face_alpha1_normalized_field;
    T_FACE_ARRAYS_SCALAR face_alpha1_normalized_dot_surface_area_field;

    if (solve_alpha1_convection_term) {
        //face_alpha1_flux_field.Resize(grid.Domain_Indices());
        // Interpolate alpha1_previous_field to face_alpha1_flux_field (to save memory)
        if (not (surface_interpolation.Interpolate_Cell_To_Face_Scalar(grid,
                    face_velocity_flux_field,
                    surface_area_field,
                    alpha1_previous_field_ghost,
                    grad_alpha1_previous_field_ghost,
                    alpha1_convection_term_surface_interpolation,
                    alpha1_convection_term_surface_interpolation_scheme_limiter_factor,
                    face_alpha1_flux_field))) {
            PHYSBAM_ERROR("ADVECTION_VOF.cpp");
        }
    }

    if (solve_alpha1_compression_term) {
        alpha1_normalized_field_ghost.Resize(grid.Domain_Indices(GLOBAL_NUM_CELLS_TRANSFER));
        face_alpha1_r_flux_field.Resize(grid, 0, false);
        face_alpha1_normalized_field.Resize(grid, 0, false);
        face_alpha1_normalized_dot_surface_area_field.Resize(grid, 0, false);

        cell_centered_operation_utilities.Normalize_Cell_Centered_Vector_Field(grid,
                grad_alpha1_previous_field_ghost,
                alpha1_normalized_field_ghost);
        // Transfer alpha1_normalized_field across processors
        T_ARRAYS_VECTOR::Put(alpha1_normalized_field_ghost,
                alpha1_normalized_field_ghost);
        vector_boundary->Fill_Ghost_Cells(grid,
                alpha1_normalized_field_ghost,
                alpha1_normalized_field_ghost,
                (T)0,
                (T)0,
                GLOBAL_NUM_CELLS_TRANSFER);
        if (mpi_grid) Exchange_Boundary_Vectors(*mpi_grid,
                alpha1_normalized_field_ghost,
                GLOBAL_NUM_CELLS_TRANSFER);

        // Interpolate alpha1_previous_field to face_alpha1_r_flux_field (to save memory)
        if (not (surface_interpolation.Interpolate_Cell_To_Face_Scalar(grid,
                    face_velocity_flux_field,
                    surface_area_field,
                    alpha1_previous_field_ghost,
                    grad_alpha1_previous_field_ghost,
                    alpha1_compression_term_surface_interpolation,
                    alpha1_compression_term_surface_interpolation_scheme_limiter_factor,
                    face_alpha1_r_flux_field))) {
            PHYSBAM_ERROR("ADVECTION_VOF.cpp");
        }

        face_centered_operation_utilities.Vector_Interpolate_Cell_To_Face_Linear(grid,
                alpha1_normalized_field_ghost,
                face_alpha1_normalized_field);

        face_centered_operation_utilities.Dot_Product_Face_Fields(grid,
                face_alpha1_normalized_field,
                surface_area_field,
                face_alpha1_normalized_dot_surface_area_field);
    }

    // Compute face_alpha1_flux_field
    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        //TV_INT first_cell = iterator.First_Cell_Index();
        //TV_INT second_cell = iterator.Second_Cell_Index();

        face_alpha1_flux_field(axis, face_index) *= face_velocity_flux_field(axis, face_index);
    }

    // Compute face_alpha1_r_flux_field
    if (solve_alpha1_compression_term) {
        T phi_r_face = (T)0;
        T one_by_surface_area_mag = (T)0;

        for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            FACE_INDEX<TV::dimension> face = iterator.Full_Index();
            int axis = face.axis;
            TV_INT face_index = face.index;
            //TV_INT first_cell = iterator.First_Cell_Index();
            //TV_INT second_cell = iterator.Second_Cell_Index();

            phi_r_face = c_alpha1 * abs(face_velocity_flux_field(axis, face_index));
            one_by_surface_area_mag = (T)1. / abs(surface_area_field(axis, face_index)); // CHECK IF WORKS
            phi_r_face *= one_by_surface_area_mag;
            phi_r_face *= face_alpha1_normalized_dot_surface_area_field(axis, face_index); // CHECK IF WORKS
            face_alpha1_r_flux_field(axis, face_index) = (phi_r_face * face_alpha1_r_flux_field(axis, face_index) *
                    ((T)1 - face_alpha1_r_flux_field(axis, face_index)));
        }
    }

    if (alpha1_use_mules) {
        // face_alpha1_upwind_flux_field
        face_alpha1_upwind_flux_field.Resize(grid.Domain_Indices());
        // Interpolate alpha1_previous_field to face_alpha1_upwind_flux_field (to save memory)
        if (not (surface_interpolation.Interpolate_Cell_To_Face_Scalar(grid,
                        face_velocity_flux_field,
                        surface_area_field,
                        alpha1_previous_field_ghost,
                        grad_alpha1_previous_field_ghost,
                        UPWIND,
                        alpha1_convection_term_surface_interpolation_scheme_limiter_factor,
                        face_alpha1_upwind_flux_field))) {
            PHYSBAM_ERROR("ADVECTION_VOF.cpp");
        }

        // Compute face_alpha1_upwind_flux_field
        for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            FACE_INDEX<TV::dimension> face = iterator.Full_Index();
            int axis = face.axis;
            TV_INT face_index = face.index;
            //TV_INT first_cell = iterator.First_Cell_Index();
            //TV_INT second_cell = iterator.Second_Cell_Index();

            face_alpha1_upwind_flux_field(axis, face_index) *= face_velocity_flux_field(axis, face_index);
        }

        // Compute face_alpha1_flux_corr_field
        if (solve_alpha1_compression_term) {
            for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
                FACE_INDEX<TV::dimension> face = iterator.Full_Index();
                int axis = face.axis;
                TV_INT face_index = face.index;
                //TV_INT first_cell = iterator.First_Cell_Index();
                //TV_INT second_cell = iterator.Second_Cell_Index();

                face_alpha1_flux_field(axis, face_index) += (face_alpha1_r_flux_field(axis, face_index) -
                        face_alpha1_upwind_flux_field(axis, face_index));
            }
        } else {
            for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
                FACE_INDEX<TV::dimension> face = iterator.Full_Index();
                int axis = face.axis;
                TV_INT face_index = face.index;
                //TV_INT first_cell = iterator.First_Cell_Index();
                //TV_INT second_cell = iterator.Second_Cell_Index();

                face_alpha1_flux_field(axis, face_index) -= face_alpha1_upwind_flux_field(axis, face_index);
            }
        }
    }

    if (alpha1_use_mules && solve_alpha1_convection_term) {

        T_FACE_ARRAYS_SCALAR lambda(grid, 0, false);
        T rDeltaT = (T)1. / dt;
        Get_MULES_Limiter(grid,
                mpi_grid,
                rDeltaT,
                alpha1_previous_field_ghost, // alpha1_field_ghost // TODO: CHECK
                alpha1_previous_field_ghost,
                face_velocity_flux_field,
                face_alpha1_upwind_flux_field,
                //const T_FACE_ARRAYS_SCALAR& face_alpha1_flux_field,
                face_alpha1_flux_field, // face_alpha1_flux_corr_field,
                zero_field, // Sp,
                zero_field, // Su,
                1, // alpha_max,
                0, // alpha_min,
                boundary,
                lambda
                );

        for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            FACE_INDEX<TV::dimension> face = iterator.Full_Index();
            int axis = face.axis;
            TV_INT face_index = face.index;
            //TV_INT first_cell = iterator.First_Cell_Index();
            //TV_INT second_cell = iterator.Second_Cell_Index();

            face_alpha1_flux_field(axis, face_index) = face_alpha1_upwind_flux_field(axis, face_index) +
                (lambda(axis, face_index) * face_alpha1_flux_field(axis, face_index));
        }

        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();
            alpha1_flux = (T)0;

            for (int axis = 1; axis <= TV::dimension; ++axis) {
                alpha1_flux += face_alpha1_flux_field(axis, iterator.Second_Face_Index(axis));
                alpha1_flux -= face_alpha1_flux_field(axis, iterator.First_Face_Index(axis));
            }
            alpha1_flux /= grid.DX().Product();

            alpha1_field(cell) = alpha1_previous_field_ghost(cell) - dt * alpha1_flux;
        }

    } else {

        for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
            const TV_INT& cell = iterator.Cell_Index();

            //const typename T_FACE_LOOKUP::LOOKUP& lookup_face_velocity_flux_field = face_velocity_flux_field_lookup.Starting_Point_Cell(cell);

            convection_term = (T)0;
            if (solve_alpha1_convection_term) {

                //T phi_face = (T)0;

                for (int axis = 1; axis <= TV::dimension; ++axis) {

                    //phi_face = lookup_face_velocity_flux_field(axis, iterator.Second_Face_Index(axis));
                    //alpha1_face = face_alpha1_field(axis, iterator.Second_Face_Index(axis));
                    //convection_term += (phi_face * alpha1_face);
                    convection_term += face_alpha1_flux_field(axis, iterator.Second_Face_Index(axis));

                    //phi_face = lookup_face_velocity_flux_field(axis, iterator.First_Face_Index(axis));
                    //alpha1_face = face_alpha1_field(axis, iterator.First_Face_Index(axis));
                    //convection_term -= (phi_face * alpha1_face);
                    convection_term -= face_alpha1_flux_field(axis, iterator.First_Face_Index(axis));

                }
                convection_term /= grid.DX().Product();
            }

            compression_term = (T)0;
            if (solve_alpha1_compression_term) {

                //T phi_r_face = (T)0;
                //T one_by_surface_area_mag = (T)0;

                for (int axis = 1; axis <= TV::dimension; ++axis) {

                    //phi_r_face = c_alpha1 * abs(lookup_face_velocity_flux_field(axis, iterator.Second_Face_Index(axis)));
                    //one_by_surface_area_mag = (T)1. / abs(surface_area_field(axis, iterator.Second_Face_Index(axis))); // CHECK IF WORKS
                    //phi_r_face *= one_by_surface_area_mag;
                    //phi_r_face *= face_alpha1_normalized_dot_surface_area_field(axis, iterator.Second_Face_Index(axis)); // CHECK IF WORKS

                    //alpha1_r_face = face_alpha1_field(axis, iterator.Second_Face_Index(axis));
                    //one_minus_alpha1_r_face = face_alpha1_field(axis, iterator.Second_Face_Index(axis));

                    //compression_term += (phi_r_face * alpha1_r_face * one_minus_alpha1_r_face);
                    compression_term += face_alpha1_r_flux_field(axis, iterator.Second_Face_Index(axis));

                    //phi_r_face = c_alpha1 * abs(lookup_face_velocity_flux_field(axis, iterator.First_Face_Index(axis)));
                    //one_by_surface_area_mag = (T)1. / abs(surface_area_field(axis, iterator.First_Face_Index(axis))); // CHECK IF WORKS
                    //phi_r_face *= one_by_surface_area_mag;
                    //phi_r_face *= face_alpha1_normalized_dot_surface_area_field(axis, iterator.First_Face_Index(axis)); // CHECK IF WORKS

                    //alpha1_r_face = face_alpha1_field(axis, iterator.First_Face_Index(axis));
                    //one_minus_alpha1_r_face = face_alpha1_field(axis, iterator.First_Face_Index(axis));

                    //compression_term -= (phi_r_face * alpha1_r_face * one_minus_alpha1_r_face);
                    compression_term -= face_alpha1_r_flux_field(axis, iterator.First_Face_Index(axis));

                }
                compression_term /= grid.DX().Product();

            }
            alpha1_flux = convection_term + compression_term;

            alpha1_field(cell) = alpha1_previous_field_ghost(cell) - dt * alpha1_flux;
        }

    }

    return result;
}
//#####################################################################
// Get_MULES_Limiter
//#####################################################################
template<class TV> bool ADVECTION_VOF<TV>::
Get_MULES_Limiter(const GRID<TV>& grid,
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
        )
{
    bool result = true;

    MULES<TV> mules;
    result = mules.Limiter(grid,
            mpi_grid,
            rDeltaT,
            alpha1_field_ghost,
            alpha1_previous_field_ghost,
            face_velocity_flux_field,
            face_alpha1_upwind_flux_field,
            //face_alpha1_flux_field,
            face_alpha1_flux_corr_field,
            Sp,
            Su,
            alpha_max,
            alpha_min,
            boundary,
            lambda
            );

    return result;
}
//#####################################################################
template class ADVECTION_VOF<VECTOR<float,1> >;
template class ADVECTION_VOF<VECTOR<float,2> >;
template class ADVECTION_VOF<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_VOF<VECTOR<double,1> >;
template class ADVECTION_VOF<VECTOR<double,2> >;
template class ADVECTION_VOF<VECTOR<double,3> >;
#endif
