//#####################################################################
// Copyright 2017, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __UTILITIES__
#define __UTILITIES__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>

#include "./constants.h"

#include <algorithm>

namespace PhysBAM {

template<class TV> struct CUSTOM_PARSE_ARGS{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;

    TV_INT cell_count;
    TV domain_min_corner,domain_max_corner;
    int use_frames,adv_temporal_scheme,apply_advection,apply_narrow_band_dirichlet_bc_for_levelset,narrow_band_extent,bound_interfacial_cells_extent;
    int reinit_temporal_scheme,use_reinitialization,use_reinitialization_at_start,reinit_iterations,root_finding_method,max_iterations_root_finding;
    int reinitialization_criterion,reinitialization_frequency;
    int use_reinitialization_primer,max_reinit_primer_iterations;
    T epsilon_root_finding;
    int is_dt_fixed,is_dt_times_dx;
    T dt_given,t_start,t_final,write_time_interval,log_data_time_interval;
    int solve_momentum,mom_convection_scheme,apply_viscosity,visc_approx_scheme,density_approx_for_visc,apply_surface_tension,solve_energy,mdot_compute_method,mdot_interface_compute_method,mdot_interface_location_compute_method,energy_temporal_scheme,energy_temperature_extrapolation_scheme,energy_convection_scheme;
    int energy_conv_discrete_treatment,energy_conv_interface_location_compute_method,energy_diff_discrete_treatment,interface_location_compute_method_diffusion,limit_normalized_interface_distance_diffusion;
    T rho_liquid,rho_vapor,mu_liquid,mu_vapor,sigma; //surface tension coefficient
    T k_liquid,k_vapor,cp_liquid,cp_vapor,h_lv;
    int initial_temperature_profile,object;
    T initial_temperature_liquid,initial_temperature_vapor;
    TV object_center,object_semi_principal_axes_length;
    T object_radius,object_radius2;
    int velocity;
    TV velocity_center,velocity_vector;
    T kappa_velocity,time_period;
    int error_band;
    int hermite_interpolation_form,number_of_ghost_cells,number_of_global_ghost_cells;
    bool reinitialized_at_start,adv_limit_grad,reinit_limit_grad,use_sharp_thermal_conductivity;
    T saturation_temperature;
    VECTOR<VECTOR<VECTOR<int,TV::dimension>,2>,TV::dimension> face_velocity_boundary_conditions;
    VECTOR<VECTOR<VECTOR<T,TV::dimension>,2>,TV::dimension> face_velocity_boundary_values;
    VECTOR<VECTOR<int,2>,TV::dimension> temperature_boundary_conditions;
    VECTOR<VECTOR<T,2>,TV::dimension> temperature_boundary_values;
    bool use_viscous_jump,is_phase_change;
    int mdot_resolve_proximity_approach,energy_conv_resolve_proximity_approach,energy_diff_resolve_proximity_approach;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_walls;

    // alpha1
    int alpha1_temporal_discretization;
    int solve_alpha1_convection_term;
    int solve_alpha1_compression_term;
    int alpha1_use_mules;
    T alpha1_convection_term_surface_interpolation_scheme_limiter_factor;
    T alpha1_compression_term_surface_interpolation_scheme_limiter_factor;
    int alpha1_convection_term_discretization;
    int alpha1_compression_term_discretization;
    int alpha1_convection_term_surface_interpolation;
    int alpha1_compression_term_surface_interpolation;
    int alpha1_velocity_convection_term_surface_interpolation;

    int compute_energy_conv_term;
    int compute_energy_diff_term;
    T mdot_interface_distance_factor_limit;
    T energy_conv_interface_distance_factor_limit;
    T energy_diff_interface_distance_factor_limit;
    T energy_diff_levelset_proximity_distance_factor;
    T energy_diff_discretization_approximation_order;
    T interface_distance_factor_limit_global;
    int proximity_scalars_treatment_when_cell_is_isolated;
    int compute_levelset_transport_everywhere;
    int gfm_jump_values_compute_method;
    int gfm_interface_identification_method;
    int gfm_compute_method_for_value_at_interface;
    int levelset_advection_algorithm;
    int levelset_advection_convection_discretization_scheme;
    int levelset_reinitialization_spatial_discretization_scheme;
    int global_temporal_scheme;

    CUSTOM_PARSE_ARGS(){cell_count=TV_INT(),domain_min_corner=TV(),domain_max_corner=TV(),use_frames=0,adv_temporal_scheme=0,apply_advection=0;
        apply_narrow_band_dirichlet_bc_for_levelset=1,narrow_band_extent=3,bound_interfacial_cells_extent=1;
        reinit_temporal_scheme=0,use_reinitialization=0,use_reinitialization_at_start=0,reinit_iterations=5,root_finding_method=0,max_iterations_root_finding=99;
        use_reinitialization_primer=1,max_reinit_primer_iterations=10,epsilon_root_finding=(T)1e-8,is_dt_fixed=1,is_dt_times_dx=1,dt_given=(T)1.,t_start=(T)0.,t_final=(T)4.,write_time_interval=(T)1,log_data_time_interval=(T)1;
        solve_momentum=0,mom_convection_scheme=5,apply_viscosity=1,visc_approx_scheme=1,density_approx_for_visc=1,apply_surface_tension=1,solve_energy=0,mdot_compute_method=2,mdot_interface_compute_method=1,mdot_interface_location_compute_method=1,energy_temporal_scheme=1,energy_temperature_extrapolation_scheme=2,energy_convection_scheme=5;
        energy_conv_discrete_treatment=1,energy_conv_interface_location_compute_method=1,energy_diff_discrete_treatment=1,interface_location_compute_method_diffusion=1,limit_normalized_interface_distance_diffusion=1,reinitialization_criterion=0,reinitialization_frequency=1;
        rho_liquid=(T)1000,rho_vapor=(T)100,mu_liquid=(T)10,mu_vapor=(T)1,sigma=(T)-24.5,k_liquid=(T).6,k_vapor=(T).026,cp_liquid=(T)4216,cp_vapor=(T)2034,h_lv=(T)2257000;
        initial_temperature_profile=0,object=0,initial_temperature_liquid=(T)363.15,initial_temperature_vapor=(T)383.15,object_center=TV(),object_semi_principal_axes_length=TV();
        object_radius=(T).15,object_radius2=(T).15,velocity=0,velocity_center=TV(),kappa_velocity=(T)2.,error_band=3,hermite_interpolation_form=1,number_of_ghost_cells=1;
        number_of_global_ghost_cells=3,reinitialized_at_start=false,adv_limit_grad=false,reinit_limit_grad=false,use_sharp_thermal_conductivity=false,saturation_temperature=(T)0.;
        for(int axis=1;axis<=TV::dimension;++axis) for(int axis_side=1;axis_side<=2;++axis_side) face_velocity_boundary_conditions[axis][axis_side]=TV_INT();
        for(int axis=1;axis<=TV::dimension;++axis) for(int axis_side=1;axis_side<=2;++axis_side) face_velocity_boundary_values[axis][axis_side]=TV();
        for(int axis=1;axis<=TV::dimension;++axis) for(int axis_side=1;axis_side<=2;++axis_side) temperature_boundary_conditions[axis][axis_side]=1;
        for(int axis=1;axis<=TV::dimension;++axis) for(int axis_side=1;axis_side<=2;++axis_side) temperature_boundary_values[axis][axis_side]=(T)0.;
        use_viscous_jump=false,is_phase_change=false;
        mdot_resolve_proximity_approach=1,energy_conv_resolve_proximity_approach=1,energy_diff_resolve_proximity_approach=1;
        for(int axis=1;axis<=TV::dimension;++axis) for(int axis_side=1;axis_side<=2;++axis_side) domain_walls[axis][axis_side]=false;

        alpha1_temporal_discretization = 1;
        solve_alpha1_convection_term = 1;
        solve_alpha1_compression_term = 1;
        alpha1_use_mules = 1;
        alpha1_convection_term_surface_interpolation_scheme_limiter_factor = (T)1;
        alpha1_compression_term_surface_interpolation_scheme_limiter_factor = (T)1;
        alpha1_convection_term_discretization = 1;
        alpha1_compression_term_discretization = 1;
        alpha1_convection_term_surface_interpolation = 5;
        alpha1_compression_term_surface_interpolation = 0;
        alpha1_velocity_convection_term_surface_interpolation = 1;

        compute_energy_conv_term = 1;
        compute_energy_diff_term = 1;
        mdot_interface_distance_factor_limit = (T)0.1;
        energy_conv_interface_distance_factor_limit = (T)0.1;
        energy_diff_interface_distance_factor_limit = (T)0.1;
        energy_diff_levelset_proximity_distance_factor = (T)1.0;
        energy_diff_discretization_approximation_order = 1;
        interface_distance_factor_limit_global = (T)0.1;
        proximity_scalars_treatment_when_cell_is_isolated = 1;
        compute_levelset_transport_everywhere = 1;
        gfm_jump_values_compute_method = 1;
        gfm_interface_identification_method = 1;
        gfm_compute_method_for_value_at_interface = 1;
        levelset_advection_algorithm = 1;
        levelset_advection_convection_discretization_scheme = 5;
        levelset_reinitialization_spatial_discretization_scheme = 1;
        global_temporal_scheme = 1;
    }
    ~CUSTOM_PARSE_ARGS(){};
};

template<class TV> struct CONVERGENCE_PARAMETERS_EXTRAPOLATION{
    typedef typename TV::SCALAR T;
    T error;
    int current_iteration;
    bool is_iterations_limit_reached,is_convergence_limit_reached,is_error_cannot_reach_limit;
    CONVERGENCE_PARAMETERS_EXTRAPOLATION(){error=(T)1.,current_iteration=0,is_iterations_limit_reached=false,is_convergence_limit_reached=false;
        is_error_cannot_reach_limit=false;// true if error matches error_previous in a limit
    };
    ~CONVERGENCE_PARAMETERS_EXTRAPOLATION(){};
};

template<typename T> static inline bool Is_Equal(const T a,const T b,const T small_num_in=SMALL_NUMBER)
{
    return abs(a-b)<small_num_in;
}

template<typename T> static inline int pos0(const T a)
{
    return (a >= 0) ? 1 : 0;
}

template<class T_GRID> static inline int Number_Of_Active_Cells(const T_GRID &grid,ARRAY<int,VECTOR<int,T_GRID::dimension> > &active_cell_array,int cell_value=1)
{
    typedef VECTOR<int,T_GRID::dimension> T_INDEX;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    int number_of_active_cells=0;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const T_INDEX& cell=iterator.Cell_Index();
        if(active_cell_array(cell)==cell_value) number_of_active_cells++;
    }
    return number_of_active_cells;
}

static inline int Get_Next_Axis(const int& axis_current,const int& dimension)
{
    PHYSBAM_ASSERT((axis_current>=1) && (axis_current<=dimension));
    return (axis_current<dimension)?axis_current+1:1;
}

template<class T> static inline T Theta(const T phi1,const T phi2)
{
    T numerator=(max(phi1,(T)0)+max(phi2,(T)0));
    T denominator=(abs(phi1)+abs(phi2));
    if(abs(denominator)) return (numerator/denominator);
    else return (T)0;
}

template<class T> static inline T Theta(const T phi1,const T phi2,const T phi3,const T phi4)
{
    T numerator=(max(phi1,(T)0)+max(phi2,(T)0)+max(phi3,(T)0)+max(phi4,(T)0));
    T denominator=(abs(phi1)+abs(phi2)+abs(phi3)+abs(phi4));
    if(abs(denominator)) return (numerator/denominator);
    else return (T)0;
}

template<class T> static inline T Mu(const T theta,const T mu_liquid, const T mu_vapor)
{
    PHYSBAM_ASSERT(theta>=(T).0 && theta<=(T)1.);
    if(theta==(T)1.) return mu_liquid;
    else if(theta==(T).0) return mu_vapor;
    else return (mu_vapor*mu_liquid)/((mu_vapor*theta)+(mu_liquid*((T)1.-theta)));
}

template<typename T> static inline T One_Dimensional_Height_Fraction(const T phi1,const T phi2)
{
    // Sussman, A sharp interface method for incompressible two-phase flows; JCP 221(2007)469-505 [pg:480]
    // Modified to use when phi<=0:liquid, phi>0:vapor
    // To be used for "k_face=k_liquid*1D_height_fraction+k_vapor*(1-1D_height_fraction)"
    // TODO phi=0 is included into liquid, which is not true for phase change

    if(phi1<=(T)0 && phi2<=(T)0) return (T)1;
    else if(phi1>(T)0 && phi2>(T)0) return (T)0;
    else{T numerator=-(min(phi1,(T)0)+min(phi2,(T)0));
        T denominator=abs(phi1)+abs(phi2);
        if(abs(denominator)) return (numerator/denominator);
        else return (T)0;}
}

template<typename T> static inline T Get_Normalized_Interface_Distance_Using_Adjacent_Levelsets(const T phi_left,const T phi_right)
{
    /**
     * Interface is split into `below result`*\DeltaX on the left and (1-`below result`)*\DeltaX on the right
     */
    return abs(phi_left)/(abs(phi_left)+abs(phi_right));
}

template<typename T> static inline T Linear_Extrapolated_Quantity_Adjacent_To_Interface_Left(const T scalar_interface,const T scalar_right,const T normalized_interface_distance)
{
    return (scalar_interface-scalar_right*normalized_interface_distance)/((T)1-normalized_interface_distance);
}

template<typename T> static inline T Linear_Extrapolated_Quantity_Adjacent_To_Interface_Right(const T scalar_interface,const T scalar_left,const T normalized_interface_distance)
{
    return (scalar_interface-scalar_left*((T)1-normalized_interface_distance))/normalized_interface_distance;
}

template<typename T> static inline T Get_Scalar_At_Interface(const T phi_left,const T phi_right,const T scalar_left,const T scalar_right)
{
    return (((scalar_left*abs(phi_right))+(scalar_right*abs(phi_left)))/(abs(phi_left)+abs(phi_right)));
}

template<typename T> static inline T Compute_Derivative_From_First_Principles(const T scalar_first,const T scalar_second,const T dx)
{
    return (scalar_second-scalar_first)/dx;
}

template<typename T> static inline T Compute_Gradient_Non_Uniform_Stencil_2O(const T scalar_left,const T scalar_middle,const T scalar_right,const T dx_left,const T dx_right)
{
    /**
     * Second order gradient
     */
    T one_by_dx_left=(T)1./dx_left,one_by_dx_right=(T)1./dx_right,one_by_dx_left_plus_dx_right=(T)1./(dx_left+dx_right);
    T result = ((-dx_right*dx_right*scalar_left)+((dx_right*dx_right-dx_left*dx_left)*scalar_middle)+(dx_left*dx_left*scalar_right))*one_by_dx_left*one_by_dx_right*one_by_dx_left_plus_dx_right;
    return result;
}

template<typename T> static inline T Compute_Gradient_Non_Uniform_Stencil_One_Sided_Left_2O(const T scalar_far_left,const T scalar_left,const T scalar_this,const T dx_far_left,const T dx_left)
{
    /**
     * Second order one-sided gradient approximation [using values to the left]
     */
    T one_by_dx_left=(T)1./dx_left,one_by_dx_far_left=(T)1./dx_far_left,one_by_dx_left_plus_dx_far_left=(T)1./(dx_left+dx_far_left);
    T dx_left_plus_dx_far_left=(dx_left+dx_far_left);
    T result = ((scalar_this*((dx_left_plus_dx_far_left*dx_left_plus_dx_far_left)-(dx_left*dx_left))*one_by_dx_left_plus_dx_far_left*one_by_dx_left*one_by_dx_far_left)-
                (scalar_left*dx_left_plus_dx_far_left*one_by_dx_left*one_by_dx_far_left)+
                (scalar_far_left*dx_left*one_by_dx_left_plus_dx_far_left*one_by_dx_far_left));
    return result;
}

template<typename T> static inline T Compute_Gradient_Non_Uniform_Stencil_One_Sided_Right_2O(const T scalar_this,const T scalar_right,const T scalar_far_right,const T dx_right,const T dx_far_right)
{
    /**
     * Second order one-sided gradient approximation [using values to the right]
     */
    T one_by_dx_right=(T)1./dx_right,one_by_dx_far_right=(T)1./dx_far_right,one_by_dx_right_plus_dx_far_right=(T)1./(dx_right+dx_far_right);
    T dx_right_plus_dx_far_right=(dx_right+dx_far_right);
    T result = (-(scalar_this*((dx_right_plus_dx_far_right*dx_right_plus_dx_far_right)-(dx_right*dx_right))*one_by_dx_right_plus_dx_far_right*one_by_dx_right*one_by_dx_far_right)+
                (scalar_right*dx_right_plus_dx_far_right*one_by_dx_right*one_by_dx_far_right)-
                (scalar_far_right*dx_right*one_by_dx_right_plus_dx_far_right*one_by_dx_far_right));
    return result;
}

template<typename T> static inline T Compute_Gradient_Non_Uniform_Stencil_Skewed_Left_2O_Not_Coinciding(const T scalar_far_left,const T scalar_left,const T scalar_right,const T h1,const T h2)
{
    /**
     * Stencil: (i-2)--h1--(i-1)--h1--(i)--h2--(i+1)
     * (d(scalar)/dx)_i
     */
    T one_by_denominator = (T)1. / (h1*(h1+h2)*(((T)2.*h1)+h2));
    T result = ( (scalar_right*((T)3.*h1*h1)) +
                 (scalar_left*( (h2*h2) - ((T)4.*h1*h1) )) +
                 (scalar_far_left*( (h1*h1) - (h2*h2) )) ) * one_by_denominator;
    return result;
}

template<typename T> static inline T Compute_Gradient_Non_Uniform_Stencil_Skewed_Right_2O_Not_Coinciding(const T scalar_left,const T scalar_right,const T scalar_far_right,const T h1,const T h2)
{
    /**
     * Stencil: (i-1)--h2--(i)--h1--(i+1)--h1--(i+2)
     * (d(scalar)/dx)_i
     */
    T one_by_denominator = (T)1. / (h1*(h1+h2)*(((T)2.*h1)+h2));
    T result = ( (-scalar_left*((T)3.*h1*h1)) +
                 (scalar_right*( -(h2*h2) + ((T)4.*h1*h1) )) +
                 (scalar_far_right*( -(h1*h1) + (h2*h2) )) ) * one_by_denominator;
    return result;
}

template<typename T> static inline T Compute_Gradient_Non_Uniform_Stencil_2O_With_Dirichlet_Condition_At_Node(const T scalar_left, const T scalar_right,
                                                                                                              const T dx_very_close, const T dx_normal)
{
    /**
     * Second order gradient simplified by imposing Dirichlet condition at nodal
     * location to the value at the interface
     */
    if (Is_Equal(scalar_left, scalar_right)) return (T)0;

    T one_over_dx_normal = (T)1. / dx_normal;
    T one_over_dx_normal_plus_dx_very_close = (T)1. / (dx_normal + dx_very_close);

    T result = (scalar_right - scalar_left) * (dx_very_close * one_over_dx_normal * one_over_dx_normal_plus_dx_very_close);

    return result;
}

template<typename T> static inline T Compute_Laplacian_Uniform_Stencil_2O(const T scalar_first,const T scalar_this,const T scalar_second,const T dx_first,const T dx_second)
{
    /**
     * Second order laplacian
     */
    T one_by_dx_first_dx_second=(T)1./(dx_first*dx_second);
    return (scalar_first-(T)2*scalar_this+scalar_second)*one_by_dx_first_dx_second;
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_1O(const T scalar_first,const T scalar_this,const T scalar_second,const T dx_first,const T dx_second)
{
    /**
     * First order laplacian
     * Gave first order results in 1D and slightly better in 2D with L1 norm
     */
    if (Is_Equal(scalar_first, scalar_this) && Is_Equal(scalar_this, scalar_second)) return (T)0;

    T one_by_dx_first = (T)1. / dx_first, one_by_dx_second = (T)1. / dx_second;
    T one_by_dx_first_plus_dx_second = (T)1. / (dx_first + dx_second);

    T result = (T)2. * dx_second * (scalar_first - scalar_this);
    result += ((T)2. * dx_first * (scalar_second - scalar_this));
    result *= (one_by_dx_first * one_by_dx_second * one_by_dx_first_plus_dx_second);

    return result;
    //return (((T)2.*scalar_first*one_by_dx_first*one_by_dx_first_plus_dx_second)-((T)2.*scalar_this*one_by_dx_first*one_by_dx_second)+((T)2.*scalar_second*one_by_dx_second*one_by_dx_first_plus_dx_second));
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_1O_One_Sided(const T scalar_this,const T scalar_side,const T scalar_far_side,const T dx_side,const T dx_far_side)
{
    /**
     * First order laplacian [one-sided]
     */
    if(Is_Equal(scalar_this,scalar_side) && Is_Equal(scalar_side,scalar_far_side)) return (T)0;
    T one_by_dx_side=(T)1./dx_side,one_by_dx_far_side=(T)1./dx_far_side,one_by_dx_side_plus_dx_far_side=(T)1./(dx_side+dx_far_side);
    return (T)2.*((scalar_this*one_by_dx_side*one_by_dx_side_plus_dx_far_side)-(scalar_side*one_by_dx_side*one_by_dx_far_side)+(scalar_far_side*one_by_dx_far_side*one_by_dx_side_plus_dx_far_side));
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_1O_One_Sided_With_Dirichlet_Condition_At_Node(
                                                                    const T scalar_first, const T scalar_second,
                                                                    const T dx_very_close, const T dx_normal)
{
    /**
     * Simplified first order laplacian [non-uniform stencil], by imposing
     * Dirichlet condition at nearest node
     */
    if (Is_Equal(scalar_first, scalar_second)) return (T)0;

    T one_over_dx_normal = (T)1. / dx_normal;
    T one_over_dx_normal_plus_dx_very_close = (T)1. / (dx_normal + dx_very_close);

    T result = (T)2. * (scalar_second - scalar_first) * one_over_dx_normal * one_over_dx_normal_plus_dx_very_close;

    return result;
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_Skewed_Second(
        const T scalar_this, const T scalar_first, const T scalar_second, const T scalar_second_second,
        const T h2, const T h1
    ) {
    /**
     * Laplacian second order using non-uniform stencil
     * Stencil arrangement: (i-1) --[h2]-- (i) --[h1]-- (i+1) --[h1]-- (i+2)
     */
    if (Is_Equal(scalar_this, scalar_first) &&
        Is_Equal(scalar_this, scalar_second) &&
        Is_Equal(scalar_this, scalar_second_second))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0;

    a = (-(T)3*h1 + h2) / (h1*h1*h2);
    b = (T)6 * h1 / ( ((T)2*h1*h1*h2) + ((T)3*h1*h2*h2) + (h2*h2*h2) );
    c = ((T)4*h1 - (T)2*h2) / (h1*h1*(h1 + h2));
    d = (-h1+h2) / (h1*h1*((T)2*h1 + h2));

    result = (a * scalar_this) + (b * scalar_first) + (c * scalar_second) + (d * scalar_second_second);

    return result;
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_Skewed_First(
        const T scalar_this, const T scalar_first, const T scalar_second, const T scalar_first_first,
        const T h1, const T h2
    ) {
    /**
     * Laplacian second order using non-uniform stencil
     * Stencil arrangement: (i-2) --[h1]-- (i-1) --[h1]-- (i) --[h2]-- (i+1)
     */
    if (Is_Equal(scalar_this, scalar_first) &&
        Is_Equal(scalar_this, scalar_second) &&
        Is_Equal(scalar_this, scalar_first_first))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0;

    a = (-(T)3*h1 + h2) / (h1*h1*h2);
    b = (-h1 + h2) / (h1*h1*((T)2*h1 + h2));
    c = ((T)4*h1 - (T)2*h2) / (h1*h1*(h1+h2));
    d = (T)6*h1 / (h2*(h1+h2)*((T)2*h1+h2));

    result = (a * scalar_this) + (b * scalar_first_first) + (c * scalar_first) + (d * scalar_second);

    return result;
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_One_Sided_Skewed_Second(
        const T scalar_this, const T scalar_second, const T scalar_second_second, const T scalar_second_second_second,
        const T h2, const T h1
    ) {
    /**
     * Laplacian second order using non-uniform stencil
     * Stencil arrangement: (i) --[h2]-- (i+1) --[h1]-- (i+2) --[h1]-- (i+3)
     */
    if (Is_Equal(scalar_this, scalar_second) &&
        Is_Equal(scalar_this, scalar_second_second) &&
        Is_Equal(scalar_this, scalar_second_second_second))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0;

    a = (T)6 / ((T)2*h1*h2 + h2*h2);
    b = -((T)3*h1 + (T)2*h2) / (h1*h1*h2);
    c = (T)4 / (h1*h1);
    d = - (h1 + (T)2*h2) / (h1*h1*((T)2*h1 + h2));

    result = (a * scalar_this) + (b * scalar_second) + (c * scalar_second_second) + (d * scalar_second_second_second);

    return result;
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_One_Sided_Skewed_First(
        const T scalar_this, const T scalar_first, const T scalar_first_first, const T scalar_first_first_first,
        const T h1, const T h2
    ) {
    /**
     * Laplacian second order using non-uniform stencil
     * Stencil arrangement: (i-3) --[h1]-- (i-2) --[h1]-- (i-1) --[h2]-- (i)
     */
    if (Is_Equal(scalar_this, scalar_first) &&
        Is_Equal(scalar_this, scalar_first_first) &&
        Is_Equal(scalar_this, scalar_first_first_first))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0;

    a = (T)6 / ((T)2*h1*h2 + h2*h2);
    b = -(h1 + (T)2*h2) / (h1*h1*((T)2*h1 + h2));
    c = (T)4 / (h1*h1);
    d = -((T)3*h1 + (T)2*h2) / (h1*h1*h2);

    result = (a * scalar_this) + (b * scalar_first_first_first) + (c * scalar_first_first) + (d * scalar_first);

    return result;
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_Skewed_Second_Not_Coinciding(
        const T scalar_first, const T scalar_second, const T scalar_second_second,
        const T scalar_second_second_second, const T h2, const T h1
    ) {
    /**
     * Laplacian second order using non-uniform stencil
     * Stencil arrangement: (i-1) --[h2]-- (i) --[h1]-- (i+1) --[h1]-- (i+2) --[h1]-- (i+3)
     */
    if (Is_Equal(scalar_first, scalar_second) &&
        Is_Equal(scalar_first, scalar_second_second) &&
        Is_Equal(scalar_first, scalar_second_second_second))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0;

    a = (T)12*h1 / ( (T)6*h1*h1*h1 + (T)11*h1*h1*h2 + (T)6*h1*h2*h2 + h2*h2*h2 );
    b = (-(T)5*h1 + h2) / (h1*h1*(h1 + h2));
    c = ((T)8*h1 - (T)2*h2) / (h1*h1*((T)2*h1 + h2));
    d = (-(T)3*h1 + h2) / (h1*h1*((T)3*h1 + h2));

    result = (a * scalar_first) + (b * scalar_second) + (c * scalar_second_second) + (d * scalar_second_second_second);

    return result;
}

template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_Skewed_First_Not_Coinciding(
        const T scalar_first, const T scalar_second, const T scalar_first_first,
        const T scalar_first_first_first, const T h1, const T h2
    ) {
    /**
     * Laplacian second order using non-uniform stencil
     * Stencil arrangement: (i-3) --[h1]-- (i-2) --[h1]-- (i-1) --[h1]-- (i) --[h2]-- (i+1)
     */
    if (Is_Equal(scalar_first, scalar_second) &&
        Is_Equal(scalar_first, scalar_first_first) &&
        Is_Equal(scalar_first, scalar_first_first_first))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0;

    a = (-(T)3*h1 + h2) / (h1*h1*((T)3*h1 + h2));
    b = ((T)8*h1 - (T)2*h2) / (h1*h1*((T)2*h1 + h2));
    c = (-(T)5*h1 + h2) / (h1*h1*(h1 + h2));
    d = (T)12*h1 / ((h1+h2)*((T)2*h1+h2)*((T)3*h1+h2));

    result = (a * scalar_first_first_first) + (b * scalar_first_first) + (c * scalar_first) + (d * scalar_second);

    return result;
}

// ----------------- Old Laplacian 2O codes ----------------- //
/**
 * Laplacian second order using non-uniform stencil
 * Stencil arrangement: (i-1) --[h1]-- (i) --[h2]-- (i+1) --[h3]-- (i+2)
 */
template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_Skewed_Second_Old(
        const T scalar_this, const T scalar_first, const T scalar_second, const T scalar_second_second,
        const T h1, const T h2, const T h3
    ) {
    if (Is_Equal(scalar_this, scalar_first) &&
        Is_Equal(scalar_this, scalar_second) &&
        Is_Equal(scalar_this, scalar_second_second))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0, tmp = (T)0;

    a = -( (T)2 * ((h3 * h1 * h1 * h1) + ((h1 + h2 + h3)* h2 * h2 * h2) - ((h1 + h2) * h2) + (h3 * h3 * h3)) );
    tmp = (h1*h1*h1) * ((h2+h3)*(h2*h2) - (h2*h2) + (h3*h3));
    tmp += (h1*h1) * ( ((h2+h3)*(h2*h2*h2)) - (h2*h2) + (h3*h3*h3) );
    tmp += h1 * ( (h2*h2*h2*h2) + (h3*h3) - (h2*h2*h2) + (h3*h3*h3) );
    a /= tmp;

    b = (T)2 * ( ((h2+h3)*(h2*h2*h2)) - (h2*h2) + (h3*h3*h3) );
    tmp = (h1*h1*h1) * ( ((h2+h3)*h2*h2) - (h2*h2) + (h3*h3) );
    tmp += (h1*h1) * ( ((h2+h3)*h2*h2*h2) - (h2*h2) + (h3*h3*h3) );
    tmp += h1 * ( (h2*h2*h2*h2) + (h3*h3) - (h2*h2*h2) + (h3*h3*h3) );
    b /= tmp;

    c = (T)2 * ( ((h2+h3)*h1*h1*h1) - (h1*h2) + (h3*h3*h3) );
    tmp = (h1*h2*h2*h2*h2) + (h3*h3) + ( h1*h1*h1*( ((h2+h3)*h2*h2) - (h2*h2) + (h3*h3) ) );
    tmp += (-h1*h2*h2*h2) + (h3*h3*h3) + (h1*h1*( ((h2+h3)*h2*h2*h2) - (h2*h2) + (h3*h3*h3) ));
    c /= tmp;

    d = -(T)2 * ( (h2*h1*h1*h1) - (h1*h2*h2*h2) );
    tmp = (h1*h2*h2*h2*h2) + (h3*h3) + (h1*h1*h1*( ((h2+h3)*h2*h2) - (h2*h2) + (h3*h3) ));
    tmp += (-h1*h2*h2*h2) + (h3*h3*h3) + (h1*h1*( ((h2+h3)*h2*h2*h2) - (h2*h2) + (h3*h3*h3) ));
    d /= tmp;

    result = (a * scalar_this) + (b * scalar_first) + (c * scalar_second) + (d * scalar_second_second);

    return result;
}

/**
 * Laplacian second order using non-uniform stencil
 * Stencil arrangement: (i-2) --[h1]-- (i-1) --[h2]-- (i) --[h3]-- (i+1)
 */
template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_Skewed_First_Old(
        const T scalar_this, const T scalar_first, const T scalar_second, const T scalar_first_first,
        const T h1, const T h2, const T h3
    ) {
    if (Is_Equal(scalar_this, scalar_first) &&
        Is_Equal(scalar_this, scalar_second) &&
        Is_Equal(scalar_this, scalar_first_first))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0, tmp = (T)0;

    a = -( (T)2 * ((h1 * h3 * h3 * h3) + ((h1 + h2 + h3)* h2 * h2 * h2) - ((h2 + h3) * h1) + (h2 * h2 * h2)) );
    tmp = (h2*h2*h2) * ((h1+h2)*(h3*h3) + (h2*h2) + (h3*h1));
    tmp += (h2*h2) * ( ((h1+h2)*(h3*h3*h3)) - (h3*h1) + (h2*h2*h2) );
    tmp += (-h2) * ( (h2*h2*h2*h3*h3) + (h1+h1) + (h2*h2*h3*h3*h3) );
    a /= tmp;

    b = (T)2 * ( ((h1+h2)*(h3*h3*h3)) - (h3*h1) + (h2*h2*h2) );
    tmp = (h2*h2*h2) * ( ((h1+h2)*h3*h3) + (h2*h2) + (h3*h1) );
    tmp += (h2*h2) * ( ((h1+h2)*h3*h3*h3) - (h3*h1) + (h2*h2*h2) );
    tmp += (-h2) * ( (h2*h2*h2*h3*h3) + (h1+h1) + (h2*h2*h3*h3*h3) );
    b /= tmp;

    c = (T)2 * ( ((h1+h2)*h2*h2*h2) - (h2*h1) + (h2*h2*h2) );
    tmp = -(h2*(h1+(h2*h2*h2*h3*h3)+h1+(h2*h2*h3*h3*h3))) + ( h2*h2*h2*( ((h1+h2)*h3*h3) + (h2*h2) + (h3*h1) ) );
    tmp += (h2*h2*( ((h1+h2)*h3*h3*h3) - (h3*h1) + (h2*h2*h2) ));
    c /= tmp;

    d = (T)2 * ( (h3*h2*h2*h2) - (h2*h3*h3*h3) );
    tmp = (-h2*( (h1+h1) + (h2*h2*h2*h3*h3) + (h2*h2*h3*h3*h3) )) + (h2*h2*h2*( ((h1+h2)*h3*h3) + (h2*h2) + (h3*h1) ));
    tmp += (h2*h2*( ((h1+h2)*h3*h3*h3) - (h3*h1) + (h2*h2*h2) ));
    d /= tmp;

    result = (a * scalar_this) + (b * scalar_first) + (c * scalar_second) + (d * scalar_first_first);

    return result;
}

/**
 * Laplacian second order using non-uniform stencil
 * Stencil arrangement: (i) --[h1]-- (i+1) --[h2]-- (i+2) --[h3]-- (i+3)
 */
template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_One_Sided_Skewed_Second_Old(
        const T scalar_this, const T scalar_second, const T scalar_second_second, const T scalar_second_second_second,
        const T h1, const T h2, const T h3
    ) {
    if (Is_Equal(scalar_this, scalar_second) &&
        Is_Equal(scalar_this, scalar_second_second) &&
        Is_Equal(scalar_this, scalar_second_second_second))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0, tmp = (T)0;

    a = ( (T)2 * ((h3 * h1 * h1 * h1) - ((h2 + h3)* h1) + (h2*h2*h2) + (h2*h1) + h2 + (h3 * h3 * h3)) );
    tmp = (-h1*h1) + (h2*h2*h2*h1) + h2 + (h3*h3);
    tmp += (-h1*h1*h1)*( ((h1+h2+h3)*h1) + (h2*h2) -((h1+h2)*h1) + h2 + (h3*h3) );
    tmp += (h1*h1) + (h2*h2*h1) + h2 + (h3*h3*h3);
    tmp += (h1*h1*( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) ));
    a /= tmp;

    b = (T)2 * ( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) );
    tmp = (-h1*h1) + (h2*h2*h2*h1) + h2 + (h3*h3) - ( (h1*h1*h1)*( ((h1+h2+h3)*h1) + (h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3) ) );
    tmp += (h1*h1) + (h2*h2*h1) + h2 + (h3*h3*h3) + ( (h1*h1)*( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) ) );
    b /= tmp;

    c = -(T)2 * ( ((h1+h2+h3)*h1*h1*h1) - (h1*h1) + h2 + (h3*h3*h3) );
    tmp = (-h1*h1) + (h2*h2*h2*h1) + h2 + (h3*h3) - ( h1*h1*h1*( ((h1+h2+h3)*h1) + (h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3) ) );
    tmp += (h1*h1) + (h2*h2*h1) + h2 + (h3*h3*h3) + (h1*h1*( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) ));
    c /= tmp;

    d = (T)2 * ( ((h1+h2)*h1*h1*h1) - (h1*h1) + (h2*h2*h2) );
    tmp = (-h1*h1) + (h2*h2*h2*h1) + h2 + (h3*h3) - (h1*h1*h1*( ((h1+h2+h3)*h1) + (h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3) ));
    tmp += (h1*h1) + (h2*h2*h1) + h2 + (h3*h3*h3) + (h1*h1*( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) ));
    d /= tmp;

    result = (a * scalar_this) + (b * scalar_second) + (c * scalar_second_second) + (d * scalar_second_second_second);

    return result;
}

/**
 * Laplacian second order using non-uniform stencil
 * Stencil arrangement: (i-3) --[h3]-- (i-2) --[h2]-- (i-1) --[h1]-- (i)
 */
template<typename T> static inline T Compute_Laplacian_Non_Uniform_Stencil_2O_One_Sided_Skewed_First_Old(
        const T scalar_this, const T scalar_first, const T scalar_first_first, const T scalar_first_first_first,
        const T h1, const T h2, const T h3
    ) {
    if (Is_Equal(scalar_this, scalar_first) &&
        Is_Equal(scalar_this, scalar_first_first) &&
        Is_Equal(scalar_this, scalar_first_first_first))
        return (T)0;

    T result = (T)0;
    T a = (T)0, b = (T)0, c = (T)0, d = (T)0, tmp = (T)0;

    a = ( (T)2 * ((h3 * h1 * h1 * h1) - ((h2 + h3)* h1) + (h2*h2*h2) + (h2*h1) + h2 + (h3 * h3 * h3)) );
    tmp = (-h1*h1) + (h2*h2*h2*h1) + h2 + (h3*h3);
    tmp += (-h1*h1*h1)*( ((h1+h2+h3)*h1) + (h2*h2) -((h1+h2)*h1) + h2 + (h3*h3) );
    tmp += (h1*h1) + (h2*h2*h1) + h2 + (h3*h3*h3);
    tmp += (h1*h1*( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) ));
    a /= tmp;

    b = (T)2 * ( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) );
    tmp = (-h1*h1) + (h2*h2*h2*h1) + h2 + (h3*h3) - ( (h1*h1*h1)*( ((h1+h2+h3)*h1) + (h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3) ) );
    tmp += (h1*h1) + (h2*h2*h1) + h2 + (h3*h3*h3) + ( (h1*h1)*( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) ) );
    b /= tmp;

    c = -(T)2 * ( ((h1+h2+h3)*h1*h1*h1) - (h1*h1) + h2 + (h3*h3*h3) );
    tmp = (-h1*h1) + (h2*h2*h2*h1) + h2 + (h3*h3) - ( h1*h1*h1*( ((h1+h2+h3)*h1) + (h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3) ) );
    tmp += (h1*h1) + (h2*h2*h1) + h2 + (h3*h3*h3) + (h1*h1*( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) ));
    c /= tmp;

    d = (T)2 * ( ((h1+h2)*h1*h1*h1) - (h1*h1) + (h2*h2*h2) );
    tmp = (-h1*h1) + (h2*h2*h2*h1) + h2 + (h3*h3) - (h1*h1*h1*( ((h1+h2+h3)*h1) + (h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3) ));
    tmp += (h1*h1) + (h2*h2*h1) + h2 + (h3*h3*h3) + (h1*h1*( ((h1+h2+h3)*h1) + (h2*h2*h2) - ((h1+h2)*h1) + h2 + (h3*h3*h3) ));
    d /= tmp;

    result = (a * scalar_this) + (b * scalar_first) + (c * scalar_first_first) + (d * scalar_first_first_first);

    return result;
}
// ----------------- Old Laplacian 2O codes ----------------- //


template<typename T> static inline bool Proximity_Warning(const T dx_computed,
                                                          const T dx_reference,
                                                          const T interface_distance_factor_limit)
{
    return dx_computed < (interface_distance_factor_limit * dx_reference);
}

}
#endif
