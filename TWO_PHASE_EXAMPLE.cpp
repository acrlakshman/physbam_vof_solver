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
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <fstream>
//#include "Poisson_Solvers/POISSON_TEMPERATURE_UNIFORM.h"
#include "TWO_PHASE_EXAMPLE.h"
//#include "Tools/VELOCITY_UTILITIES.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TWO_PHASE_EXAMPLE<TV>::
TWO_PHASE_EXAMPLE(const STREAM_TYPE& stream_type):
    BASE(stream_type),grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),
    scale(0),time_step_count_for_reinitialization(0),
    trigger_reinitialization(false),has_log_data_started(false),
    boundary(0),vector_boundary(0),global_time((T)0),
    alpha1_field(grid.Domain_Indices()),
    alpha1_initial_field(grid.Domain_Indices()),
    alpha1_previous_field_ghost(grid.Domain_Indices()),
    alpha1_current_field_ghost(grid.Domain_Indices()),
    curvature_field_ghost(grid.Domain_Indices()),
    grad_alpha1_field(grid.Domain_Indices()),
    grad_alpha1_previous_field_ghost(grid.Domain_Indices()),
    grad_alpha1_current_field_ghost(grid.Domain_Indices()),
    cell_centered_velocity_field(grid.Domain_Indices()),
    cell_centered_velocity_previous_field_ghost(grid.Domain_Indices()),
    surface_area_field(grid, 0, false),
    face_centered_velocity_field(grid, 0, false),
    face_centered_velocity_previous_field_ghost(grid, 0, false),
    face_velocity_flux_field(grid, 0, false),
    face_velocity_flux_previous_field(grid, 0, false),
    face_centered_dummy_field(grid, 0, false),
    face_alpha1_flux_field(grid, 0, false),
    face_centered_rho_phi_field(grid, 0, false),
    face_centered_rho_phi_previous_field(grid, 0, false),
    levelset_alpha1(grid, *new ARRAY<T,TV_INT>()),
    levelset_reinitialization(0)
{
    LOG::Initialize_Logging(false,false,1<<30,true,1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TWO_PHASE_EXAMPLE<TV>::
~TWO_PHASE_EXAMPLE()
{
    if (levelset_reinitialization) delete levelset_reinitialization;
    if (mpi_grid) delete mpi_grid;
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    LOG::SCOPE scope("Write_Output_Files",1);
    std::string f=STRING_UTILITIES::string_sprintf("/%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Write_To_Text_File(output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),frame_title);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",grid);
    FILE_UTILITIES::Write_To_File(stream_type, output_directory + "/" + f + "/alpha1", alpha1_field);
    //FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",two_phase_momentum->face_velocities);
    //FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/centered_velocities", grad_alpha1_field);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/centered_velocities", cell_centered_velocity_field);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",levelset_alpha1);

    // Restart data
    if(frame==first_frame){FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
        FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/first_frame",frame,"\n");
        if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);}

    if(write_last_frame) FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=write_substeps_level){
        frame_title=title;
        std::stringstream ss;
        ss<<"Writing substep ["<<frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        LOG::filecout(ss.str());
        Write_Output_Files(++output_number);frame_title="";}
}
//#####################################################################
// Initialize_Grids
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Initialize_Grids()
{
}
//#####################################################################
// Initialize_Fields
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Initialize_Fields()
{
    advection_vof.Set_Options(custom_parse_args);
    SURFACE_AREA<TV> surface_area;
    surface_area_field.Resize(grid, 0, false);
    surface_area.Compute(grid, surface_area_field);

    velocity_field.Compute(grid, custom_parse_args, 1, 0, cell_centered_velocity_field);

    Initialize_Alpha1();

    // Initializing standard level set class object
    levelset_reinitialization = new LEVELSET_REINITIALIZATION<TV>(grid, mpi_grid, boundary, custom_parse_args);

}
//#####################################################################
// Initialize_Alpha1
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Initialize_Alpha1()
{
    INITIALIZE_ALPHA1<TV> initialize_alpha1;
    initialize_alpha1.Initialize(grid, custom_parse_args, alpha1_field);

    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();

        alpha1_initial_field(cell) = alpha1_field(cell);
        alpha1_previous_field_ghost(cell) = alpha1_field(cell);
        alpha1_current_field_ghost(cell) = alpha1_field(cell);
    }
}
//#####################################################################
// Write_All_Log_Data_Files
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Write_All_Log_Data_Files()
{
}
//#####################################################################
// Read_Output_Files
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("/%d",frame);
    // Restart data
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/common/grid",grid);
    if(mpi_grid) FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
}
//#####################################################################
// Limit_Dt
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Limit_Dt(T& dt,const T time)
{

    if (custom_parse_args.is_dt_fixed == 1) {

        if (custom_parse_args.is_dt_times_dx == 1) {
            dt = custom_parse_args.dt_given * grid.DX().Min();
        } else {
            dt = custom_parse_args.dt_given;
        }

    } else {

    }

}
//#####################################################################
// Log_Parameters
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Log_Parameters() const
{
    BASE::Log_Parameters();
    LOG::SCOPE scope("TWO_PHASE_EXAMPLE parameters");
    LOG::cout<<"CFL Number = "<<cfl<<std::endl;
    LOG::cout<<"Scale = "<<scale<<std::endl;
}
//#####################################################################
// Register_Options
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Register_Options()
{
    BASE::Register_Options();
    if(!parse_args) return;

    // Custom Stuff
    parse_args->Add_Double_Argument("-cfl",(T).5,"Compressible Flow CFL number.");
    parse_args->Add_Double_Argument("-gravity",(T).98,"Gravity.");
    parse_args->Add_Integer_Argument("-scale",0,"Resolution of the fluid mesh.");
    parse_args->Add_Option_Argument("-write_debug_data","Write debug data.");

    // More custom stuff
    parse_args->Add_Integer_Argument("-x_cells",40,"x_cells","Number of cells along x-direction");
    parse_args->Add_Integer_Argument("-y_cells",80,"y_cells","Number of cells along y-direction");
    parse_args->Add_Integer_Argument("-z_cells",40,"z_cells","Number of cells along z-direction");
    parse_args->Add_Double_Argument("-x_min",.0,"x_min","Minimum value for x-component");
    parse_args->Add_Double_Argument("-x_max",1.,"x_max","Maximum value for x-component");
    parse_args->Add_Double_Argument("-y_min",.0,"y_min","Minimum value for y-component");
    parse_args->Add_Double_Argument("-y_max",2.,"y_max","Maximum value for y-component");
    parse_args->Add_Double_Argument("-z_min",.0,"z_min","Minimum value for z-component");
    parse_args->Add_Double_Argument("-z_max",1.,"z_max","Maximum value for z-component");
    parse_args->Add_Integer_Argument("-use_frames",0,"use_frames","switch between frames and time steps");
    parse_args->Add_Integer_Argument("-adv_temporal_scheme",0,"adv_temporal_scheme","temporal scheme for advection");
    parse_args->Add_Integer_Argument("-apply_advection",1,"apply_advection","Advects levelset field");
    parse_args->Add_Integer_Argument("-adv_limit_grad",0,"adv_limit_grad","Implement gradient limiting in advection");
    parse_args->Add_Integer_Argument("-apply_narrow_band_dirichlet_bc_for_levelset",1,"apply_narrow_band_dirichlet_bc_for_levelset","Dirichlet bc is implemented for levelset and its gradient beyond a band");
    parse_args->Add_Integer_Argument("-narrow_band_extent",3,"narrow_band_extent","Narrow band extent");
    parse_args->Add_Integer_Argument("-bound_interfacial_cells_extent",1,"bound_interfacial_cells_extent","Bounding region for interfacial cells");
    parse_args->Add_Integer_Argument("-reinit_temporal_scheme",1,"reinit_temporal_scheme","temporal scheme for reinitialization");
    parse_args->Add_Integer_Argument("-reinit_limit_grad",0,"reinit_limit_grad","Implement gradient limiting in reinitialization");
    parse_args->Add_Integer_Argument("-use_reinitialization",1,"use_reinitialization","Reinitializes levelset field");
    parse_args->Add_Integer_Argument("-reinitialization_criterion",0,"reinitialization_criterion","Condition for reinitialization");
    parse_args->Add_Integer_Argument("-reinitialization_frequency",1,"reinitialization_frequency","Frequency of reinitialization (if taking from user)");
    parse_args->Add_Integer_Argument("-use_reinitialization_at_start",0,"use_reinitialization_at_start","Reinitializes levelset field before the start of sim");
    parse_args->Add_Integer_Argument("-reinit_iterations",4,"reinit_iterations","No. of times reinitialization equation is solved");
    parse_args->Add_Integer_Argument("-root_finding_method",0,"root_finding_method","Method to locate interface");
    parse_args->Add_Double_Argument("-epsilon_root_finding",1e-8,"epsilon_root_finding","Error tolerance to stop iterations");
    parse_args->Add_Integer_Argument("-max_iterations_root_finding",99,"max_iterations_root_finding","Maximum iterations allowed to find interface location");
    parse_args->Add_Integer_Argument("-use_reinitialization_primer",1,"use_reinitialization_primer","Uses standard reinitialization technique before triggering GALS reinitialization");
    parse_args->Add_Integer_Argument("-max_reinit_primer_iterations",10,"max_reinit_primer_iterations","Maximum standard reinitialization iterations");
    parse_args->Add_Integer_Argument("-is_dt_fixed",1,"is_dt_fixed","Use fixed time step");
    parse_args->Add_Integer_Argument("-is_dt_times_dx",1,"is_dt_times_dx","Is dt provided multiplied by dx?");
    parse_args->Add_Double_Argument("-dt_given",(T)1.,"dt_given","Provided time step");
    parse_args->Add_Double_Argument("-t_start",(T)0.,"final time");
    parse_args->Add_Double_Argument("-t_final",(T)3.,"final time");
    parse_args->Add_Double_Argument("-write_time_interval",(T)1.,"write_time_interval","Write interval for output");
    parse_args->Add_Double_Argument("-log_data_time_interval",(T)1.,"log_data_time_interval","Write interval for logging computed data (e.g. KE)");
    parse_args->Add_Integer_Argument("-solve_momentum",1,"solve_momentum","Solve two phase momentum equation");
    parse_args->Add_Integer_Argument("-mom_convection_scheme",5,"mom_convection_scheme","Discretization scheme for convection term in momemtum equation");
    parse_args->Add_Integer_Argument("-apply_viscosity",1,"apply_viscosity","Compute viscous forces");
    parse_args->Add_Integer_Argument("-visc_approx_scheme",1,"visc_approx_scheme","Viscous term approximation scheme");
    parse_args->Add_Integer_Argument("-density_approx_for_visc",1,"density_approx_for_visc","Density approximation in diffusion term");
    parse_args->Add_Integer_Argument("-apply_surface_tension",1,"apply_surface_tension","Include surface tension affects");
    parse_args->Add_Double_Argument("-rho_liquid",(T)1000,"rho_liquid","Density of liquid phase");
    parse_args->Add_Double_Argument("-rho_vapor",(T)100,"rho_vapor","Density of gas phase");
    parse_args->Add_Double_Argument("-mu_liquid",(T)10,"mu_liquid","Dynamic viscosity of liquid phase");
    parse_args->Add_Double_Argument("-mu_vapor",(T)1,"mu_vapor","Dynamic viscosity of gas phase");
    parse_args->Add_Double_Argument("-sigma",(T)24.5,"sigma","Surface tension coefficient");
    parse_args->Add_Integer_Argument("-solve_energy",1,"solve_energy","Solve two phase energy equation");
    parse_args->Add_Integer_Argument("-use_viscous_jump",0,"use_viscous_jump","Include viscosity jump in pressure");
    parse_args->Add_Integer_Argument("-is_phase_change",0,"is_phase_change","Solve with phase change");
    parse_args->Add_Double_Argument("-saturation_temperature",(T)0.,"saturation_temperature","Saturation temperature");
    parse_args->Add_Integer_Argument("-mdot_compute_method",2,"mdot_compute_method","Method to compute mass transfer rate");
    parse_args->Add_Integer_Argument("-mdot_interface_compute_method",1,"mdot_interface_compute_method","Method to compute mass transfer rate at interface");
    parse_args->Add_Integer_Argument("-mdot_interface_location_compute_method",1,"mdot_interface_location_compute_method","Method to compute interface location while computing mass transfer rate at interface");
    parse_args->Add_Integer_Argument("-energy_temporal_scheme",1,"energy_temporal_scheme","Temporal discretization scheme for energy equation");
    parse_args->Add_Double_Argument("-k_liquid",(T).6,"k_liquid","Thermal conductivity of liquid (W/m-K)");
    parse_args->Add_Double_Argument("-k_vapor",(T).026,"k_vapor","Thermal conductivity of vapor (W/m-K)");
    parse_args->Add_Double_Argument("-cp_liquid",(T)4216,"cp_liquid","Heat capacity per unit mass for liquid (J/kg-K)");
    parse_args->Add_Double_Argument("-cp_vapor",(T)2034,"cp_vapor","Heat capacity per unit mass for vapor (J/kg-K)");
    parse_args->Add_Double_Argument("-h_lv",(T)2257000,"h_lv","Latent heat of vaporization (Jg/kg)");
    parse_args->Add_Integer_Argument("-mdot_resolve_proximity_approach",1,"mdot_resolve_proximity_approach","Resolve proximity approach: while computing mdot");
    parse_args->Add_Integer_Argument("-energy_conv_resolve_proximity_approach",1,"energy_conv_resolve_proximity_approach","Resolve proximity approach: while computing energy convection term");
    parse_args->Add_Integer_Argument("-energy_diff_resolve_proximity_approach",1,"energy_diff_resolve_proximity_approach","Resolve proximity approach: while computing energy diffusiont term");

    parse_args->Add_Integer_Argument("-u_boundary_condition_left",1,"u_boundary_condition_left","Horizontal velocity boundary condition on left side");
    parse_args->Add_Integer_Argument("-u_boundary_condition_right",1,"u_boundary_condition_right","Horizontal velocity boundary condition on right side");
    parse_args->Add_Integer_Argument("-u_boundary_condition_bottom",1,"u_boundary_condition_bottom","Horizontal velocity boundary condition on bottom side");
    parse_args->Add_Integer_Argument("-u_boundary_condition_top",1,"u_boundary_condition_top","Horizontal velocity boundary condition on top side");
    parse_args->Add_Integer_Argument("-u_boundary_condition_back",1,"u_boundary_condition_back","Horizontal velocity boundary condition on back side");
    parse_args->Add_Integer_Argument("-u_boundary_condition_front",1,"u_boundary_condition_front","Horizontal velocity boundary condition on front side");
    parse_args->Add_Double_Argument("-u_boundary_value_left",(T)0.,"u_boundary_value_left","Horizontal velocity boundary value on left side");
    parse_args->Add_Double_Argument("-u_boundary_value_right",(T)0.,"u_boundary_value_right","Horizontal velocity boundary value on right side");
    parse_args->Add_Double_Argument("-u_boundary_value_bottom",(T)0.,"u_boundary_value_bottom","Horizontal velocity boundary value on bottom side");
    parse_args->Add_Double_Argument("-u_boundary_value_top",(T)0.,"u_boundary_value_top","Horizontal velocity boundary value on top side");
    parse_args->Add_Double_Argument("-u_boundary_value_back",(T)0.,"u_boundary_value_back","Horizontal velocity boundary value on back side");
    parse_args->Add_Double_Argument("-u_boundary_value_front",(T)0.,"u_boundary_value_front","Horizontal velocity boundary value on front side");

    parse_args->Add_Integer_Argument("-v_boundary_condition_left",1,"v_boundary_condition_left","Vertical velocity boundary condition on left side");
    parse_args->Add_Integer_Argument("-v_boundary_condition_right",1,"v_boundary_condition_right","Vertical velocity boundary condition on right side");
    parse_args->Add_Integer_Argument("-v_boundary_condition_bottom",1,"v_boundary_condition_bottom","Vertical velocity boundary condition on bottom side");
    parse_args->Add_Integer_Argument("-v_boundary_condition_top",1,"v_boundary_condition_top","Vertical velocity boundary condition on top side");
    parse_args->Add_Integer_Argument("-v_boundary_condition_back",1,"v_boundary_condition_back","Vertical velocity boundary condition on back side");
    parse_args->Add_Integer_Argument("-v_boundary_condition_front",1,"v_boundary_condition_front","Vertical velocity boundary condition on front side");
    parse_args->Add_Double_Argument("-v_boundary_value_left",(T)0.,"v_boundary_value_left","Vertical velocity boundary value on left side");
    parse_args->Add_Double_Argument("-v_boundary_value_right",(T)0.,"v_boundary_value_right","Vertical velocity boundary value on right side");
    parse_args->Add_Double_Argument("-v_boundary_value_bottom",(T)0.,"v_boundary_value_bottom","Vertical velocity boundary value on bottom side");
    parse_args->Add_Double_Argument("-v_boundary_value_top",(T)0.,"v_boundary_value_top","Vertical velocity boundary value on top side");
    parse_args->Add_Double_Argument("-v_boundary_value_back",(T)0.,"v_boundary_value_back","Vertical velocity boundary value on back side");
    parse_args->Add_Double_Argument("-v_boundary_value_front",(T)0.,"v_boundary_value_front","Vertical velocity boundary value on front side");

    parse_args->Add_Integer_Argument("-w_boundary_condition_left",1,"w_boundary_condition_left","Z velocity boundary condition on left side");
    parse_args->Add_Integer_Argument("-w_boundary_condition_right",1,"w_boundary_condition_right","Z velocity boundary condition on right side");
    parse_args->Add_Integer_Argument("-w_boundary_condition_bottom",1,"w_boundary_condition_bottom","Z velocity boundary condition on bottom side");
    parse_args->Add_Integer_Argument("-w_boundary_condition_top",1,"w_boundary_condition_top","Z velocity boundary condition on top side");
    parse_args->Add_Integer_Argument("-w_boundary_condition_back",1,"w_boundary_condition_back","Z velocity boundary condition on back side");
    parse_args->Add_Integer_Argument("-w_boundary_condition_front",1,"w_boundary_condition_front","Z velocity boundary condition on front side");
    parse_args->Add_Double_Argument("-w_boundary_value_left",(T)0.,"w_boundary_value_left","Z velocity boundary value on left side");
    parse_args->Add_Double_Argument("-w_boundary_value_right",(T)0.,"w_boundary_value_right","Z velocity boundary value on right side");
    parse_args->Add_Double_Argument("-w_boundary_value_bottom",(T)0.,"w_boundary_value_bottom","Z velocity boundary value on bottom side");
    parse_args->Add_Double_Argument("-w_boundary_value_top",(T)0.,"w_boundary_value_top","Z velocity boundary value on top side");
    parse_args->Add_Double_Argument("-w_boundary_value_back",(T)0.,"w_boundary_value_back","Z velocity boundary value on back side");
    parse_args->Add_Double_Argument("-w_boundary_value_front",(T)0.,"w_boundary_value_front","Z velocity boundary value on front side");

    parse_args->Add_Integer_Argument("-domain_wall_left",1,"domain_wall_left","Is wall on left side");
    parse_args->Add_Integer_Argument("-domain_wall_right",1,"domain_wall_right","Is wall on right side");
    parse_args->Add_Integer_Argument("-domain_wall_bottom",1,"domain_wall_bottom","Is wall on bottom side");
    parse_args->Add_Integer_Argument("-domain_wall_top",0,"domain_wall_top","Is wall on top side");
    parse_args->Add_Integer_Argument("-domain_wall_back",1,"domain_wall_back","Is wall on back side");
    parse_args->Add_Integer_Argument("-domain_wall_front",1,"domain_wall_front","Is wall on front side");

    parse_args->Add_Integer_Argument("-temperature_boundary_condition_left",1,"temperature_boundary_condition_left","Temperature boundary condition on left side");
    parse_args->Add_Integer_Argument("-temperature_boundary_condition_right",1,"temperature_boundary_condition_right","Temperature boundary condition on right side");
    parse_args->Add_Integer_Argument("-temperature_boundary_condition_bottom",1,"temperature_boundary_condition_bottom","Temperature boundary condition on bottom side");
    parse_args->Add_Integer_Argument("-temperature_boundary_condition_top",1,"temperature_boundary_condition_top","Temperature boundary condition on top side");
    parse_args->Add_Integer_Argument("-temperature_boundary_condition_back",1,"temperature_boundary_condition_back","Temperature boundary condition on back side");
    parse_args->Add_Integer_Argument("-temperature_boundary_condition_front",1,"temperature_boundary_condition_front","Temperature boundary condition on front side");
    parse_args->Add_Double_Argument("-temperature_boundary_value_left",(T)0.,"temperature_boundary_value_left","Temperature boundary value on left side");
    parse_args->Add_Double_Argument("-temperature_boundary_value_right",(T)0.,"temperature_boundary_value_right","Temperature boundary value on right side");
    parse_args->Add_Double_Argument("-temperature_boundary_value_bottom",(T)0.,"temperature_boundary_value_bottom","Temperature boundary value on bottom side");
    parse_args->Add_Double_Argument("-temperature_boundary_value_top",(T)0.,"temperature_boundary_value_top","Temperature boundary value on top side");
    parse_args->Add_Double_Argument("-temperature_boundary_value_back",(T)0.,"temperature_boundary_value_back","Temperature boundary value on back side");
    parse_args->Add_Double_Argument("-temperature_boundary_value_front",(T)0.,"temperature_boundary_value_front","Temperature boundary value on front side");
    parse_args->Add_Integer_Argument("-energy_temperature_extrapolation_scheme",2,"energy_temperature_extrapolation_scheme","Extrapolation scheme to employ for temperature, to use in convection and diffusion terms");
    parse_args->Add_Integer_Argument("-energy_convection_scheme",5,"energy_convection_scheme","Discretization scheme for convection term in energy equation");
    parse_args->Add_Integer_Argument("-energy_conv_discrete_treatment",1,"energy_conv_discrete_treatment","Method to discretize convection term in energy equation");
    parse_args->Add_Integer_Argument("-energy_conv_interface_location_compute_method",1,"energy_conv_interface_location_compute_method","Method to locate interface");
    parse_args->Add_Integer_Argument("-energy_diff_discrete_treatment",1,"energy_diff_discrete_treatment","Method to populate ghost temperature values");
    parse_args->Add_Integer_Argument("-interface_location_compute_method_diffusion",1,"interface_location_compute_method_diffusion","Method to identify interface between adjacent cell centers");
    parse_args->Add_Integer_Argument("-limit_normalized_interface_distance_diffusion",1,"limit_normalized_interface_distance_diffusion","Use constant extrapolation depending on the threshold");
    parse_args->Add_Integer_Argument("-use_sharp_thermal_conductivity",0,"use_sharp_thermal_conductivity","Sharp or weighted thermal conductivity");
    parse_args->Add_Integer_Argument("-initial_temperature_profile",0,"1D KH planar");
    parse_args->Add_Double_Argument("-initial_temperature_liquid",(T)363.15,"Initial temperature of liquid");
    parse_args->Add_Double_Argument("-initial_temperature_vapor",(T)383.15,"Initial temperature of vapor");
    parse_args->Add_Integer_Argument("-object",1,"circle");
    parse_args->Add_Double_Argument("-xo",(T).5,"x-coordinate of object center.");
    parse_args->Add_Double_Argument("-yo",(T).5,"y-coordinate of object center.");
    parse_args->Add_Double_Argument("-zo",(T).5,"z-coordinate of object center.");
    parse_args->Add_Double_Argument("-ro",(T).25,"Object radius.");
    parse_args->Add_Double_Argument("-p1",(T).75,"x-coordinate of object center.");
    parse_args->Add_Double_Argument("-p2",(T).3333,"y-coordinate of object center.");
    parse_args->Add_Double_Argument("-p3",(T).3333,"z-coordinate of object center.");
    parse_args->Add_Double_Argument("-pr",(T).25,"Object radius.");
    parse_args->Add_Integer_Argument("-velocity",3,"velocity","spatially constant");
    parse_args->Add_Double_Argument("-xc",(T).5,"x-coordinate of velocity center.");
    parse_args->Add_Double_Argument("-yc",(T).5,"y-coordinate of velocity center.");
    parse_args->Add_Double_Argument("-zc",(T).5,"z-coordinate of velocity center.");
    parse_args->Add_Double_Argument("-u_given",(T).0,"x-component of spatially constant velocity vector.");
    parse_args->Add_Double_Argument("-v_given",(T).0,"y-component of spatially constant velocity vector.");
    parse_args->Add_Double_Argument("-w_given",(T).0,"z-component of spatially constant velocity vector.");
    parse_args->Add_Double_Argument("-kappa_velocity",(T)2.,"kappa_velocity","Velocity deformation factor");
    parse_args->Add_Double_Argument("-T",(T)2.,"Time period","Time period for single vortex test case");
    parse_args->Add_Integer_Argument("-error_band",3,"error_band","Band with respect to dx with in which error is computed");
    parse_args->Add_Integer_Argument("-hermite_interpolation_form",1,"hermite_interpolation_form","Choose between B-form and standard form");
    parse_args->Add_Integer_Argument("-number_of_ghost_cells",3,"number_of_ghost_cells","Number of ghost cells for MPI communication");
    parse_args->Add_Integer_Argument("-number_of_global_ghost_cells",3,"number_of_global_ghost_cells","Number of ghost cells for variables in general");
    parse_args->Add_Double_Argument("-mdot_interface_distance_factor_limit", (T)0.1, "mdot_interface_distance_factor_limit",
                                    "mdot equation - interface distance factor limit to impose isothermal or appropriate approximations");
    parse_args->Add_Double_Argument("-energy_conv_interface_distance_factor_limit", (T)0.1, "energy_conv_interface_distance_factor_limit",
                                    "energy convection term - interface distance factor limit to impose isothermal or appropriate approximations");
    parse_args->Add_Double_Argument("-energy_diff_interface_distance_factor_limit", (T)0.1, "energy_diff_interface_distance_factor_limit",
                                    "energy diffusion term - interface distance factor limit to impose isothermal or appropriate approximations");
    parse_args->Add_Double_Argument("-energy_diff_levelset_proximity_distance_factor", (T)0.1, "energy_diff_levelset_proximity_distance_factor",
                                     "energy diffusion term - interface proximity limit to impose isothermal or appropriate approximations");
    parse_args->Add_Integer_Argument("-energy_diff_discretization_approximation_order", 1, "energy_diff_discretization_approximation_order",
                                    "energy diffusion term - discretization approximation order");
    parse_args->Add_Double_Argument("-interface_distance_factor_limit_global", (T)0.1, "interface_distance_factor_limit_global",
                                    "Interface distance factor limit to impose isothermal or appropriate approximations");
    parse_args->Add_Integer_Argument("-proximity_scalars_treatment_when_cell_is_isolated", 1, "proximity_scalars_treatment_when_cell_is_isolated",
                                     "How to generate ghost quantities, when interface is on both sides");

    parse_args->Add_Integer_Argument("-alpha1_temporal_discretization", 1, "alpha1_temporal_discretization", "alpha1 temporal discretization");
    parse_args->Add_Integer_Argument("-solve_alpha1_convection_term", 1, "solve_alpha1_convection_term", "solve alpha1 convection term");
    parse_args->Add_Integer_Argument("-solve_alpha1_compression_term", 1, "solve_alpha1_compression_term", "solve alpha1 compression term");
    parse_args->Add_Integer_Argument("-alpha1_use_mules", 1, "alpha1_use_mules", "MULES");
    parse_args->Add_Double_Argument("-alpha1_convection_term_surface_interpolation_scheme_limiter_factor", 1,
            "alpha1_convection_term_surface_interpolation_scheme_limiter_factor",
            "alpha1 convection term surface interpolation scheme limiter factor");
    parse_args->Add_Double_Argument("-alpha1_compression_term_surface_interpolation_scheme_limiter_factor", 1,
            "alpha1_compression_term_surface_interpolation_scheme_limiter_factor",
            "alpha1 compression term surface interpolation scheme limiter factor");
    parse_args->Add_Integer_Argument("-alpha1_convection_term_discretization", 1,
            "alpha1_convection_term_discretization", "alpha1 convection term discretization");
    parse_args->Add_Integer_Argument("-alpha1_compression_term_discretization", 1,
            "alpha1_compression_term_discretization", "alpha1 compression term discretization");
    parse_args->Add_Integer_Argument("-alpha1_convection_term_surface_interpolation", 5,
            "alpha1_convection_term_surface_interpolation", "alpha1 convection term surface interpolation");
    parse_args->Add_Integer_Argument("-alpha1_compression_term_surface_interpolation", 0,
            "alpha1_compression_term_surface_interpolation", "alpha1 compression term surface interpolation");
    parse_args->Add_Integer_Argument("-alpha1_velocity_convection_term_surface_interpolation", 1,
            "alpha1_velocity_convection_term_surface_interpolation", "alpha1 velocity convection term surface interpolation");

    parse_args->Add_Integer_Argument("-compute_energy_conv_term", 1, "compute_energy_conv_term", "Compute energy convection term");
    parse_args->Add_Integer_Argument("-compute_energy_diff_term", 1, "compute_energy_diff_term", "Compute energy diffusion term");
    parse_args->Add_Integer_Argument("-compute_levelset_transport_everywhere", 1, "compute_levelset_transport_everywhere", "Compute level set transport equation everywhere");
    parse_args->Add_Integer_Argument("-gfm_jump_values_compute_method", 1, "gfm_jump_values_compute_method", "Jump values compute method for GFM");
    parse_args->Add_Integer_Argument("-gfm_interface_identification_method", 1, "gfm_interface_identification_method", "Interface identification method for GFM");
    parse_args->Add_Integer_Argument("-gfm_compute_method_for_value_at_interface", 1, "gfm_compute_method_for_value_at_interface", "Compute method for value at interface (GFM)");
    parse_args->Add_Integer_Argument("-levelset_advection_algorithm", 1, "levelset_advection_algorithm", "Algorithm to advect level set");
    parse_args->Add_Integer_Argument("-levelset_advection_convection_discretization_scheme", 5, "levelset_advection_convection_discretization_scheme", "Discrete approximation for convection term in levelset advection equation");
    parse_args->Add_Integer_Argument("-levelset_reinitialization_spatial_discretization_scheme", 1, "levelset_reinitialization_spatial_discretization_scheme", "Spatial approximation for level set reinitialization equation");
    parse_args->Add_Integer_Argument("-global_temporal_scheme", 1, "global_temporal_scheme", "Global temporal scheme");
}
//#####################################################################
// Parse_Options
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Parse_Options()
{
    BASE::Parse_Options();
    scale=parse_args->Get_Integer_Value("-scale");
    cfl=(T)parse_args->Get_Double_Value("-cfl");
    gravity=(T)parse_args->Get_Double_Value("-gravity");
    write_debug_data=parse_args->Is_Value_Set("-write_debug_data");

    // More custom stuff
    T xo=(T)parse_args->Get_Double_Value("-xo");
    T yo=(T)parse_args->Get_Double_Value("-yo");
    T zo=(T)parse_args->Get_Double_Value("-zo");
    T ro=(T)parse_args->Get_Double_Value("-ro");
    T p1=(T)parse_args->Get_Double_Value("-p1");
    T p2=(T)parse_args->Get_Double_Value("-p2");
    T p3=(T)parse_args->Get_Double_Value("-p3");
    T pr=(T)parse_args->Get_Double_Value("-pr");
    T xc=(T)parse_args->Get_Double_Value("-xc");
    T yc=(T)parse_args->Get_Double_Value("-yc");
    T zc=(T)parse_args->Get_Double_Value("-zc");

    custom_parse_args.use_frames=parse_args->Get_Integer_Value("-use_frames");
    custom_parse_args.adv_temporal_scheme=parse_args->Get_Integer_Value("-adv_temporal_scheme");
    custom_parse_args.apply_advection=parse_args->Get_Integer_Value("-apply_advection");
    custom_parse_args.adv_limit_grad=(bool)parse_args->Get_Integer_Value("-adv_limit_grad");
    custom_parse_args.apply_narrow_band_dirichlet_bc_for_levelset=parse_args->Get_Integer_Value("-apply_narrow_band_dirichlet_bc_for_levelset");
    custom_parse_args.narrow_band_extent=parse_args->Get_Integer_Value("-narrow_band_extent");
    custom_parse_args.bound_interfacial_cells_extent=parse_args->Get_Integer_Value("-bound_interfacial_cells_extent");
    custom_parse_args.reinit_temporal_scheme=parse_args->Get_Integer_Value("-reinit_temporal_scheme");
    custom_parse_args.reinit_limit_grad=(bool)parse_args->Get_Integer_Value("-reinit_limit_grad");
    custom_parse_args.use_reinitialization=parse_args->Get_Integer_Value("-use_reinitialization");
    custom_parse_args.reinitialization_criterion=parse_args->Get_Integer_Value("-reinitialization_criterion");
    custom_parse_args.reinitialization_frequency=parse_args->Get_Integer_Value("-reinitialization_frequency");
    custom_parse_args.use_reinitialization_at_start=parse_args->Get_Integer_Value("-use_reinitialization_at_start");
    custom_parse_args.reinit_iterations=parse_args->Get_Integer_Value("-reinit_iterations");
    custom_parse_args.root_finding_method=parse_args->Get_Integer_Value("-root_finding_method");
    custom_parse_args.epsilon_root_finding=parse_args->Get_Double_Value("-epsilon_root_finding");
    custom_parse_args.max_iterations_root_finding=parse_args->Get_Integer_Value("-max_iterations_root_finding");
    custom_parse_args.use_reinitialization_primer=parse_args->Get_Integer_Value("-use_reinitialization_primer");
    custom_parse_args.max_reinit_primer_iterations=parse_args->Get_Integer_Value("-max_reinit_primer_iterations");
    custom_parse_args.is_dt_fixed=parse_args->Get_Integer_Value("-is_dt_fixed");
    custom_parse_args.is_dt_times_dx=parse_args->Get_Integer_Value("-is_dt_times_dx");
    custom_parse_args.dt_given=parse_args->Get_Double_Value("-dt_given");
    custom_parse_args.t_start=parse_args->Get_Double_Value("-t_start");
    custom_parse_args.t_final=parse_args->Get_Double_Value("-t_final");
    custom_parse_args.write_time_interval=parse_args->Get_Double_Value("-write_time_interval");
    custom_parse_args.log_data_time_interval=parse_args->Get_Double_Value("-log_data_time_interval");
    custom_parse_args.solve_momentum=parse_args->Get_Integer_Value("-solve_momentum");
    custom_parse_args.mom_convection_scheme=parse_args->Get_Integer_Value("-mom_convection_scheme");
    custom_parse_args.apply_viscosity=parse_args->Get_Integer_Value("-apply_viscosity");
    custom_parse_args.visc_approx_scheme=parse_args->Get_Integer_Value("-visc_approx_scheme");
    custom_parse_args.density_approx_for_visc=parse_args->Get_Integer_Value("-density_approx_for_visc");
    custom_parse_args.apply_surface_tension=parse_args->Get_Integer_Value("-apply_surface_tension");
    custom_parse_args.rho_liquid=(T)parse_args->Get_Double_Value("-rho_liquid");
    custom_parse_args.rho_vapor=(T)parse_args->Get_Double_Value("-rho_vapor");
    custom_parse_args.mu_liquid=(T)parse_args->Get_Double_Value("-mu_liquid");
    custom_parse_args.mu_vapor=(T)parse_args->Get_Double_Value("-mu_vapor");
    custom_parse_args.sigma=-(T)parse_args->Get_Double_Value("-sigma");
    custom_parse_args.solve_energy=parse_args->Get_Integer_Value("-solve_energy");
    custom_parse_args.use_viscous_jump=(bool)parse_args->Get_Integer_Value("-use_viscous_jump");
    custom_parse_args.is_phase_change=(bool)parse_args->Get_Integer_Value("-is_phase_change");
    custom_parse_args.saturation_temperature=(T)parse_args->Get_Double_Value("-saturation_temperature");
    custom_parse_args.mdot_compute_method=parse_args->Get_Integer_Value("-mdot_compute_method");
    custom_parse_args.mdot_interface_compute_method=parse_args->Get_Integer_Value("-mdot_interface_compute_method");
    custom_parse_args.mdot_interface_location_compute_method=parse_args->Get_Integer_Value("-mdot_interface_location_compute_method");
    custom_parse_args.energy_temporal_scheme=parse_args->Get_Integer_Value("-energy_temporal_scheme");
    custom_parse_args.k_liquid=parse_args->Get_Double_Value("-k_liquid");
    custom_parse_args.k_vapor=parse_args->Get_Double_Value("-k_vapor");
    custom_parse_args.cp_liquid=parse_args->Get_Double_Value("-cp_liquid");
    custom_parse_args.cp_vapor=parse_args->Get_Double_Value("-cp_vapor");
    custom_parse_args.h_lv=parse_args->Get_Double_Value("-h_lv");
    custom_parse_args.mdot_resolve_proximity_approach=parse_args->Get_Integer_Value("-mdot_resolve_proximity_approach");
    custom_parse_args.energy_conv_resolve_proximity_approach=parse_args->Get_Integer_Value("-energy_conv_resolve_proximity_approach");
    custom_parse_args.energy_diff_resolve_proximity_approach=parse_args->Get_Integer_Value("-energy_diff_resolve_proximity_approach");
    custom_parse_args.energy_temperature_extrapolation_scheme=parse_args->Get_Integer_Value("-energy_temperature_extrapolation_scheme");
    custom_parse_args.energy_convection_scheme=parse_args->Get_Integer_Value("-energy_convection_scheme");
    custom_parse_args.energy_conv_discrete_treatment=parse_args->Get_Integer_Value("-energy_conv_discrete_treatment");
    custom_parse_args.energy_conv_interface_location_compute_method=parse_args->Get_Integer_Value("-energy_conv_interface_location_compute_method");
    custom_parse_args.energy_diff_discrete_treatment=parse_args->Get_Integer_Value("-energy_diff_discrete_treatment");
    custom_parse_args.interface_location_compute_method_diffusion=parse_args->Get_Integer_Value("-interface_location_compute_method_diffusion");
    custom_parse_args.limit_normalized_interface_distance_diffusion=parse_args->Get_Integer_Value("-limit_normalized_interface_distance_diffusion");
    custom_parse_args.use_sharp_thermal_conductivity=(bool)parse_args->Get_Integer_Value("-use_sharp_thermal_conductivity");
    custom_parse_args.initial_temperature_profile=parse_args->Get_Integer_Value("-initial_temperature_profile");
    custom_parse_args.initial_temperature_liquid=parse_args->Get_Double_Value("-initial_temperature_liquid");
    custom_parse_args.initial_temperature_vapor=parse_args->Get_Double_Value("-initial_temperature_vapor");
    custom_parse_args.object=parse_args->Get_Integer_Value("-object");
    if(TV::dimension>=1){
        custom_parse_args.cell_count(1)=parse_args->Get_Integer_Value("-x_cells");
        custom_parse_args.domain_min_corner(1)=parse_args->Get_Double_Value("-x_min");
        custom_parse_args.domain_max_corner(1)=parse_args->Get_Double_Value("-x_max");
        custom_parse_args.object_center(1)=xo;
        custom_parse_args.object_semi_principal_axes_length(1)=p1;
        custom_parse_args.velocity_center(1)=xc;
        custom_parse_args.velocity_vector(1)=(T)parse_args->Get_Double_Value("-u_given");

        custom_parse_args.face_velocity_boundary_conditions[1][1][1]=parse_args->Get_Integer_Value("-u_boundary_condition_left");
        custom_parse_args.face_velocity_boundary_conditions[1][2][1]=parse_args->Get_Integer_Value("-u_boundary_condition_right");
        custom_parse_args.face_velocity_boundary_values[1][1][1]=(T)parse_args->Get_Double_Value("-u_boundary_value_left");
        custom_parse_args.face_velocity_boundary_values[1][2][1]=(T)parse_args->Get_Double_Value("-u_boundary_value_right");

        custom_parse_args.domain_walls[1][1]=(bool)parse_args->Get_Integer_Value("-domain_wall_left");
        custom_parse_args.domain_walls[1][2]=(bool)parse_args->Get_Integer_Value("-domain_wall_right");

        custom_parse_args.temperature_boundary_conditions[1][1]=parse_args->Get_Integer_Value("-temperature_boundary_condition_left");
        custom_parse_args.temperature_boundary_conditions[1][2]=parse_args->Get_Integer_Value("-temperature_boundary_condition_right");
        custom_parse_args.temperature_boundary_values[1][1]=(T)parse_args->Get_Double_Value("-temperature_boundary_value_left");
        custom_parse_args.temperature_boundary_values[1][2]=(T)parse_args->Get_Double_Value("-temperature_boundary_value_right");
    }
    if((int)TV::dimension>1){
        custom_parse_args.cell_count(2)=parse_args->Get_Integer_Value("-y_cells");
        custom_parse_args.domain_min_corner(2)=parse_args->Get_Double_Value("-y_min");
        custom_parse_args.domain_max_corner(2)=parse_args->Get_Double_Value("-y_max");
        custom_parse_args.object_center(2)=yo;
        custom_parse_args.object_semi_principal_axes_length(2)=p2;
        custom_parse_args.velocity_center(2)=yc;
        custom_parse_args.velocity_vector(2)=(T)parse_args->Get_Double_Value("-v_given");

        custom_parse_args.face_velocity_boundary_conditions[1][1][2]=parse_args->Get_Integer_Value("-v_boundary_condition_left");
        custom_parse_args.face_velocity_boundary_conditions[1][2][2]=parse_args->Get_Integer_Value("-v_boundary_condition_right");
        custom_parse_args.face_velocity_boundary_values[1][1][2]=(T)parse_args->Get_Double_Value("-v_boundary_value_left");
        custom_parse_args.face_velocity_boundary_values[1][2][2]=(T)parse_args->Get_Double_Value("-v_boundary_value_right");

        custom_parse_args.face_velocity_boundary_conditions[2][1][1]=parse_args->Get_Integer_Value("-u_boundary_condition_bottom");
        custom_parse_args.face_velocity_boundary_conditions[2][2][1]=parse_args->Get_Integer_Value("-u_boundary_condition_top");
        custom_parse_args.face_velocity_boundary_values[2][1][1]=(T)parse_args->Get_Double_Value("-u_boundary_value_bottom");
        custom_parse_args.face_velocity_boundary_values[2][2][1]=(T)parse_args->Get_Double_Value("-u_boundary_value_top");

        custom_parse_args.face_velocity_boundary_conditions[2][1][2]=parse_args->Get_Integer_Value("-v_boundary_condition_bottom");
        custom_parse_args.face_velocity_boundary_conditions[2][2][2]=parse_args->Get_Integer_Value("-v_boundary_condition_top");
        custom_parse_args.face_velocity_boundary_values[2][1][2]=(T)parse_args->Get_Double_Value("-v_boundary_value_bottom");
        custom_parse_args.face_velocity_boundary_values[2][2][2]=(T)parse_args->Get_Double_Value("-v_boundary_value_top");

        custom_parse_args.domain_walls[2][1]=(bool)parse_args->Get_Integer_Value("-domain_wall_bottom");
        custom_parse_args.domain_walls[2][2]=(bool)parse_args->Get_Integer_Value("-domain_wall_top");

        custom_parse_args.temperature_boundary_conditions[2][1]=parse_args->Get_Integer_Value("-temperature_boundary_condition_bottom");
        custom_parse_args.temperature_boundary_conditions[2][2]=parse_args->Get_Integer_Value("-temperature_boundary_condition_top");
        custom_parse_args.temperature_boundary_values[2][1]=(T)parse_args->Get_Double_Value("-temperature_boundary_value_bottom");
        custom_parse_args.temperature_boundary_values[2][2]=(T)parse_args->Get_Double_Value("-temperature_boundary_value_top");
    }
    if((int)TV::dimension>2){
        custom_parse_args.cell_count(3)=parse_args->Get_Integer_Value("-z_cells");
        custom_parse_args.domain_min_corner(3)=parse_args->Get_Double_Value("-z_min");
        custom_parse_args.domain_max_corner(3)=parse_args->Get_Double_Value("-z_max");
        custom_parse_args.object_center(3)=zo;
        custom_parse_args.object_semi_principal_axes_length(3)=p3;
        custom_parse_args.velocity_center(3)=zc;
        custom_parse_args.velocity_vector(3)=(T)parse_args->Get_Double_Value("-w_given");

        custom_parse_args.face_velocity_boundary_conditions[1][1][3]=parse_args->Get_Integer_Value("-w_boundary_condition_left");
        custom_parse_args.face_velocity_boundary_conditions[1][2][3]=parse_args->Get_Integer_Value("-w_boundary_condition_right");
        custom_parse_args.face_velocity_boundary_values[1][1][3]=(T)parse_args->Get_Double_Value("-w_boundary_value_left");
        custom_parse_args.face_velocity_boundary_values[1][2][3]=(T)parse_args->Get_Double_Value("-w_boundary_value_right");

        custom_parse_args.face_velocity_boundary_conditions[2][1][3]=parse_args->Get_Integer_Value("-w_boundary_condition_bottom");
        custom_parse_args.face_velocity_boundary_conditions[2][2][3]=parse_args->Get_Integer_Value("-w_boundary_condition_top");
        custom_parse_args.face_velocity_boundary_values[2][1][3]=(T)parse_args->Get_Double_Value("-w_boundary_value_bottom");
        custom_parse_args.face_velocity_boundary_values[2][2][3]=(T)parse_args->Get_Double_Value("-w_boundary_value_top");

        custom_parse_args.face_velocity_boundary_conditions[3][1][1]=parse_args->Get_Integer_Value("-u_boundary_condition_back");
        custom_parse_args.face_velocity_boundary_conditions[3][2][1]=parse_args->Get_Integer_Value("-u_boundary_condition_front");
        custom_parse_args.face_velocity_boundary_values[3][1][1]=(T)parse_args->Get_Double_Value("-u_boundary_value_back");
        custom_parse_args.face_velocity_boundary_values[3][2][1]=(T)parse_args->Get_Double_Value("-u_boundary_value_front");

        custom_parse_args.face_velocity_boundary_conditions[3][1][2]=parse_args->Get_Integer_Value("-v_boundary_condition_back");
        custom_parse_args.face_velocity_boundary_conditions[3][2][2]=parse_args->Get_Integer_Value("-v_boundary_condition_front");
        custom_parse_args.face_velocity_boundary_values[3][1][2]=(T)parse_args->Get_Double_Value("-v_boundary_value_back");
        custom_parse_args.face_velocity_boundary_values[3][2][2]=(T)parse_args->Get_Double_Value("-v_boundary_value_front");

        custom_parse_args.face_velocity_boundary_conditions[3][1][3]=parse_args->Get_Integer_Value("-w_boundary_condition_back");
        custom_parse_args.face_velocity_boundary_conditions[3][2][3]=parse_args->Get_Integer_Value("-w_boundary_condition_front");
        custom_parse_args.face_velocity_boundary_values[3][1][3]=(T)parse_args->Get_Double_Value("-w_boundary_value_back");
        custom_parse_args.face_velocity_boundary_values[3][2][3]=(T)parse_args->Get_Double_Value("-w_boundary_value_front");

        custom_parse_args.domain_walls[3][1]=(bool)parse_args->Get_Integer_Value("-domain_wall_back");
        custom_parse_args.domain_walls[3][2]=(bool)parse_args->Get_Integer_Value("-domain_wall_front");

        custom_parse_args.temperature_boundary_conditions[3][1]=parse_args->Get_Integer_Value("-temperature_boundary_condition_back");
        custom_parse_args.temperature_boundary_conditions[3][2]=parse_args->Get_Integer_Value("-temperature_boundary_condition_front");
        custom_parse_args.temperature_boundary_values[3][1]=(T)parse_args->Get_Double_Value("-temperature_boundary_value_back");
        custom_parse_args.temperature_boundary_values[3][2]=(T)parse_args->Get_Double_Value("-temperature_boundary_value_front");
    }
    custom_parse_args.object_radius=ro;
    custom_parse_args.object_radius2=pr;
    custom_parse_args.velocity=parse_args->Get_Integer_Value("-velocity");
    custom_parse_args.kappa_velocity=parse_args->Get_Double_Value("-kappa_velocity");
    custom_parse_args.time_period=parse_args->Get_Double_Value("-T");
    custom_parse_args.error_band=parse_args->Get_Integer_Value("-error_band");
    custom_parse_args.hermite_interpolation_form=parse_args->Get_Integer_Value("-hermite_interpolation_form");
    custom_parse_args.number_of_ghost_cells=parse_args->Get_Integer_Value("-number_of_ghost_cells");
    custom_parse_args.number_of_global_ghost_cells=parse_args->Get_Integer_Value("-number_of_global_ghost_cells");
    custom_parse_args.mdot_interface_distance_factor_limit = parse_args->Get_Double_Value("-mdot_interface_distance_factor_limit");
    custom_parse_args.energy_conv_interface_distance_factor_limit = parse_args->Get_Double_Value("-energy_conv_interface_distance_factor_limit");
    custom_parse_args.energy_diff_interface_distance_factor_limit = parse_args->Get_Double_Value("-energy_diff_interface_distance_factor_limit");
    custom_parse_args.energy_diff_levelset_proximity_distance_factor = parse_args->Get_Double_Value("-energy_diff_levelset_proximity_distance_factor");
    custom_parse_args.energy_diff_discretization_approximation_order = parse_args->Get_Integer_Value("-energy_diff_discretization_approximation_order");
    custom_parse_args.interface_distance_factor_limit_global = parse_args->Get_Double_Value("-interface_distance_factor_limit_global");
    custom_parse_args.proximity_scalars_treatment_when_cell_is_isolated = parse_args->Get_Integer_Value("-proximity_scalars_treatment_when_cell_is_isolated");

    custom_parse_args.alpha1_temporal_discretization = parse_args->Get_Integer_Value("-alpha1_temporal_discretization");
    custom_parse_args.solve_alpha1_convection_term = parse_args->Get_Integer_Value("-solve_alpha1_convection_term");
    custom_parse_args.solve_alpha1_compression_term = parse_args->Get_Integer_Value("-solve_alpha1_compression_term");
    custom_parse_args.alpha1_use_mules = parse_args->Get_Integer_Value("-alpha1_use_mules");
    custom_parse_args.alpha1_convection_term_surface_interpolation_scheme_limiter_factor =
        parse_args->Get_Double_Value("-alpha1_convection_term_surface_interpolation_scheme_limiter_factor");
    custom_parse_args.alpha1_compression_term_surface_interpolation_scheme_limiter_factor =
        parse_args->Get_Double_Value("-alpha1_compression_term_surface_interpolation_scheme_limiter_factor");
    custom_parse_args.alpha1_convection_term_discretization =
        parse_args->Get_Integer_Value("-alpha1_convection_term_discretization");
    custom_parse_args.alpha1_compression_term_discretization =
        parse_args->Get_Integer_Value("-alpha1_compression_term_discretization");
    custom_parse_args.alpha1_convection_term_surface_interpolation =
        parse_args->Get_Integer_Value("-alpha1_convection_term_surface_interpolation");
    custom_parse_args.alpha1_compression_term_surface_interpolation =
        parse_args->Get_Integer_Value("-alpha1_compression_term_surface_interpolation");
    custom_parse_args.alpha1_velocity_convection_term_surface_interpolation =
        parse_args->Get_Integer_Value("-alpha1_velocity_convection_term_surface_interpolation");

    custom_parse_args.compute_energy_conv_term = parse_args->Get_Integer_Value("-compute_energy_conv_term");
    custom_parse_args.compute_energy_diff_term = parse_args->Get_Integer_Value("-compute_energy_diff_term");
    custom_parse_args.compute_levelset_transport_everywhere = parse_args->Get_Integer_Value("-compute_levelset_transport_everywhere");
    custom_parse_args.gfm_jump_values_compute_method = parse_args->Get_Integer_Value("-gfm_jump_values_compute_method");
    custom_parse_args.gfm_interface_identification_method = parse_args->Get_Integer_Value("-gfm_interface_identification_method");
    custom_parse_args.gfm_compute_method_for_value_at_interface = parse_args->Get_Integer_Value("-gfm_compute_method_for_value_at_interface");
    custom_parse_args.levelset_advection_algorithm = parse_args->Get_Integer_Value("-levelset_advection_algorithm");
    custom_parse_args.levelset_advection_convection_discretization_scheme = parse_args->Get_Integer_Value("-levelset_advection_convection_discretization_scheme");
    custom_parse_args.levelset_reinitialization_spatial_discretization_scheme = parse_args->Get_Integer_Value("-levelset_reinitialization_spatial_discretization_scheme");
    custom_parse_args.global_temporal_scheme = parse_args->Get_Integer_Value("-global_temporal_scheme");
}
//#####################################################################
// Parse_Late_Options
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Parse_Late_Options()
{
    BASE::Parse_Late_Options();
    Initialize_Grids();
}
//#####################################################################
// Allocate_Space
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Allocate_Space()
{
    alpha1_field.Resize(grid.Domain_Indices());
    alpha1_initial_field.Resize(grid.Domain_Indices());
    alpha1_previous_field_ghost.Resize(grid.Domain_Indices(GLOBAL_NUM_CELLS_TRANSFER));
    alpha1_current_field_ghost.Resize(grid.Domain_Indices(GLOBAL_NUM_CELLS_TRANSFER));
    curvature_field_ghost.Resize(grid.Domain_Indices(GLOBAL_NUM_CELLS_TRANSFER));
    grad_alpha1_field.Resize(grid.Domain_Indices());
    grad_alpha1_previous_field_ghost.Resize(grid.Domain_Indices(GLOBAL_NUM_CELLS_TRANSFER));
    grad_alpha1_current_field_ghost.Resize(grid.Domain_Indices(GLOBAL_NUM_CELLS_TRANSFER));
    cell_centered_velocity_field.Resize(grid.Domain_Indices());
    cell_centered_velocity_previous_field_ghost.Resize(grid.Domain_Indices(GLOBAL_NUM_CELLS_TRANSFER));
    face_centered_velocity_field.Resize(grid, 0, false);
    face_centered_velocity_previous_field_ghost.Resize(grid, GLOBAL_NUM_CELLS_TRANSFER, false);
    face_velocity_flux_field.Resize(grid, 0, false);
    face_velocity_flux_previous_field.Resize(grid, 0, false);
    face_centered_dummy_field.Resize(grid, 0, false);
    face_alpha1_flux_field.Resize(grid, 0, false);
    face_centered_rho_phi_field.Resize(grid, 0, false);
    face_centered_rho_phi_previous_field.Resize(grid, 0, false);
    levelset_alpha1.phi.Resize(grid.Domain_Indices());
}
//#####################################################################
// Bind_Variables
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Bind_Variables()
{
}
//#####################################################################
// Write_Variables_To_Temporary_Buckets
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Write_Variables_To_Temporary_Buckets()
{
}
//#####################################################################
// Advect_Alpha1
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Advect_Alpha1(const T dt,const T time)
{
    T_ARRAYS_SCALAR alpha1_previous_field(grid.Domain_Indices());
    alpha1_previous_field = alpha1_field;

    // Transfer alpha1
    T_ARRAYS_SCALAR::Put(alpha1_previous_field,
            alpha1_previous_field_ghost);
    boundary->Fill_Ghost_Cells(grid,
            alpha1_previous_field,
            alpha1_previous_field_ghost,
            (T)0,
            (T)0,
            GLOBAL_NUM_CELLS_TRANSFER);
    if (mpi_grid) Exchange_Boundary_Scalars(*mpi_grid,
            alpha1_previous_field_ghost,
            GLOBAL_NUM_CELLS_TRANSFER);


    T_ARRAYS_VECTOR::Put(cell_centered_velocity_field,
            cell_centered_velocity_previous_field_ghost);
    vector_boundary->Fill_Ghost_Cells(grid,
            cell_centered_velocity_field,
            cell_centered_velocity_previous_field_ghost,
            (T)0,
            (T)0,
            GLOBAL_NUM_CELLS_TRANSFER);
    if (mpi_grid) Exchange_Boundary_Vectors(*mpi_grid,
            cell_centered_velocity_previous_field_ghost,
            GLOBAL_NUM_CELLS_TRANSFER);

    // Compute face_velocities
    T_FACE_ARRAYS_SCALAR face_centered_velocity_previous_field(grid, 0, false);
    face_centered_operation_utilities.Vector_Interpolate_Cell_To_Face_Linear(grid,
            cell_centered_velocity_previous_field_ghost,
            face_centered_velocity_previous_field);

    // Compute phi
    compute_phi.Compute(grid,
            face_centered_velocity_previous_field,
            surface_area_field,
            face_velocity_flux_previous_field);

    //LOG::cout << "alpha1_previous_field = " << alpha1_previous_field << std::endl;
    //LOG::cout << "alpha1_field = " << alpha1_field << std::endl;
    //LOG::cout << "alpha1_previous_field_ghost = " << alpha1_previous_field_ghost << std::endl;

    // Compute grad(alpha1)
    /*
    cell_centered_operation_utilities.Compute_Scalar_Gradient(grid,
        mpi_grid,
        alpha1_previous_field,
        alpha1_previous_field_ghost,
        2, //gradient_approximation_order,
        boundary,
        grad_alpha1_field);
        */
    if (not (surface_interpolation.Interpolate_Cell_To_Face_Scalar(grid,
                face_velocity_flux_previous_field,
                surface_area_field,
                alpha1_previous_field_ghost,
                grad_alpha1_previous_field_ghost, // DUMMY FOR LINEAR
                LINEAR,
                1,
                face_centered_dummy_field))) {
        PHYSBAM_ERROR("TWO_PHASE_EXAMPLE.cpp");
    }
    gradient_fv.Compute_For_Scalar_Field(grid,
        face_centered_dummy_field,
        surface_area_field,
        grad_alpha1_field);

    //LOG::cout << "grad_alpha1_field = " << grad_alpha1_field << std::endl;

    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const TV_INT cell = iterator.Cell_Index();
        alpha1_previous_field_ghost(cell) = alpha1_field(cell);
        grad_alpha1_previous_field_ghost(cell) = grad_alpha1_field(cell);
        cell_centered_velocity_previous_field_ghost(cell) = cell_centered_velocity_field(cell);
    }

    // Transfer values
    T_ARRAYS_VECTOR::Put(grad_alpha1_field,
            grad_alpha1_previous_field_ghost);
    vector_boundary->Fill_Ghost_Cells(grid,
            grad_alpha1_field,
            grad_alpha1_previous_field_ghost,
            (T)0,
            (T)0,
            GLOBAL_NUM_CELLS_TRANSFER);
    if (mpi_grid) Exchange_Boundary_Vectors(*mpi_grid,
            grad_alpha1_previous_field_ghost,
            GLOBAL_NUM_CELLS_TRANSFER);
    // End transfer

    advection_vof.Explicit_Solve(grid,
            mpi_grid,
            GLOBAL_NUM_CELLS_TRANSFER,
            dt,
            cell_centered_velocity_previous_field_ghost,
            face_centered_velocity_previous_field,
            face_velocity_flux_previous_field,
            T_FACE_LOOKUP(face_velocity_flux_previous_field),
            alpha1_previous_field_ghost,
            grad_alpha1_previous_field_ghost,
            surface_area_field,
            boundary,
            vector_boundary,
            face_alpha1_flux_field,
            alpha1_field);

    // Update current ghost field values
    // Update and transfer alpha1_current_field_ghost
    T_ARRAYS_SCALAR::Put(alpha1_field,
            alpha1_current_field_ghost);
    boundary->Fill_Ghost_Cells(grid,
            alpha1_field,
            alpha1_current_field_ghost,
            (T)0,
            (T)0,
            GLOBAL_NUM_CELLS_TRANSFER);
    if (mpi_grid) Exchange_Boundary_Scalars(*mpi_grid,
            alpha1_current_field_ghost,
            GLOBAL_NUM_CELLS_TRANSFER);

    // Compute grad(alpha1)
    // TODO: Replace below gradient calculation with GRADIENT_FV
    T_ARRAYS_VECTOR grad_alpha1_current_field(grid.Domain_Indices());
    /*
    cell_centered_operation_utilities.Compute_Scalar_Gradient(grid,
        mpi_grid,
        alpha1_field,
        alpha1_current_field_ghost,
        2, //gradient_approximation_order,
        boundary,
        grad_alpha1_current_field);
        */
    if (not (surface_interpolation.Interpolate_Cell_To_Face_Scalar(grid,
                face_velocity_flux_previous_field,
                surface_area_field,
                alpha1_current_field_ghost,
                grad_alpha1_previous_field_ghost, // DUMMY FOR LINEAR
                LINEAR,
                1,
                face_centered_dummy_field))) {
        PHYSBAM_ERROR("TWO_PHASE_EXAMPLE.cpp");
    }
    gradient_fv.Compute_For_Scalar_Field(grid,
        face_centered_dummy_field,
        surface_area_field,
        grad_alpha1_current_field);

    // Transfer grad_alpha1_current_field_ghost
    T_ARRAYS_VECTOR::Put(grad_alpha1_current_field,
            grad_alpha1_current_field_ghost);
    vector_boundary->Fill_Ghost_Cells(grid,
            grad_alpha1_current_field,
            grad_alpha1_current_field_ghost,
            (T)0,
            (T)0,
            GLOBAL_NUM_CELLS_TRANSFER);
    if (mpi_grid) Exchange_Boundary_Vectors(*mpi_grid,
            grad_alpha1_current_field_ghost,
            GLOBAL_NUM_CELLS_TRANSFER);
}
//#####################################################################
// Reinitialize_Levelset
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Reinitialize_Levelset()
{
    // Standard reinitialization
    int reinit_temporal_scheme = 2;
    if (custom_parse_args.reinit_temporal_scheme == 2)
        reinit_temporal_scheme = 1;
    else if (custom_parse_args.reinit_temporal_scheme == 3)
        reinit_temporal_scheme = 2;

    for (CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();
        // step function with alpha1=0.5 as interface
        levelset_alpha1.phi(cell) = tanh(alpha1_field(cell) - (T)0.5);
    }

    levelset_reinitialization->Solve_Levelset_Reinitialization(
            levelset_alpha1.phi,
            grid.DX().Min()/(T)4,
            reinit_temporal_scheme, // temporal scheme
            custom_parse_args.levelset_reinitialization_spatial_discretization_scheme,
            custom_parse_args.levelset_reinitialization_spatial_discretization_scheme,
            max(custom_parse_args.reinit_iterations, 10*grid.Counts().Max()),
            custom_parse_args.number_of_global_ghost_cells
            );
}
//#####################################################################
// Update_Cell_Centered_Velocity
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Update_Cell_Centered_Velocity_Field(const T dt,const T time)
{
    if (custom_parse_args.apply_advection == 1) {
        velocity_field.Compute(grid, custom_parse_args, dt, time, cell_centered_velocity_field);
    }
}
//#####################################################################
// Compute_Rho_Phi
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Compute_Rho_Phi()
{
    for (FACE_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        FACE_INDEX<TV::dimension> face = iterator.Full_Index();
        int axis = face.axis;
        TV_INT face_index = face.index;
        //TV_INT first_cell = iterator.First_Cell_Index();
        //TV_INT second_cell = iterator.Second_Cell_Index();

        face_centered_rho_phi_field(axis, face_index) =
            face_alpha1_flux_field(axis, face_index) * (custom_parse_args.rho_liquid - custom_parse_args.rho_vapor) +
            face_velocity_flux_previous_field(axis, face_index) * custom_parse_args.rho_vapor;
    }
}
//#####################################################################
// Solve_Momentum_Predictor
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Solve_Momentum_Predictor(const T dt,const T time)
{
}
//#####################################################################
// Solve_Momentum_Projection
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Solve_Momentum_Projection(const T dt,const T time)
{
}
//#####################################################################
// Compute_Alpha1_Curvature
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Compute_Alpha1_Curvature()
{
    alpha1_utilities.Compute_Curvature(grid,
                                       grad_alpha1_current_field_ghost,
                                       surface_area_field,
                                       curvature_field_ghost
                                      );
}
//#####################################################################
// Log_Output_Data_To_Files
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Log_Output_Data_To_Files(std::string filename,const typename TV::SCALAR data_) const
{
    std::ofstream ofs;
    const char *log_file_name=filename.c_str();
    if(has_log_data_started) ofs.open(log_file_name,std::ofstream::out | std::ofstream::app);
    else ofs.open(log_file_name,std::ofstream::out);
    ofs<<data_<<"\n";
    ofs.close();
}
//#####################################################################
// Test_Functions
//#####################################################################
template<class TV> void TWO_PHASE_EXAMPLE<TV>::
Test_Functions()
{
    LOG::cout << "Test Functions" << std::endl;
    SURFACE_AREA<TV> surface_area;
    INTERPOLATION_WEIGHTS<TV> interpolation_weights;

    T_FACE_ARRAYS_SCALAR surface_area_field;
    surface_area_field.Resize(grid, 0, false);

    T_FACE_ARRAYS_SCALAR geometric_weights_field;
    geometric_weights_field.Resize(grid, 0, false);

    surface_area.Compute(grid, surface_area_field);
    interpolation_weights.Compute_Geometric_Weights(grid, surface_area_field, geometric_weights_field);

    T l1_error = (T)0;
    int num_cells = 0;
    for(CELL_ITERATOR iterator(grid); iterator.Valid(); iterator.Next()) {
        const TV_INT& cell = iterator.Cell_Index();
        l1_error += fabs(alpha1_initial_field(cell) - alpha1_field(cell));
        num_cells++;
    }
    l1_error = mpi_grid ? mpi_grid->Reduce_Add(l1_error) : l1_error;
    num_cells = mpi_grid ? (int)mpi_grid->Reduce_Add((T)num_cells) : num_cells;

    l1_error /= num_cells;

    LOG::cout << "l1_error (alpha1) = " << l1_error << std::endl;

    // Curvature test for static circle
    T_ARRAYS_SCALAR curvature_field(grid.Domain_Indices());
    alpha1_utilities.Compute_Curvature(grid,
            grad_alpha1_current_field_ghost,
            surface_area_field,
            curvature_field
            );
    T_ARRAYS_SCALAR::Put(curvature_field, curvature_field_ghost);
}
//#####################################################################
template class TWO_PHASE_EXAMPLE<VECTOR<float,1> >;
template class TWO_PHASE_EXAMPLE<VECTOR<float,2> >;
template class TWO_PHASE_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TWO_PHASE_EXAMPLE<VECTOR<double,1> >;
template class TWO_PHASE_EXAMPLE<VECTOR<double,2> >;
template class TWO_PHASE_EXAMPLE<VECTOR<double,3> >;
#endif
