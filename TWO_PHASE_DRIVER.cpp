//#####################################################################
// Copyright 2017, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include "TWO_PHASE_DRIVER.h"
//#include "Tools/utilities.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TWO_PHASE_DRIVER<TV>::
TWO_PHASE_DRIVER(TWO_PHASE_EXAMPLE<TV>& example_input)
    : BASE(example_input),example(example_input)
{}
//#####################################################################
// Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR TWO_PHASE_DRIVER<TV>::
Compute_Dt(const T time,const T target_time)
{
    T dt=target_time-time;
    example.Limit_Dt(dt,time);
    return dt;
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void TWO_PHASE_DRIVER<TV>::
Initialize()
{
    BASE::Initialize();
    example.Parse_Late_Options();
    example.Log_Parameters();

    example.output_number=output_number;
    example.current_frame=current_frame;

    if(example.mpi_grid){example.mpi_grid->Initialize(example.custom_parse_args.domain_walls);
        example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.boundary_scalar);
        example.vector_boundary=new BOUNDARY_MPI<GRID<TV>,TV>(example.mpi_grid,example.vector_boundary_scalar);
    }
    else{example.boundary=&example.boundary_scalar;
        example.vector_boundary=&example.vector_boundary_scalar;
    }

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.custom_parse_args.domain_walls);
    example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.vector_boundary->Set_Constant_Extrapolation(domain_open_boundaries);

    example.Allocate_Space();

    // Initialize Fields
    example.Initialize_Fields();

    if (example.restart) {
        example.Read_Output_Files(example.restart_frame);
    }
}
//#####################################################################
// Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class TV> void TWO_PHASE_DRIVER<TV>::
Advance_One_Time_Step_Explicit_Part(const T dt,const T time)
{
    example.Advect_Alpha1(dt, time);
    example.Compute_Rho_Phi();
}
//#####################################################################
// Advance_One_Time_Step_Forces
//#####################################################################
template<class TV> void TWO_PHASE_DRIVER<TV>::
Advance_One_Time_Step_Forces(const T dt,const T time)
{
}
//#####################################################################
// Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class TV> void TWO_PHASE_DRIVER<TV>::
Advance_One_Time_Step_Implicit_Part(const T dt,const T time)
{
}
//#####################################################################
// Advance_One_Time_Step_Solve_Momentum
//#####################################################################
template<class TV> void TWO_PHASE_DRIVER<TV>::
Advance_One_Time_Step_Solve_Momentum(const T dt,const T time)
{
}
//#####################################################################
// Advance_One_Time_Step_Solve_Energy
//#####################################################################
template<class TV> void TWO_PHASE_DRIVER<TV>::
Advance_One_Time_Step_Solve_Energy(const T dt,const T time)
{
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void TWO_PHASE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep,1);
        T dt=Compute_Dt(time,target_time);
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);
        LOG::cout<<"Time step="<<dt<<std::endl;

        Advance_One_Time_Step_Explicit_Part(dt, time);

        example.Update_Cell_Centered_Velocity_Field(dt, time);

        // Write variables into temporarily introduced buckets
        example.Write_Variables_To_Temporary_Buckets();

        example.Postprocess_Substep(dt);

        if(!done) example.Write_Substep("END Substep",substep,0);
        time+=dt;
    }

}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void TWO_PHASE_DRIVER<TV>::
Simulate_To_Frame(const int target_frame)
{
    T write_time_interval=example.custom_parse_args.write_time_interval,next_write_time=write_time_interval;
    T log_data_time_interval=example.custom_parse_args.log_data_time_interval,next_log_data_time=log_data_time_interval;
    example.frame_title=STRING_UTILITIES::string_sprintf("Frame %d",example.current_frame);

    if (!example.restart) {
        Write_Output_Files(example.current_frame);
    }
    example.Write_All_Log_Data_Files();

    if(example.custom_parse_args.use_frames) while(example.current_frame<target_frame){
        LOG::SCOPE scope("FRAME","Frame %d",++example.current_frame,1);

        Advance_To_Target_Time(example.Time_At_Frame(example.current_frame));

        example.frame_title=STRING_UTILITIES::string_sprintf("Frame %d",example.current_frame);
        Write_Output_Files(++example.output_number);
        LOG::cout<<"TIME = "<<time<<std::endl;
    } else while(time<example.custom_parse_args.t_final) {
        LOG::SCOPE scope("FRAME","Frame %d",++example.current_frame,1);

        T dt=example.grid.DX().Min();
        if(example.custom_parse_args.is_dt_fixed){
            if(example.custom_parse_args.is_dt_times_dx) dt=example.custom_parse_args.dt_given*example.grid.dX(1);
            else dt=example.custom_parse_args.dt_given;}

        // Correct dt
        bool write_output=false;
        if(write_time_interval>(T)0){if ((time+dt > next_write_time || Is_Equal(time+dt, next_write_time, dt*(T)SMALL_NUMBER)) && time < next_write_time) {
            dt=next_write_time-time,write_output=true,next_write_time+=write_time_interval;
            next_write_time=min(next_write_time, example.custom_parse_args.t_final);
            if(example.mpi_grid){example.mpi_grid->Synchronize_Dt(dt);
                int local_write_output=(int)write_output;
                write_output=(bool)((int)example.mpi_grid->Reduce_Max(abs((T)local_write_output)));}}}
        else write_output=true;

        bool log_data_output = false;
        if (log_data_time_interval > (T)0) {
            if ((time+dt > next_log_data_time || Is_Equal(time+dt, next_log_data_time, dt*(T)SMALL_NUMBER)) && time < next_log_data_time) {
                dt = next_log_data_time - time;
                log_data_output = true;
                next_log_data_time += log_data_time_interval;

                if (example.mpi_grid) {
                    example.mpi_grid->Synchronize_Dt(dt);
                    int local_log_data_output = (int)log_data_output;
                    log_data_output = (bool)((int)example.mpi_grid->Reduce_Max(abs((T)local_log_data_output)));
                }
            }
        }

        Advance_To_Target_Time(time+dt);

        example.global_time=time;
        example.Test_Functions();

        example.frame_title=STRING_UTILITIES::string_sprintf("Frame %d",example.current_frame);
        if (write_output) {
            if (example.custom_parse_args.use_reinitialization)
                example.Reinitialize_Levelset();
            Write_Output_Files(++example.output_number);
        }
        if(log_data_output) {
            // Log
        }

        LOG::cout<<"TIME = "<<time<<std::endl;
        if(next_write_time>example.custom_parse_args.t_final) break;
    }
}
//#####################################################################
template class TWO_PHASE_DRIVER<VECTOR<float,1> >;
template class TWO_PHASE_DRIVER<VECTOR<float,2> >;
template class TWO_PHASE_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TWO_PHASE_DRIVER<VECTOR<double,1> >;
template class TWO_PHASE_DRIVER<VECTOR<double,2> >;
template class TWO_PHASE_DRIVER<VECTOR<double,3> >;
#endif
