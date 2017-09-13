//#####################################################################
// Copyright 2017, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TWO_PHASE_EXAMPLE
//#####################################################################
#ifndef __TWO_PHASE_EXAMPLE__
#define __TWO_PHASE_EXAMPLE__

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
#include "ADVECTION_VOF.h"
#include "FACE_CENTERED_OPERATION_UTILITIES.h"
#include "COMPUTE_PHI.h"
#include "ALPHA1_UTILITIES.h"
#include "SURFACE_AREA.h"
#include "GRADIENT_FV.h"
#include "SURFACE_INTERPOLATION.h"
#include "INTERPOLATION_WEIGHTS.h"
#include "Parallel_Computation/MPI_UNIFORM_COMMUNICATION_HELPER.h"
#include "Tools/utilities.h"
#include "Helpers/INITIALIZE_ALPHA1.h"
#include "Helpers/VELOCITY_FIELD.h"
#include "LEVELSET_REINITIALIZATION.h"

namespace PhysBAM {

template<class TV>
class TWO_PHASE_EXAMPLE:public EXAMPLE<TV>
{
    typedef EXAMPLE<TV> BASE;typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;
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
    using BASE::Write_Frame_Title;

    using BASE::output_directory;using BASE::first_frame;using BASE::restart;using BASE::stream_type;using BASE::parse_args;using BASE::write_substeps_level;using BASE::frame_title;using BASE::write_last_frame;
    using BASE::test_number;

    T_GRID grid;
    MPI_UNIFORM_GRID<T_GRID> *mpi_grid;

    T cfl,gravity;
    int scale,output_number,current_frame,time_step_count_for_reinitialization;
    bool trigger_reinitialization,write_debug_data,has_log_data_started;
    T_BOUNDARY_SCALAR boundary_scalar;
    T_BOUNDARY_SCALAR *boundary;
    T_BOUNDARY_VECTOR vector_boundary_scalar;
    T_BOUNDARY_VECTOR *vector_boundary;
    CUSTOM_PARSE_ARGS<TV> custom_parse_args;
    T global_time;

    T_ARRAYS_SCALAR alpha1_field;
    T_ARRAYS_SCALAR alpha1_initial_field;
    T_ARRAYS_SCALAR alpha1_previous_field_ghost;
    T_ARRAYS_SCALAR alpha1_current_field_ghost;
    T_ARRAYS_SCALAR curvature_field_ghost;
    T_ARRAYS_VECTOR grad_alpha1_field;
    T_ARRAYS_VECTOR grad_alpha1_previous_field_ghost;
    T_ARRAYS_VECTOR grad_alpha1_current_field_ghost;
    T_ARRAYS_VECTOR cell_centered_velocity_field;
    T_ARRAYS_VECTOR cell_centered_velocity_previous_field_ghost;
    T_FACE_ARRAYS_SCALAR surface_area_field;
    T_FACE_ARRAYS_SCALAR face_centered_velocity_field;
    T_FACE_ARRAYS_SCALAR face_centered_velocity_previous_field_ghost;
    T_FACE_ARRAYS_SCALAR face_velocity_flux_field;
    T_FACE_ARRAYS_SCALAR face_velocity_flux_previous_field;
    T_FACE_ARRAYS_SCALAR face_centered_dummy_field;
    T_FACE_ARRAYS_SCALAR face_alpha1_flux_field;
    T_FACE_ARRAYS_SCALAR face_centered_rho_phi_field;
    T_FACE_ARRAYS_SCALAR face_centered_rho_phi_previous_field;
    T_LEVELSET levelset_alpha1;

    ADVECTION_VOF<TV> advection_vof;
    VELOCITY_FIELD<TV> velocity_field;
    CELL_CENTERED_OPERATION_UTILITIES<TV> cell_centered_operation_utilities;
    FACE_CENTERED_OPERATION_UTILITIES<TV> face_centered_operation_utilities;
    COMPUTE_PHI<TV> compute_phi;
    GRADIENT_FV<TV> gradient_fv;
    SURFACE_INTERPOLATION<TV> surface_interpolation;
    ALPHA1_UTILITIES<TV> alpha1_utilities;
    LEVELSET_REINITIALIZATION<TV> *levelset_reinitialization;

    TWO_PHASE_EXAMPLE(const STREAM_TYPE& stream_type);
    virtual ~TWO_PHASE_EXAMPLE();

    int Number_Of_Active_Cells_Near_Interface(ARRAY<int,TV_INT> &active_cell_array,int cell_value=1)
    {
        int number_of_active_cells=0;
        int number_of_active_cells_local=Number_Of_Active_Cells<T_GRID>(grid,active_cell_array,cell_value);
        number_of_active_cells=(mpi_grid)?(int)mpi_grid->Reduce_Add((T)number_of_active_cells_local):number_of_active_cells_local;
        return number_of_active_cells;
    }
//#####################################################################
    virtual void Initialize_Grids();
    virtual void Set_Boundary_Conditions(const T dt,const T time) = 0;
    virtual void Postprocess_Substep(const T dt) {}
    virtual void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
//#####################################################################
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Log_Parameters() const PHYSBAM_OVERRIDE;
    void Register_Options() PHYSBAM_OVERRIDE;
    void Parse_Options() PHYSBAM_OVERRIDE;
    void Parse_Late_Options() PHYSBAM_OVERRIDE;
    void Initialize_Fields();
    void Initialize_Alpha1();
    void Write_All_Log_Data_Files();
    void Read_Output_Files(const int frame);
    void Limit_Dt(T& dt,const T time);
    void Allocate_Space();
    void Bind_Variables();
    void Write_Variables_To_Temporary_Buckets();
    void Advect_Alpha1(const T dt, const T time);
    void Reinitialize_Levelset();
    void Update_Cell_Centered_Velocity_Field(const T dt, const T time);
    void Compute_Rho_Phi();
    void Solve_Momentum_Predictor(const T dt,const T time);
    void Solve_Momentum_Projection(const T dt,const T time);
    void Compute_Alpha1_Curvature();
    void Log_Output_Data_To_Files(std::string filename,const typename TV::SCALAR data_) const;
    void Test_Functions();
//#####################################################################
};
}
#endif // __TWO_PHASE_EXAMPLE__
