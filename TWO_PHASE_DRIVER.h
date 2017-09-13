//#####################################################################
// Copyright 2014-2015, Mridul Aanjaneya, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TWO_PHASE_DRIVER
//#####################################################################
#ifndef __TWO_PHASE_DRIVER__
#define __TWO_PHASE_DRIVER__

#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include "TWO_PHASE_EXAMPLE.h"

namespace PhysBAM {

template<class TV>
class TWO_PHASE_DRIVER:public DRIVER<TV>
{
    typedef DRIVER<TV> BASE;
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef GRID<TV> T_GRID;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;

  public:
    using BASE::output_number;using BASE::current_frame;using BASE::time;using BASE::Write_Output_Files;using BASE::Write_Substep;using BASE::Read_Time;

    TWO_PHASE_EXAMPLE<TV>& example;

    TWO_PHASE_DRIVER(TWO_PHASE_EXAMPLE<TV>& example_input);

//#####################################################################
    T Compute_Dt(const T time,const T target_time);
    void Initialize() PHYSBAM_OVERRIDE;
    void Advance_One_Time_Step_Explicit_Part(const T dt,const T time);
    void Advance_One_Time_Step_Forces(const T dt,const T time);
    void Advance_One_Time_Step_Implicit_Part(const T dt,const T time);
    void Advance_One_Time_Step_Solve_Momentum(const T dt,const T time);
    void Advance_One_Time_Step_Solve_Energy(const T dt,const T time);
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Simulate_To_Frame(const int frame) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif  // __TWO_PHASE_DRIVER__
