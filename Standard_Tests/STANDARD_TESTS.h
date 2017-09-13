//#####################################################################
// Copyright 2014-2015, Mridul Aanjaneya, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include "../TWO_PHASE_EXAMPLE.h"
#include <fstream>

namespace PhysBAM {

template<class TV>
class STANDARD_TESTS : public TWO_PHASE_EXAMPLE<TV>
{
    typedef TWO_PHASE_EXAMPLE<TV> BASE;
    typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;typedef VECTOR<int,TV::dimension> TV_INT;

  public:
    using BASE::Write_Substep;

    using BASE::initial_time;using BASE::frame_rate;using BASE::last_frame;using BASE::grid;using BASE::mpi_grid;using BASE::scale;using BASE::cfl;
    using BASE::test_number;using BASE::output_directory;using BASE::parse_args;using BASE::boundary;using BASE::custom_parse_args;

    /****************************
     * example explanation:
     * **************************/

    STANDARD_TESTS(const STREAM_TYPE& stream_type)
        : BASE(stream_type)
    {}

//#####################################################################
    void Parse_Options() PHYSBAM_OVERRIDE
    {
        BASE::Parse_Options();
        output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d_Resolution_%d_CFL_%lf",test_number,scale,cfl,1);

        initial_time=(T)0.;
        TV_INT counts=TV_INT::All_Ones_Vector();

        grid.Initialize(custom_parse_args.cell_count,RANGE<TV>(custom_parse_args.domain_min_corner,custom_parse_args.domain_max_corner),true);
    }
//#####################################################################
    void Postprocess_Substep(const T dt)
    {LOG::cout<<"Time Step: "<<dt<<std::endl;}
//#####################################################################
    void Set_Boundary_Conditions(const T dt,const T time) PHYSBAM_OVERRIDE
    {
    }
//#####################################################################
};
}

#endif  // __STANDARD_TESTS__
