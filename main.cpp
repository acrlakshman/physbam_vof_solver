//#####################################################################
// Copyright 2017, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include "TWO_PHASE_DRIVER.h"
#include "TWO_PHASE_EXAMPLE.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Tools/PARSE_INPUT_FILE.h"
#include <vector>
#include <cstring>

using namespace PhysBAM;

template<class TV>
void Execute_Main_Program(int argc,char** argv)
{
    typedef typename TV::SCALAR RW;
    STREAM_TYPE stream_type((RW()));

    TWO_PHASE_EXAMPLE<TV>* example=new STANDARD_TESTS<TV>(stream_type);
    example->want_mpi_world=true;
    example->Parse(argc,argv);
    if(example->mpi_world->initialized){
        example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("/%d",(example->mpi_world->rank+1));}

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false,1);

    TWO_PHASE_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;
}

int main(int argc,char** argv)
{
    typedef double T;

    typedef VECTOR<T,2> TV;

    std::vector<std::string> input_arguments;
    if(argc>=2){PARSE_INPUT_FILE<TV> parse_input_file;
        if(parse_input_file.Parse_Input_File(argv[1],argc,input_arguments)){
            char *argv_tmp[argc+1];
            argv_tmp[0]=argv[0];
            for(int i=1;i<=argc;++i){
                char *arg=new char[input_arguments[i-1].length()+1];
                std::strcpy(arg,input_arguments[i-1].c_str());
                argv_tmp[i]=arg;
            }
            argc++;
            argv=argv_tmp;
            Execute_Main_Program<TV>(argc,argv);
        }
        else{LOG::cout<<"Invalid file/file-format"<<std::endl;}}
    else Execute_Main_Program<TV>(argc,argv);
    return 0;
}
