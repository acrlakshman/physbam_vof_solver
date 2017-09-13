//#####################################################################
// Copyright 2014-2015, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "PARSE_INPUT_FILE.h"

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARSE_INPUT_FILE<TV>::
PARSE_INPUT_FILE()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARSE_INPUT_FILE<TV>::
~PARSE_INPUT_FILE()
{}
//#####################################################################
// Parse_Input_File
//#####################################################################
template<class TV> bool PARSE_INPUT_FILE<TV>::
Parse_Input_File(const char* filename,int& argc,std::vector<std::string>& input_arguments)
{
    std::ifstream ifs(filename);
    while(!ifs.eof()){
        std::string str;
        getline(ifs,str);
        if(str[0]=='#' || str=="\n" || (int)str[0]==0){continue;}
        else{
            int pos=str.find(' ');
            if((int)std::string::npos==pos){return false;}
            std::string arg,val;
            arg = str.substr(0,pos);
            val = str.substr(pos+1,str.length()-pos);

            input_arguments.push_back(arg);
            input_arguments.push_back(val);
        }
    }
    ifs.close();

    argc=input_arguments.size();

    return true;
}
//#####################################################################
template class PARSE_INPUT_FILE<VECTOR<float,1> >;
template class PARSE_INPUT_FILE<VECTOR<float,2> >;
template class PARSE_INPUT_FILE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARSE_INPUT_FILE<VECTOR<double,1> >;
template class PARSE_INPUT_FILE<VECTOR<double,2> >;
template class PARSE_INPUT_FILE<VECTOR<double,3> >;
#endif
