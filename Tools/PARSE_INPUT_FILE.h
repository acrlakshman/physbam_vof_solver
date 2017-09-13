//#####################################################################
// Copyright 2014-2015, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARSE_INPUT_FILE
//#####################################################################
#ifndef __PARSE_INPUT_FILE__
#define __PARSE_INPUT_FILE__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <fstream>
#include <vector>
#include <cstring>

namespace PhysBAM {

template<class TV>
class PARSE_INPUT_FILE
{
    typedef typename TV::SCALAR T;
  public:
    PARSE_INPUT_FILE();
    ~PARSE_INPUT_FILE();

    bool Parse_Input_File(const char* filename,int& argc,std::vector<std::string>& input_arguments);
};
}
#endif
