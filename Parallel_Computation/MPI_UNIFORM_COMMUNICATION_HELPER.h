//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_UNIFORM_COMMUNICATION_HELPER
//#####################################################################
#ifndef __MPI_UNIFORM_COMMUNICATION_HELPER__
#define __MPI_UNIFORM_COMMUNICATION_HELPER__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>

namespace PhysBAM{

template<class T_GRID,class T_MIXED_DERIVATIVES_ARRAYS>
void Recv_Boundary_Derivatives(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_MIXED_DERIVATIVES_ARRAYS& derivatives,const int tag,const MPI::Status& probe_status,int bandwidth);

template<class T_GRID,class T_MIXED_DERIVATIVES_ARRAYS>
void Exchange_Boundary_Derivatives(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_MIXED_DERIVATIVES_ARRAYS& derivatives,int bandwidth);

template<class T_GRID,class TV_ARRAYS>
void Recv_Boundary_Vectors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,TV_ARRAYS& vectors,const int tag,const MPI::Status& probe_status,int bandwidth);

template<class T_GRID,class TV_ARRAYS>
void Exchange_Boundary_Vectors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,TV_ARRAYS& vectors,int bandwidth);

template<class T_GRID,class T_ARRAYS>
void Recv_Boundary_Scalars(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_ARRAYS& scalars,const int tag,const MPI::Status& probe_status,int bandwidth);

template<class T_GRID,class T_ARRAYS>
void Exchange_Boundary_Scalars(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_ARRAYS& scalars,int bandwidth);

}
#endif
