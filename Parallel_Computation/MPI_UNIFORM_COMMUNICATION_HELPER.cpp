//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include "MPI_UNIFORM_COMMUNICATION_HELPER.h"
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#endif
namespace PhysBAM{

#ifdef USE_MPI

//#####################################################################
// Function ISend_Derivatives
//#####################################################################
template<class T_GRID> MPI::Request
ISend_Derivatives(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const ARRAY<PAIR<typename T_GRID::VECTOR_T_MIXED_DERIVATIVE,typename T_GRID::VECTOR_INT> >& derivatives,const int destination_rank,
    const typename T_GRID::VECTOR_INT& destination_direction,const int tag,ARRAY<char>& buffer)
{
    typedef typename T_GRID::VECTOR_T_MIXED_DERIVATIVE TV_MIXED_DERIVATIVE;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    int position=0,pack_size=0;
    MPI::Comm& comm=*mpi_grid.comm;
    // compute pack size
    pack_size+=MPI_UTILITIES::Pack_Size(destination_direction,comm)+MPI_UTILITIES::Pack_Size(derivatives.m,comm);
    for(int i=1;i<=derivatives.m;i++){const TV_MIXED_DERIVATIVE& derivative=derivatives(i).x;const TV_INT& index=derivatives(i).y;
        pack_size+=MPI_UTILITIES::Pack_Size(derivative,comm);
        pack_size+=MPI_UTILITIES::Pack_Size(index,comm);}
    buffer.Resize(pack_size);
    MPI_UTILITIES::Pack(destination_direction,buffer,position,comm);
    MPI_UTILITIES::Pack(derivatives.m,buffer,position,comm);
    for(int i=1;i<=derivatives.m;++i){const TV_MIXED_DERIVATIVE& derivative=derivatives(i).x;const TV_INT& index=derivatives(i).y;
        MPI_UTILITIES::Pack(derivative,buffer,position,comm);
        MPI_UTILITIES::Pack(index,buffer,position,comm);}
    return comm.Isend(buffer.Get_Array_Pointer(),position,MPI::PACKED,destination_rank,tag);
}
//#####################################################################
// Function Recv_Boundary_Derivatives
//#####################################################################
template<class T_GRID,class T_MIXED_DERIVATIVES_ARRAYS>
void Recv_Boundary_Derivatives(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_MIXED_DERIVATIVES_ARRAYS& derivatives,const int tag,const MPI::Status& probe_status,int bandwidth)
{
    typedef typename T_GRID::VECTOR_T_MIXED_DERIVATIVE TV_MIXED_DERIVATIVE;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;
    MPI::Comm& comm=*mpi_grid.comm;
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    for(int i=1;i<=m;++i){
        TV_MIXED_DERIVATIVE derivative;MPI_UTILITIES::Unpack(derivative,buffer,position,comm);
        TV_INT index;MPI_UTILITIES::Unpack(index,buffer,position,comm);
        TV global_location=mpi_grid.global_grid.Center(index);
        TV_INT local_cell=mpi_grid.local_grid.Clamp_To_Cell(global_location,bandwidth);
        derivatives(local_cell)=derivative;}
}
//#####################################################################
// Function Exchange_Boundary_Derivatives
//#####################################################################
template<class T_GRID,class T_MIXED_DERIVATIVES_ARRAYS>
void Exchange_Boundary_Derivatives(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_MIXED_DERIVATIVES_ARRAYS& derivatives,int bandwidth)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    //typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_T_MIXED_DERIVATIVE TV_MIXED_DERIVATIVE;
    STATIC_ASSERT((IS_SAME<TV_MIXED_DERIVATIVE,typename T_MIXED_DERIVATIVES_ARRAYS::ELEMENT>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);

    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(T_GRID::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<TV_MIXED_DERIVATIVE,TV_INT> > > exchange_derivatives(T_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        for(CELL_ITERATOR iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            TV local_location=mpi_grid.local_grid.Center(cell);
            TV_INT global_cell_index=mpi_grid.global_grid.Clamp_To_Cell(local_location);
            exchange_derivatives(n).Append(PAIR<TV_MIXED_DERIVATIVE,TV_INT>(derivatives(cell),global_cell_index));}
        requests.Append(ISend_Derivatives(mpi_grid,exchange_derivatives(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}

    // probe and receive
    for(int message=1;message<=requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Boundary_Derivatives(mpi_grid,derivatives,tag,probe_status,bandwidth);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function ISend_Vectors
//#####################################################################
template<class T_GRID> MPI::Request
ISend_Vectors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const ARRAY<PAIR<typename T_GRID::VECTOR_T,typename T_GRID::VECTOR_INT> >& derivatives,const int destination_rank,
    const typename T_GRID::VECTOR_INT& destination_direction,const int tag,ARRAY<char>& buffer)
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    int position=0,pack_size=0;
    MPI::Comm& comm=*mpi_grid.comm;
    // compute pack size
    pack_size+=MPI_UTILITIES::Pack_Size(destination_direction,comm)+MPI_UTILITIES::Pack_Size(derivatives.m,comm);
    for(int i=1;i<=derivatives.m;i++){const TV& derivative=derivatives(i).x;const TV_INT& index=derivatives(i).y;
        pack_size+=MPI_UTILITIES::Pack_Size(derivative,comm);
        pack_size+=MPI_UTILITIES::Pack_Size(index,comm);}
    buffer.Resize(pack_size);
    MPI_UTILITIES::Pack(destination_direction,buffer,position,comm);
    MPI_UTILITIES::Pack(derivatives.m,buffer,position,comm);
    for(int i=1;i<=derivatives.m;++i){const TV& derivative=derivatives(i).x;const TV_INT& index=derivatives(i).y;
        MPI_UTILITIES::Pack(derivative,buffer,position,comm);
        MPI_UTILITIES::Pack(index,buffer,position,comm);}
    return comm.Isend(buffer.Get_Array_Pointer(),position,MPI::PACKED,destination_rank,tag);
}
//#####################################################################
// Function Recv_Boundary_Vectors
//#####################################################################
template<class T_GRID,class TV_ARRAYS>
void Recv_Boundary_Vectors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,TV_ARRAYS& vectors,const int tag,const MPI::Status& probe_status,int bandwidth)
{
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;
    MPI::Comm& comm=*mpi_grid.comm;
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    for(int i=1;i<=m;++i){
        TV derivative;MPI_UTILITIES::Unpack(derivative,buffer,position,comm);
        TV_INT index;MPI_UTILITIES::Unpack(index,buffer,position,comm);
        TV global_location=mpi_grid.global_grid.Center(index);
        TV_INT local_cell=mpi_grid.local_grid.Clamp_To_Cell(global_location,bandwidth);
        vectors(local_cell)=derivative;}
}
//#####################################################################
// Function Exchange_Boundary_Vectors
//#####################################################################
template<class T_GRID,class TV_ARRAYS>
void Exchange_Boundary_Vectors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,TV_ARRAYS& vectors,int bandwidth)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    //typedef typename T_GRID::SCALAR T;
    STATIC_ASSERT((IS_SAME<TV,typename TV_ARRAYS::ELEMENT>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);

    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(T_GRID::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<TV,TV_INT> > > exchange_vectors(T_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        for(CELL_ITERATOR iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            TV local_location=mpi_grid.local_grid.Center(cell);
            TV_INT global_cell_index=mpi_grid.global_grid.Clamp_To_Cell(local_location);
            exchange_vectors(n).Append(PAIR<TV,TV_INT>(vectors(cell),global_cell_index));}
        requests.Append(ISend_Vectors(mpi_grid,exchange_vectors(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}

    // probe and receive
    for(int message=1;message<=requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Boundary_Vectors(mpi_grid,vectors,tag,probe_status,bandwidth);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}
//#####################################################################
// Function ISend_Scalars
//#####################################################################
template<class T_GRID> MPI::Request
ISend_Scalars(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,const ARRAY<PAIR<typename T_GRID::SCALAR,typename T_GRID::VECTOR_INT> >& scalars,const int destination_rank,
    const typename T_GRID::VECTOR_INT& destination_direction,const int tag,ARRAY<char>& buffer)
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    int position=0,pack_size=0;
    MPI::Comm& comm=*mpi_grid.comm;
    // compute pack size
    pack_size+=MPI_UTILITIES::Pack_Size(destination_direction,comm)+MPI_UTILITIES::Pack_Size(scalars.m,comm);
    for(int i=1;i<=scalars.m;i++){const T& scalar=scalars(i).x;const TV_INT& index=scalars(i).y;
        pack_size+=MPI_UTILITIES::Pack_Size(scalar,comm);
        pack_size+=MPI_UTILITIES::Pack_Size(index,comm);}
    buffer.Resize(pack_size);
    MPI_UTILITIES::Pack(destination_direction,buffer,position,comm);
    MPI_UTILITIES::Pack(scalars.m,buffer,position,comm);
    for(int i=1;i<=scalars.m;++i){const T& scalar=scalars(i).x;const TV_INT& index=scalars(i).y;
        MPI_UTILITIES::Pack(scalar,buffer,position,comm);
        MPI_UTILITIES::Pack(index,buffer,position,comm);}
    return comm.Isend(buffer.Get_Array_Pointer(),position,MPI::PACKED,destination_rank,tag);
}
//#####################################################################
// Function Recv_Boundary_Scalars
//#####################################################################
template<class T_GRID,class T_ARRAYS>
void Recv_Boundary_Scalars(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_ARRAYS& scalars,const int tag,const MPI::Status& probe_status,int bandwidth)
{
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    MPI::Comm& comm=*mpi_grid.comm;
    ARRAY<char> buffer(probe_status.Get_count(MPI::PACKED));int position=0;
    MPI::Status status;
    comm.Recv(buffer.Get_Array_Pointer(),buffer.m,MPI::PACKED,probe_status.Get_source(),tag,status);
    TV_INT direction;MPI_UTILITIES::Unpack(direction,buffer,position,comm);
    TV wrap_offset=-mpi_grid.Wrap_Offset(-direction);
    int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);
    for(int i=1;i<=m;++i){
        T scalar;MPI_UTILITIES::Unpack(scalar,buffer,position,comm);
        TV_INT index;MPI_UTILITIES::Unpack(index,buffer,position,comm);
        TV global_location=mpi_grid.global_grid.Center(index);
        TV_INT local_cell=mpi_grid.local_grid.Clamp_To_Cell(global_location,bandwidth);
        scalars(local_cell)=scalar;}
}
//#####################################################################
// Function Exchange_Boundary_Scalars
//#####################################################################
template<class T_GRID,class T_ARRAYS>
void Exchange_Boundary_Scalars(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_ARRAYS& scalars,int bandwidth)
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::SCALAR T;
    STATIC_ASSERT((IS_SAME<T,typename T_ARRAYS::ELEMENT>::value));
    int tag=mpi_grid.Get_Unique_Tag();
    ARRAY<RANGE<TV_INT> > send_regions;
    RANGE<TV_INT> sentinels=RANGE<TV_INT>(TV_INT(),TV_INT());
    mpi_grid.Find_Boundary_Regions(send_regions,sentinels,false,RANGE<VECTOR<int,1> >(0,bandwidth-1),true,false);

    ARRAY<MPI::Request> requests;
    ARRAY<ARRAY<char> > buffers(T_GRID::number_of_one_ring_neighbors_per_cell);
    ARRAY<ARRAY<PAIR<T,TV_INT> > > exchange_scalars(T_GRID::number_of_one_ring_neighbors_per_cell);
    for(int n=1;n<=send_regions.m;n++) if(mpi_grid.all_neighbor_ranks(n)!=MPI::PROC_NULL){
        for(CELL_ITERATOR iterator(mpi_grid.local_grid,send_regions(n));iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            TV local_location=mpi_grid.local_grid.Center(cell);
            TV_INT global_cell_index=mpi_grid.global_grid.Clamp_To_Cell(local_location);
            exchange_scalars(n).Append(PAIR<T,TV_INT>(scalars(cell),global_cell_index));}
        requests.Append(ISend_Scalars(mpi_grid,exchange_scalars(n),mpi_grid.all_neighbor_ranks(n),mpi_grid.all_neighbor_directions(n),tag,buffers(n)));}

    // probe and receive
    for(int message=1;message<=requests.m;message++){
        MPI::Status probe_status;
        mpi_grid.comm->Probe(MPI::ANY_SOURCE,tag,probe_status);
        Recv_Boundary_Scalars(mpi_grid,scalars,tag,probe_status,bandwidth);}
    // wait for sends to complete
    MPI_UTILITIES::Wait_All(requests);
}


#else

//#####################################################################
// Function Recv_Boundary_Derivatives
//#####################################################################
template<class T_GRID,class T_MIXED_DERIVATIVES_ARRAYS>
void Recv_Boundary_Derivatives(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_MIXED_DERIVATIVES_ARRAYS& derivatives,const int tag,const MPI::Status& probe_status,int bandwidth)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Exchange_Boundary_Derivatives
//#####################################################################
template<class T_GRID,class T_MIXED_DERIVATIVES_ARRAYS>
void Exchange_Boundary_Derivatives(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_MIXED_DERIVATIVES_ARRAYS& derivatives,int bandwidth)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Recv_Boundary_Vectors
//#####################################################################
template<class T_GRID,class TV_ARRAYS>
void Recv_Boundary_Vectors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,TV_ARRAYS& vectors,const int tag,const MPI::Status& probe_status,int bandwidth)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Exchange_Boundary_Vectors
//#####################################################################
template<class T_GRID,class TV_ARRAYS>
void Exchange_Boundary_Vectors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,TV_ARRAYS& vectors,int bandwidth)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Recv_Boundary_Scalars
//#####################################################################
template<class T_GRID,class T_ARRAYS>
void Recv_Boundary_Scalars(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_ARRAYS& scalars,const int tag,const MPI::Status& probe_status,int bandwidth)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Exchange_Boundary_Scalars
//#####################################################################
template<class T_GRID,class T_ARRAYS>
void Exchange_Boundary_Scalars(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,T_ARRAYS& scalars,int bandwidth)
{PHYSBAM_NOT_IMPLEMENTED();}

#endif

//#####################################################################
#define INSTANTIATION_HELPER(T_GRID) \
    template void Exchange_Boundary_Derivatives(const MPI_UNIFORM_GRID<T_GRID >&,GRID_ARRAYS_POLICY<T_GRID >::ARRAYS_SCALAR::REBIND<T_GRID::VECTOR_T_MIXED_DERIVATIVE>::TYPE&,const int); \
    template void Exchange_Boundary_Vectors(const MPI_UNIFORM_GRID<T_GRID >&,GRID_ARRAYS_POLICY<T_GRID >::ARRAYS_SCALAR::REBIND<T_GRID::VECTOR_T>::TYPE&,const int); \
    template void Exchange_Boundary_Scalars(const MPI_UNIFORM_GRID<T_GRID >&,GRID_ARRAYS_POLICY<T_GRID >::ARRAYS_SCALAR::REBIND<T_GRID::SCALAR>::TYPE&,const int);
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(P(GRID<VECTOR<float,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,3> >));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(P(GRID<VECTOR<double,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,3> >));
#endif
}
