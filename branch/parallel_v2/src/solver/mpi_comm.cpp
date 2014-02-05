//==================================================================================================
// Name        : mpi_comm.cpp
// Authors     : Abhishek Y. Deshmukh
// Version     : 1.3
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the main functions to solve the fem problem.
//==================================================================================================

#include "mpi_comm.h"

//==================================================================================================
// MPI Communication
//==================================================================================================
void mpiComm::communicate(triMesh* mesh, double* M, double* RHS)
{    
    int nn = mesh->getNn();
    int conn_g, conn_l;
    int nprocb, source;

    // Allocate packets of send/recv buffer
    Sendbuf = new double[3*mesh->nbnmax]();

    // Prepare Sendbuffer
    int count = 0;
    int nbn_l = mesh->nbn;
    for(int i=0;i<mesh->nbn;i++){
	conn_g = mesh->procb[6*i+0];
	conn_l = mesh->getMapConn_G_L()[conn_g];
	Sendbuf[0+3*count] = conn_g;
	Sendbuf[1+3*count] = M[conn_l];
	Sendbuf[2+3*count] = RHS[conn_l];
	count++;
    }

	// Determine previous and next proc rank for communication
	int myrank = MPI::COMM_WORLD.Get_rank();
	int nprocs = MPI::COMM_WORLD.Get_size();
	int prev = (myrank+nprocs-1)%nprocs;
	int next = (myrank+1)%nprocs;

/*	if(myrank == 0){
		cout<<"nbn = "<<mesh->nbn<<endl;
		for(int i=0;i<mesh->nbn;i++)
			cout<<Sendbuf[3*i+0]<<"\t"<<Sendbuf[3*i+1]<<"\t"<<Sendbuf[3*i+2]<<endl;
	}
*/

	// Round-robin communication loop
	for(int iprocs=0;iprocs<nprocs-1;iprocs++){

		MPI::COMM_WORLD.Sendrecv_replace(Sendbuf, 3*mesh->nbnmax, MPI::DOUBLE, next, 0, prev, 0);
		MPI::COMM_WORLD.Sendrecv_replace(&nbn_l, 1, MPI::INT, next, 0, prev, 0);
		for(int i=0;i<nbn_l;i++){
			conn_l = mesh->getMapConn_G_L()[int(Sendbuf[3*i+0])];
			if(conn_l!=-1){
				M[conn_l] += Sendbuf[3*i+1];
				RHS[conn_l] += Sendbuf[3*i+2];
			}
		}	
	}
	
/*	if(myrank == 0){
		cout<<"nbn = "<<mesh->nbn<<endl;
		for(int i=0;i<mesh->nbn;i++)
			cout<<Sendbuf[3*i+0]<<"\t"<<Sendbuf[3*i+1]<<"\t"<<Sendbuf[3*i+2]<<endl;
	}
*/

    ///MPI Barrier
    //MPI::COMM_WORLD.Barrier();

return;
}
