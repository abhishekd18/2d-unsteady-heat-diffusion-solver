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

    MPI::Status status;

    int nn = mesh->getNn();
    int conn_g, conn_l;
    int nprocb, source, destination;

    // Loop through nodes
    for(int i=0;i<mesh->nbn;i++){
	nprocb = mesh->procb[6*i+1];

	for(int j=0;j<nprocb;j++){
		// Get destination and source
		destination = mesh->procb[6*i+j+2];
		source = destination;

		// Prepare send buffer
		conn_g = mesh->procb[6*i+0];
		conn_l = mesh->getMapConn_G_L()[conn_g];
		Sendbuf[0] = M[conn_l];
		Sendbuf[1] = RHS[conn_l];
		
		// Send/recv		
		MPI::COMM_WORLD.Sendrecv(&Sendbuf,2,MPI::DOUBLE,destination,0,&Recvbuf,2,MPI::DOUBLE,source,0,status);

		// Add the received values to proper location
		M[conn_l] += Recvbuf[0];
		RHS[conn_l] += Recvbuf[1];

	}
    }

    ///MPI Barrier
    MPI::COMM_WORLD.Barrier();

return;
}
