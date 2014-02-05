//==================================================================================================
// Name        : mpi_comm.h
// Author      : A. Emre Ongut
// Version     : 1.3
// Copyright   : See the copyright notice in the README file.
// Description : Class definitions for triangular elements. There are different classes for storing
//               variables at node and element level.
//==================================================================================================

#ifndef MPI_COMM_H_
#define MPI_COMM_H_

#include "mpi.h"
#include "tri.h"

/*!
 * \brief This class defines the mpi communication between processors.
 * 
 * This class contains the definition of datatypes used for send/recv buffers and the methods
 * to operate on them.
 */
class mpiComm
{
    private:
        /// PRIVATE VARIABLES
    	double Sendbuf[2];
    	double Recvbuf[2]; 
        
    protected:

    public:
        /// DEFAULT CONSTRUCTOR
        mpiComm(){};

        /// DESTRUCTOR
        ~mpiComm(){};

        /// GETTERS  

        /// PUBLIC INTERFACE METHODS
        void communicate(triMesh*, double*, double*);
};

#endif
