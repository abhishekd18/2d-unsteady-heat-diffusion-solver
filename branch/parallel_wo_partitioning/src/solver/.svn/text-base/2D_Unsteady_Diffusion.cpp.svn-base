//==================================================================================================
// Name        : 2D_Unsteady_Diffusion.cpp
// Author      : 
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : This is an educational code desgined for solving 2-dimensional, transient heat
//               diffusion using finite element method. See README file for further details.
//==================================================================================================

#include "settings.h"
#include "tri.h"
#include "solver.h"
#include "postProcessor.h"

using namespace std;

int main(int argc, char **argv)
{
//==================================================================================================
//  MAIN PROGRAM FLOW
//  1. Pre-Processing Stage
//      1.1. Settings
//      1.2. Mesh
//  2. Solution Stage
//  3. Post-Processing Stage
//==================================================================================================

    /// MPI Initialize 

    MPI::Init(argc, argv);
    MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
    try {
	int myrank = MPI::COMM_WORLD.Get_rank();
	std::cout << "I am " << myrank << std::endl;
    }
    catch (MPI::Exception e) {
    std::cout << "MPI ERROR: " << e.Get_error_code() << " - " << e.Get_error_string() << std::endl;
    }

    inputSettings*  settings    = new inputSettings;
    triMesh*        mesh        = new triMesh;
    femSolver*      solver      = new femSolver;
    postProcessor*  postP       = new postProcessor;

    /// Pre-Processing Stage
    settings->readSettingsFile();

    mesh->readMeshFiles(settings);

	double start, end;
	start = MPI::Wtime();

    /// Solution Stage
    solver->solverControl(settings, mesh);

	end = MPI::Wtime();

	double local_time = end-start;

	double MAX_time;
	MPI::COMM_WORLD.Reduce(&local_time,&MAX_time,1,MPI::DOUBLE,MPI::MAX,0);

    /// Write a data file with field distribution for initial condition
    if(settings->getRestart() == "yes")
	mesh->writeDataFile(settings);

    if(MPI::COMM_WORLD.Get_rank()==0)
	cout<<"\nTime for solver = "<<MAX_time;

    // Post-Processing Stage integrated with solver
    //if(settings->getNprocs()==1)
    	//postP->postProcessorControl(settings, mesh,1,0.0);

    /// Cleanup
    delete settings;
    //delete mesh;
    delete solver;
    delete postP;

    cout << endl << "Ciao :)" << endl;

    MPI::Finalize();
    return 0;
}

