//==================================================================================================
// Name        : 2D_Unsteady_Diffusion.cpp
// Author      : Raghavan Lakshmanan
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : This is an educational code desgined for solving 2-dimensional, transient heat
//               diffusion using finite element method. See README file for further details.
//==================================================================================================

#include "settings.h"
#include "tri.h"
#include "solver.h"
#include "postProcessor.h"
#include "mpi.h"

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

    inputSettings*  settings    = new inputSettings;
    triMesh*        mesh        = new triMesh;
    femSolver*      solver      = new femSolver;
    postProcessor*  postP       = new postProcessor;

    int num_procs,my_rank,prev,next;
    double start,end,total_time,avg_time;

//==================================================================================================
      // 1.Pre-Processing Stage
 
      // 1.1 Reading settings file
         settings->readSettingsFile();         
    
      // MPI Initialisation
         MPI::Init(argc,argv);  

      // Getting total number of processors                       
         num_procs = MPI::COMM_WORLD.Get_size();

      // Finding rank of the processor   
         my_rank =  MPI::COMM_WORLD.Get_rank();

         MPI::Status status;
  
      // Calculating rank of previous and next processors

         prev=(my_rank+num_procs-1)%num_procs;
         next=(my_rank+1)%num_procs;

       
      // 1.2 Parallel reading of mesh data
 
         mesh->readMeshFiles(settings, my_rank, num_procs, prev, next);

//==================================================================================================

 

//==================================================================================================


        
     // 2. Solution Stage  

        start = MPI::Wtime();       //  clock start
    
        solver->solverControl(settings, mesh, my_rank, num_procs, prev, next);
 
        end = MPI::Wtime();        // clock end

//==================================================================================================



//==================================================================================================

    // 3. Post-Processing Stage
 
    //  postP->postProcessorControl(settings, mesh);


//==================================================================================================   
 

   // Calculating average time of solver

      total_time = end - start;

      MPI::COMM_WORLD.Reduce(&total_time, &avg_time, 1, MPI::DOUBLE, MPI::MAX, 0);
      if(my_rank==0)
          printf("\n Average time taken for execution of program is %lf\n", avg_time);
    
  // Cleanup

     delete settings;
     delete mesh;
     delete solver;
     delete postP;

     cout << endl << "Ciao :)" << endl;
  
 // MPI Finalize    
  
    MPI::Finalize();

    return 0;
}

