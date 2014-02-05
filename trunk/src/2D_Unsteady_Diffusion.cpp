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

    inputSettings*  settings    = new inputSettings;
    triMesh*        mesh        = new triMesh;
    femSolver*      solver      = new femSolver;
    postProcessor*  postP       = new postProcessor;

    /// Pre-Processing Stage
    settings->readSettingsFile();
    mesh->readMeshFiles(settings);

    clock_t start, end;
    start = clock();

    /// Solution Stage
    solver->solverControl(settings, mesh);

    end = clock();
    cout<<"\nTime for solver = "<<((end-start)/(double)CLOCKS_PER_SEC);

    /// Write a data file with field distribution for initial condition
    if(settings->getRestart() == "yes")
	mesh->writeDataFile(settings);

    // Post-Processing Stage integrated with solver
    // postP->postProcessorControl(settings, mesh);

    /// Cleanup
    delete settings;
    delete mesh;
    delete solver;
    delete postP;

    cout << endl << "Ciao :)" << endl;
    return 0;
}

