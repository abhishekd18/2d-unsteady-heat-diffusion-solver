//==================================================================================================
// Name        : decompose.cpp
// Author      : 
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : This is an educational code desgined for solving 2-dimensional, transient heat
//               diffusion using finite element method. See README file for further details.
//==================================================================================================

#include "settings.h"
#include "tri.h"
#include "decomposePar.h"

using namespace std;

int main(int argc, char **argv)
{
//==================================================================================================
//  Decompose utility
//==================================================================================================

    inputSettings*  settings    = new inputSettings;
    triMesh*        mesh        = new triMesh;
    decomposePar*   dcompPar    = new decomposePar;

    /// Pre-Processing Stage
    settings->readSettingsFile();
    mesh->readMeshFiles(settings);

    /// Decompose into nprocs domains
    //if(settings->getNprocs()>1)
	dcompPar->decompose(settings, mesh);

    /// Cleanup
    delete settings;
    delete mesh;
    delete dcompPar;

    cout << endl << "Partitioning Complete!" << endl;
    return 0;
}

