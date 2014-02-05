//==================================================================================================
// Name        : decomposePar.cpp
// Authors     : Abhishek Y. Deshmukh
// Version     : 1.3
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the main functions to solve the fem problem.
//==================================================================================================

#include "decomposePar.h"

//==================================================================================================
// Decomposition
//==================================================================================================
void decomposePar::decompose(inputSettings* settings, triMesh* meshG)
{
    // Get the number of processors    
    int nprocs = settings->getNprocs();

    // Local mesh structure for each processor 
    meshL = new triMesh[nprocs];

    // Allocate array for local number of elements
    nel = new int[nprocs]();

    // Global number of elements
    int neg = meshG->getNe();

/*    // Calculate the number of elements for local mesh
    int sum_ne = 0;
    for(int i=0;i<nprocs;i++){
	nel[i] = ((neg - 1) / nprocs) + 1;
	sum_ne = sum_ne + nel[i];
	if (sum_ne>neg) nel[i] = MAXI(neg-(sum_ne-nel[i]), 0);
    }*/
    for(int i=0;i<nprocs;i++)
	nel[i] = meshG->getNel()[i];

    // Renumbering of elements
    int* newElemN = new int[neg];

    for(int i=0;i<neg;i++)
	newElemN[meshG->getMprm()[i]] = i;

    meshG->setNewElemN(newElemN,neg);

    int off = 0;
    // Distribute and Map global mesh to local mesh
    for(int i=0;i<nprocs;i++){
	meshL[i].setMesh(i,off,meshG,settings);
	off = off + nel[i];
	//cout<<off<<endl;
    }

    // Generate Processor boundary file
    generateProcB(nprocs, settings);
/*    for(int ni=0;ni<nprocs;ni++){
	for(int i=0;i<meshL[ni].getNn();i++){
		cout<<meshL[ni].getNode(i)->getProcB()[0]<<"\t"<<meshL[ni].getNode(i)->getProcB()[1]<<"\t";
		for(int k=0;k<meshL[ni].getNode(i)->getProcB()[1];k++)			
			cout<<meshL[ni].getNode(i)->getProcB()[k+2]<<"\t";
		cout<<endl;
	}cout<<endl;
    }
*/
    delete[] newElemN;
return;
}

//==================================================================================================
// Generate "procb" which contains node level information for the nodes which are on boundaries 
// proc is the processor for which the file will be written
//==================================================================================================
void decomposePar::generateProcB(const int nprocs, inputSettings* settings){

int nn1, nn2, nprocb, GNn;

for(int ni=0;ni<nprocs;ni++){
    for(int nj=0;nj<nprocs;nj++){
	if(ni!=nj){
		nn1 = meshL[ni].getNn();
		nn2 = meshL[nj].getNn();
		nprocb=0;
		// Count the number of nodes on boundary of mesh_1 and mesh_2
		for(int i=0;i<nn1;i++){
			nprocb = meshL[ni].getNode(i)->getProcB()[1];
			for(int j=0;j<nn2;j++){
				GNn = meshL[ni].getMapConn_L_G()[i];
				meshL[ni].getNode(i)->setProcBi(0,GNn);
				if(GNn==meshL[nj].getMapConn_L_G()[j]){
					nprocb++;
					meshL[ni].getNode(i)->setProcBi(1,nprocb);
					meshL[ni].getNode(i)->setProcBi(nprocb+1,nj);
				}
			}

		}
	}
    }

    ofstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       writeStream;    // temperory var used for strings read from files
    double      dummyDouble;    // temperory var used for double values read from files
    std::ostringstream ostr; 	// output string stream
    string	dir;		// directory

    //==============================================================================================
    // WRITE THE PROCB FILE
    // This file contains the element connectivity
    //==============================================================================================
    cout << "====== Writing procb file for processor_"<<ni<< endl;
    ostr << ni; 	     //use the string stream just like cout,except the stream prints not to stdout but to a string.
    dir = settings->getWdir();
    dir = dir.append("proc_").append(ostr.str());

    dummy = dir.append("/procb");

    file.open(dummy.c_str(), ios::out|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }

    int MAXPROCB = 6;
    writeStream = new char [MAXPROCB*sizeof(int)];
    for(int i=0; i<nn1; i++)
    {
       	for(int j=0; j<(MAXPROCB); j++)
       	    *((int*)writeStream+j) = meshL[ni].getNode(i)->getProcB()[j];
       	//swapBytes(writeStream, (nprocb+2), sizeof(int));
    	file.write (writeStream, (MAXPROCB)*sizeof(int));
    }

/*    writeStream = new char [(nprocb+2)*sizeof(int)];
    for(int i=0; i<nn1; i++)
    {
        for(int j=0; j<(nprocb+2); j++)
            *((int*)writeStream+j) = meshL[ni].getNode(i)->getProcB()[j];
        //swapBytes(writeStream, (nprocb+2), sizeof(int));
    	file.write (writeStream, (nprocb+2)*sizeof(int));
    }*/
    cout << "> File write complete: " << dummy << endl;
    file.close();
    delete [] writeStream;



}

return;
}

