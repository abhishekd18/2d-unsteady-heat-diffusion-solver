//==================================================================================================
// Name        : tri.cpp
// Author      : A. Emre Ongut
// Version     : 1.3
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the routines for triangular mesh manipulation such as reading
//               the mesh info from file or shape functions values for triangular elements.
//==================================================================================================

#include "tri.h"
#include <sys/types.h>
#include <sys/stat.h>

//==================================================================================================
// void triMesh::readMeshFiles()
//==================================================================================================
/* File read procedure :
 * 1- Name of the file to be opened is retrieved from the inputSetting obj.
 * 2- File is opened in appropriate format, this is ascii format for minf and binary format for
 *    binary mesh files.
 * 3- Read operation for minf file is straight forward. Binary files are read as size of a double or
 *    int and stored in readStream. Then swapbytes function is called to swap the bytes for the 
 *    correct endianness.
 * 4- Finally obtained data is deep-copied to the mesh data structure. 
 */
//==================================================================================================
void triMesh::readMeshFiles(inputSettings* settings)
{
    ifstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       readStream;     // temperory var used for strings read from files
    double      dummyDouble;    // temperory var used for double values read from files

    //==============================================================================================
    // READ THE MINF FILE
    // This file should hold the number of elements and nodes.
    //==============================================================================================
    cout << "====== Mesh =====" << endl;
    //dummy = settings->getMinfFile();
    dummy = settings->getWdir();
    std::ostringstream myrank_str; 	// output string stream
    myrank_str << MPI::COMM_WORLD.Get_rank();
    dummy.append("proc_").append(myrank_str.str()).append("/minf");

    file.open(dummy.c_str(), ios::in);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    file >> dummy >> ne;
    file >> dummy >> nn;
    cout << "> Number of mesh elements : " << ne << endl;
    cout << "> Number of nodes : " << nn << endl;
    cout << "> File read complete: minf" << endl;
    file.close();

    //Allocation of memeory for the mesh data structure
    node = new triNode[nn];
    elem = new triElement[ne];
    ME   = new triMasterElement[nGQP];
    ME->setupGaussQuadrature();
    ME->evaluateShapeFunctions();
    cout << "> Mesh data structure is created." << endl;

    //==============================================================================================
    // READ THE MXYZ FILE
    // This file contains the node coordinates
    //==============================================================================================
    //dummy = settings->getMxyzFile();
    dummy = settings->getWdir();
    dummy.append("proc_").append(myrank_str.str()).append("/mxyz");

    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }

    double scaleF = settings->getScale();
    readStream = new char [nsd*sizeof(double)];
    file.seekg (0, ios::beg);
    for(int i=0; i<nn; i++)
    {
        file.read (readStream, nsd*sizeof(double));
        swapBytes(readStream, nsd, sizeof(double));
        node[i].setX(*((double*)readStream)*scaleF);
        node[i].setY(*((double*)readStream+1)*scaleF);
    }
    cout << "> File read complete: " << dummy << endl;
    file.close();

    //==============================================================================================
    // READ THE MIEN FILE
    // This file contains the element connectivity
    //==============================================================================================
    //dummy = settings->getMienFile();
    dummy = settings->getWdir();

    dummy.append("proc_").append(myrank_str.str()).append("/mien");

    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    readStream = new char [nen*sizeof(int)];
    file.seekg (0, ios::beg);
    for(int i=0; i<ne; i++)
    {
        file.read (readStream, nen*sizeof(int));
        swapBytes(readStream, nen, sizeof(int));
        for(int j=0; j<nen; j++)
            elem[i].setConn(j, *((int*)readStream+j)-1);
    }
    cout << "> File read complete: " << dummy << endl;
    file.close();

    //==============================================================================================
    // READ THE MRNG FILE
    // This file contains the boundry information
    //==============================================================================================
    //dummy = settings->getMrngFile();
    dummy = settings->getWdir();
    dummy.append("proc_").append(myrank_str.str()).append("/mrng");
    
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    readStream = new char [nef*sizeof(int)];
    file.seekg (0, ios::beg);
    for(int i=0; i<ne; i++)
    {
        file.read (readStream, nef*sizeof(int));
        swapBytes(readStream, nef, sizeof(int));
        for(int j=0; j<nef; j++)
            elem[i].setFG(j, *((int*)readStream+j));
    }
    cout << "> File read complete: " << dummy << endl;
    file.close();

    //==============================================================================================
    // READ THE PROCB FILE
    // This file contains the element connectivity
    //==============================================================================================
    //dummy = settings->getMienFile();
    dummy = settings->getWdir();
    dummy.append("proc_").append(myrank_str.str()).append("/procb");

    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }

    readStream = new char [6*sizeof(int)];
    file.seekg (0, ios::beg);
    int nprocb, nng;

    MPI::COMM_WORLD.Allreduce(&nn,&nng,1,MPI::INT,MPI::SUM);

    map_conn_g_l = new int[nng]();
    for(int i=0; i<nng; i++)
	map_conn_g_l[i] = -1;

    map_conn_l_g = new int[nn]();

    nbn=0;
    for(int i=0; i<nn; i++)
    {
        file.read (readStream, 6*sizeof(int));

	// Set the global to local and vice versa mappings
	map_conn_g_l[*((int*)readStream)] = i;
	map_conn_l_g[i] = *((int*)readStream);

	nprocb = *((int*)readStream+1);
	if(nprocb>0) nbn++;
    }

    procb = new int[6*nbn];

    MPI::COMM_WORLD.Allreduce(&nbn,&nbnmax,1,MPI::INT,MPI::MAX);

    file.seekg (0, ios::beg);

    int ii=0;
    for(int i=0; i<nn; i++)
    {
        file.read (readStream, 6*sizeof(int));

	nprocb = *((int*)readStream+1);
	*((int*)readStream-1);
        //swapBytes(readStream, nen, sizeof(int));
	if(nprocb>0){
	        for(int j=0; j<6; j++){
        	    procb[6*ii+j] = *((int*)readStream+j);
		    if(j>nprocb+1) procb[6*ii+j] = -1;
		}
		ii++;
	}
    }
    cout << "> File read complete: " << dummy << endl;
    file.close();

  
    //==============================================================================================
    // READ THE INITIAL FILE OR INITIALISE
    // This file contains initial field distribution
    //==============================================================================================
    dummy = settings->getWdir();
    dummy.append("proc_").append(myrank_str.str()).append("/data");

    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);

    int initdata = 1; 
    if (file.is_open()==false){
        cout << "> Initial Distribution file is not present.\n> Initializing temperature field to a constant value: " << settings->getInitT() <<" K"<< endl;
        initdata = 0;
    }else{	
	cout<<"> Setting temperature field from initial distribution file..."<<endl;
	readStream = new char [sizeof(double)];
	file.seekg (0, ios::beg);
	for(int i=0; i<nn; i++){
	        file.read (readStream, sizeof(double));
        	swapBytes(readStream, 1, sizeof(double));
        	node[i].setT(*((double*)readStream));
	}
	cout << "> File read complete: " << dummy << endl;
	file.close();
    }

    if(initdata == 0){
	dummyDouble = settings->getInitT();
	for(int i=0; i<nn; i++)
	        node[i].setT(dummyDouble);
    }

    return;
}

/* File write procedure :
 * Write data file so that it can be used for furthur processing 
 * or as initial distribution file for next simulation
 */
void triMesh::writeDataFile(inputSettings* settings){
    ofstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       writeStream;    // temperory var used for strings write to files
    double      dummyDouble;    // temperory var used for double values write to files

    dummy = settings->getWdir();
    std::ostringstream myrank_str; 	// output string stream
    myrank_str << MPI::COMM_WORLD.Get_rank();
    dummy.append("proc_").append(myrank_str.str()).append("/data");

    file.open(dummy.c_str(), ios::out|ios::binary|ios::ate);

    if (file.is_open()==false){
        cout << "Unable to open file : " << dummy << endl;
	exit(0);
    }

    cout<<"> Writing temperature field distribution file."<<endl;
    writeStream = new char [sizeof(double)];

    for(int i=0; i<nn; i++){
       	*((double*)writeStream) = node[i].getT();
       	swapBytes(writeStream, 1, sizeof(double));
	file.write (writeStream, sizeof(double));
    }
	
    cout << "> File write complete: " << dummy << endl;
    file.close();

return;
}

void triMesh::swapBytes (char *array, int nelem, int elsize)
{
    register int sizet, sizem, i, j;
    char *bytea, *byteb;
    sizet = elsize;
    sizem = sizet - 1;
    bytea = new char [sizet];
    byteb = new char [sizet];
    for (i = 0; i < nelem; i++)
    {
        memcpy((void *)bytea, (void *)(array+i*sizet), sizet);
        for (j = 0; j < sizet; j++) 
            byteb[j] = bytea[sizem - j];
        memcpy((void *)(array+i*sizet), (void *)byteb, sizet);
    }
    free(bytea); 
    free(byteb);

    return;
}


//==================================================================================================
// GAUSS QUADRATURE POINTS AND WEIGHTS ARE SET FOR 7 POINT QUADRATURE FORMULA
//==================================================================================================
void triMasterElement::setupGaussQuadrature()
{
    this[0].point[0] = 0.333333333333333;   
    this[0].point[1] = 0.333333333333333;
    this[0].weight   = 0.225 / 2.0;
    
    this[1].point[0] = 0.059715871789770;   
    this[1].point[1] = 0.470142064105115;
    this[1].weight   = 0.132394152788 / 2.0;
    
    this[2].point[0] = 0.470142064105115;   
    this[2].point[1] = 0.059715871789770;
    this[2].weight   = 0.132394152788 / 2.0;
    
    this[3].point[0] = 0.470142064105115;   
    this[3].point[1] = 0.470142064105115;
    this[3].weight   = 0.132394152788 / 2.0;
    
    this[4].point[0] = 0.101286507323456;   
    this[4].point[1] = 0.797426985353087;
    this[4].weight   = 0.125939180544 / 2.0;
    
    this[5].point[0] = 0.101286507323456;   
    this[5].point[1] = 0.101286507323456;
    this[5].weight   = 0.125939180544 / 2.0;
    
    this[6].point[0] = 0.797426985353087;   
    this[6].point[1] = 0.101286507323456;
    this[6].weight   = 0.125939180544 / 2.0;

    return;
}

//==================================================================================================
// EVALUATES SHAPE FUNCTIONS FOR LINEAR TRIANGULAR ELEMENT
//==================================================================================================
void triMasterElement::evaluateShapeFunctions()
{
    double ksi;
    double eta;
    
    for(int i=0; i<nGQP; i++)
    {
        ksi  = this[i].point[0];
        eta  = this[i].point[1];

        this[i].S[0] = 1.0-ksi-eta;
        this[i].S[1] = ksi;
        this[i].S[2] = eta;
        
        this[i].dSdKsi[0] = -1.0;
        this[i].dSdKsi[1] =  1.0;
        this[i].dSdKsi[2] =  0.0;

        this[i].dSdEta[0] = -1.0;
        this[i].dSdEta[1] =  0.0;
        this[i].dSdEta[2] =  1.0;
    }

    return;
}

