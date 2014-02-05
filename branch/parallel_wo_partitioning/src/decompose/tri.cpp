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
    dummy = settings->getMinfFile();
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
    dummy = settings->getMxyzFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    readStream = new char [nsd*sizeof(double)];
    file.seekg (0, ios::beg);
    for(int i=0; i<nn; i++)
    {
        file.read (readStream, nsd*sizeof(double));
        swapBytes(readStream, nsd, sizeof(double));
        node[i].setX(*((double*)readStream));
        node[i].setY(*((double*)readStream+1));
    }
    cout << "> File read complete: " << dummy << endl;
    file.close();

    //==============================================================================================
    // READ THE MIEN FILE
    // This file contains the element connectivity
    //==============================================================================================
    dummy = settings->getMienFile();
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
    dummy = settings->getMrngFile();
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
    // READ THE MIEN FILE
    // This file contains the element connectivity
    //==============================================================================================
    dummy = settings->getMienFile();
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
    // READ THE MPRM FILE
    // This file contains the element connectivity
    //==============================================================================================
    dummy = settings->getMienFile();

    int np = settings->getNprocs();
    std::ostringstream ostr; 	// output string stream
    ostr.fill( '0' );    
    ostr.width( 5 );
    ostr << np;
    dummy = dummy.substr(0,dummy.size()-4).append("mprm.").append(ostr.str());

    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }

    mprm = new int[ne]();
    newElemN = new int[ne];

    nel = new int[np];
    readStream = new char [sizeof(int)];
    file.seekg (0, ios::beg);
    for(int i=0; i<ne+np; i++)
    {
        file.read (readStream, sizeof(int));
        swapBytes(readStream, 1, sizeof(int));
        if(i<ne) mprm[i] = *((int*)readStream) - 1;
	else	nel[i-ne] = *((int*)readStream);
    }

    cout << "> File read complete: " << dummy << endl;
    file.close();

    //==============================================================================================
    // READ THE INITIAL FILE OR INITIALISE IF PRESENT
    // This file contains initial field distribution
    //==============================================================================================
    dummy = settings->getDataFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    int initdata = 1; 

    if (file.is_open()==false){
        cout << "> Initial Distribution file is not present."<< endl;
        initdata = 0;
    }else{	cout<<"> Reading temperature field from initial distribution file..."<<endl;
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



void triMesh::setMesh(const int npe, const int off, triMesh* meshG, inputSettings* settings)
{
    // Global numbers of elements and nodes
    int nng = meshG->getNn();
    int neg = meshG->getNe();

    // Set local number of elements;
    ne = meshG->getNel()[npe];

    elem = new triElement[ne]();
    int prev = -1;
    int conn;

    int k = 0;
    conn_g = new int[ne*nen]();

    // Sort connectivity array
    for(int i=off;i<off+ne;i++){
       	for(int j=0; j<nen; j++){
	    conn_g[k] = meshG->getElem(i)->getConn(j);
	    k++;
	}
    }  

    //bubble sort [?]
    int temp;
    int swap;
    do{
	swap = 0;
    	for(int i=0;i<ne*nen-1;i++){
		if(conn_g[i+1] < conn_g[i]){
			//swap them
			temp = conn_g[i];
			conn_g[i] = conn_g[i+1];
			conn_g[i+1] = temp;
			swap = 1;
		}
    	}
    }while(swap!=0);
	
    // Count unique entries
    int nnl=0;
    for(int i=0;i<ne*nen;i++){
	conn = conn_g[i];
	if(prev!=conn){
		nnl++;
	}
	prev = conn;
    }

    // Set local number of nodes
    nn = nnl;
    //cout<<nnl<<endl;
   
    // Mapping from local to global connectivity
    map_conn_l_g = new int[nnl]();

    // Mapping from global to local connectivity
    map_conn_g_l = new int[nng]();

    prev = -1;

    k = 0;
 
    // Assign both mappings
    for(int i=0;i<ne*nen;i++){
	conn = conn_g[i];
	if(prev!=conn){
		map_conn_l_g[k] = conn;
		map_conn_g_l[conn] = k;
		k++;
	}
	prev = conn;
    }

    int i_l=0, mprm_l;

    // Set local mesh element connectivity using global to local mapping
    for(int i=off;i<off+ne;i++){
	mprm_l = i;
       	for(int j=0; j<nen; j++){
       	    elem[i_l].setConn(j, map_conn_g_l[meshG->getElem(mprm_l)->getConn(j)]);
	}
       	for(int j=0; j<nef; j++)
      	    elem[i_l].setFG(j, meshG->getElem(mprm_l)->getFG(j));
  	i_l++;
    }

    //Allocation of memeory for the node data structure
    node = new triNode[nnl];
    prev = -1;
    int conn_l=0;
    for(int i=0; i<ne*nen; i++){
	    conn = conn_g[i];

    	    // Set X, Y coordinates for unique nodes
	    if(prev!=conn){
	    	conn_l = map_conn_g_l[conn];
		node[conn_l].setX(meshG->getNode(conn)->getX());
		node[conn_l].setY(meshG->getNode(conn)->getY());
		node[conn_l].setT(meshG->getNode(conn)->getT());
	    }
	    prev = conn;
    }

    ofstream    file;           // file name obj
    string      dummy;          // dummy string to hold names
    char*       writeStream;    // temperory var used for strings read from files
    double      dummyDouble;    // temperory var used for double values read from files
    std::ostringstream ostr; 	// output string stream
    string	dir;		// directory

    //==============================================================================================
    // CREATE SEPARATE DIRECTORY FOR EACH PROCESSOR
    // This directory should hold the corresponding files. 
    //==============================================================================================
    cout << "====== Creating processor_"<<npe<<" directory=====" << endl;
    ostr << npe; 	     //use the string stream just like cout,except the stream prints not to stdout but to a string.
    dir = settings->getWdir();
    dir = dir.append("proc_").append(ostr.str());
    //dir = dummy.substr(0,dummy.size()-4).append("proc_").append(ostr.str());
    mkdir(dir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    //mkdir(dummy.su)

    //==============================================================================================
    // WRITE THE MINF FILE
    // This file should hold the number of elements and nodes.
    //==============================================================================================
    cout << "====== Mesh_"<<npe<<" =====" << endl;
    //dummy = settings->getMinfFile();
    dummy = dir.append("/minf");
    //dummy.append(ostr.str());

    file.open(dummy.c_str(), ios::out);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    cout << "> Number of mesh elements : " << ne << endl;
    file << "ne " << ne << endl;
    cout << "> Number of nodes : " << nnl << endl;
    file << "nn " << nnl << endl;
    cout << "> File write complete: "<< dummy << endl;
    file.close();

    //==============================================================================================
    // WRITE THE MXYZ FILE
    // This file contains the node coordinates
    //==============================================================================================
    //dummy = settings->getMxyzFile();
    //dummy.append(ostr.str());
    dir = dir.substr(0,dir.size()-5);
    dummy = dir.append("/mxyz");

    file.open(dummy.c_str(), ios::out|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    writeStream = new char [nsd*sizeof(double)];
    for(int i=0; i<nnl; i++)
    {
        *((double*)writeStream) = node[i].getX();
        *((double*)writeStream+1) = node[i].getY();
        swapBytes(writeStream, nsd, sizeof(double));
    	file.write (writeStream, nsd*sizeof(double));
    }

    cout << "> File write complete: " << dummy << endl;
    file.close();
    delete [] writeStream;

    //==============================================================================================
    // WRITE THE MIEN FILE
    // This file contains the element connectivity
    //==============================================================================================
    //dummy = settings->getMienFile();
    //dummy.append(ostr.str());
    dir = dir.substr(0,dir.size()-5);
    dummy = dir.append("/mien");

    file.open(dummy.c_str(), ios::out|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    writeStream = new char [nen*sizeof(int)];
    for(int i=0; i<ne; i++)
    {
        for(int j=0; j<nen; j++)
            *((int*)writeStream+j) = elem[i].getConn(j)+1;
        swapBytes(writeStream, nen, sizeof(int));
    	file.write (writeStream, nen*sizeof(int));
    }
    cout << "> File write complete: " << dummy << endl;
    file.close();
    delete [] writeStream;

    //==============================================================================================
    // WRITE THE MRNG FILE
    // This file contains the boundry information
    //==============================================================================================
    //dummy = settings->getMrngFile();
    //dummy.append(ostr.str());
    dir = dir.substr(0,dir.size()-5);
    dummy = dir.append("/mrng");

    file.open(dummy.c_str(), ios::out|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    writeStream = new char [nef*sizeof(int)];

    for(int i=0; i<ne; i++)
    {
        for(int j=0; j<nef; j++)
            *((int*)writeStream+j) = elem[i].getFG(j);
        swapBytes(writeStream, nef, sizeof(int));
        file.write (writeStream, nef*sizeof(int));
    }
    cout << "> File write complete: " << dummy << endl;
    file.close();
    delete [] writeStream;
 
    //==============================================================================================
    // WRITE THE INITIAL FILE OR INITIALISE IF PRESENT
    // This file contains initial field distribution
    //==============================================================================================
    dir = dir.substr(0,dir.size()-5);
    dummy = dir.append("/data");
 
    ifstream dfile;
    dfile.open(settings->getDataFile().c_str(), ios::in|ios::binary|ios::ate);

    if (dfile.is_open()==false){
        cout << "> Initial Distribution file is not present.\n"<< endl;
    }else{
	
	file.open(dummy.c_str(), ios::out|ios::binary|ios::ate);
	cout<<"> Writing initial distribution file per processor..."<<endl;
	writeStream = new char [sizeof(double)];

	for(int i=0; i<nn; i++){
	       	*((double*)writeStream) = node[i].getT();
	       	swapBytes(writeStream, 1, sizeof(double));
		file.write (writeStream, sizeof(double));
	}
	
	cout << "> File write complete: " << dummy << endl;
	file.close();
    }

    return;
}

