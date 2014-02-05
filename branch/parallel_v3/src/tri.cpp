//==================================================================================================
// Name        : tri.cpp
// Author      : Raghavan Lakshmanan	
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the routines for triangular mesh manipulation such as reading
//               the mesh info from file or shape functions values for triangular elements.
//==================================================================================================

#include "tri.h"
#include "mpi.h"
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
void triMesh::readMeshFiles(inputSettings* settings,int my_rank, int num_procs, int prev, int next)
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

    
    //==============================================================================================
    // READ THE MPRM FILE
    // This file contains element permutation data for mesh partitioning.
    //==============================================================================================

    int elementperm[ne+num_procs];                          
    
    std::ostringstream nprocs;
    nprocs.fill( '0' );
    nprocs.width( 5 );
    nprocs << num_procs;
    dummy = settings->getMprmFile();
    dummy.append(".").append(nprocs.str());

    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
   
    file.seekg (0, ios::beg);
    file.read ((char*)elementperm, (ne+num_procs)*sizeof(int)); 
    swapBytes((char*)elementperm, (ne+num_procs), sizeof(int));  
  
    cout << "> File read complete: " << dummy << endl;
    file.close();
  
    
    //==============================================================================================
    // READ THE NPRM FILE
    // This file contains node permutation data for mesh partitioning.
    //==============================================================================================

    int nodeperm[nn+num_procs];                         
   
    dummy = settings->getNprmFile();
    dummy.append(".").append(nprocs.str());
    
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
   
    file.seekg (0, ios::beg);
    file.read ((char*)nodeperm, (nn+num_procs)*sizeof(int));  
    swapBytes((char*)nodeperm, (nn+num_procs), sizeof(int));  
  
    cout << "> File read complete: " << dummy << endl;
    file.close();
   
    //Allocation of memory for the mesh data structure

    ne_pro = elementperm[ne+my_rank];      // no. of elements per processor = elementperm[ne+my_rank]
    nn_pro = nodeperm[nn+my_rank];         // no. of nodes per processor    = nodeperm[nn+my_rank]
    node = new triNode[nn];
    elem = new triElement[ne_pro];
  
    ME   = new triMasterElement[nGQP];
    ME->setupGaussQuadrature();
    ME->evaluateShapeFunctions();

    cout << "> Mesh data structure is created." << endl;

    //==============================================================================================
    // READ THE MXYZ FILE 
    // Each processor reads parallely the coordinates of nodes it contains 
    // This file contains the node coordinates
    //==============================================================================================
  
  
    int off_xyz, buff_xyz,dest_rank,dest_pos;
    int prev_nodes[num_procs]; 
    double *xyz;
    double *xyz_p;
   
    dummy = settings->getMxyzFile();
    MPI::File file_xyz = MPI::File::Open(MPI::COMM_WORLD,dummy.c_str(),MPI::MODE_RDONLY,MPI::INFO_NULL);
    buff_xyz = nodeperm[nn+my_rank]*nsd;                       // buffer size of node coordinates = no.of nodes in each processor * no.of node dimensions    
    xyz = (double *)malloc(buff_xyz*sizeof(double));

    // Calculate previous nodes in each processor

    prev_nodes[0]=0;
    for(int i=1;i<=num_procs;i++)
          prev_nodes[i] = prev_nodes[i-1] + nodeperm[nn+i-1];
   
    // Calcuate beginning node index and offset for reading file in each processor 
 
    if(my_rank==0)
    {
      node_index = 0;
      off_xyz =0;
    }   
    else
    {
      node_index = prev_nodes[my_rank];
      off_xyz = prev_nodes[my_rank]*nsd*sizeof(double);
    }  
    
    // Reading the node coordinates data 

    file_xyz.Set_view( off_xyz, MPI::DOUBLE, MPI::DOUBLE, "native" , MPI::INFO_NULL);    
    file_xyz.Read( xyz, buff_xyz, MPI::DOUBLE );
    swapBytes((char*)xyz, buff_xyz, sizeof(double)); 

    // Node permutation using nprm file

    // ( A dummy array of size equal to buffer size of xyz in each processor is created to hold permuted node coordinates data

    xyz_p = (double *)malloc(buff_xyz*sizeof(double));
 
    // Creating window of dummy array in each processor where  permuted node coordinates to be received   

    MPI::Win window_1 = MPI::Win::Create(xyz_p,buff_xyz*8,8,MPI::INFO_NULL,MPI::COMM_WORLD);
    window_1.Fence(0);
  
    int j=0;
    int nodeperm_value;
   
    // Sending coordinates of every node to destination process in destination position calculated based on node permutation data

    for(int i=0; i< nn_pro; i++)
    {
    
       nodeperm_value = nodeperm[node_index + i]-1;  // Calculating permuted node index of every node

       // Calculating destination rank and position of permuted node index

       for(int k=0;k<num_procs;k++)
       {
         if(nodeperm_value>=prev_nodes[k] && nodeperm_value<prev_nodes[k+1])
         {
            dest_rank = k;
            dest_pos = nodeperm_value - prev_nodes[k];
         }
       }
    
       // put node coordinates data from each processor to destination rank in destination position

       window_1.Put( &xyz[i+j], 2, MPI::DOUBLE, dest_rank, 2*dest_pos, 2, MPI::DOUBLE);
      
       j = j+1;
    
     } 

    window_1.Fence(0);
    window_1.Free();

    cout << "> File read complete: " << dummy << endl;

    file_xyz.Close();
    delete [] xyz;
 
    //==============================================================================================
    // READ THE MIEN FILE
    // Each processor reads parallely the element connectivity data of elements it contains 
    // This file contains the element connectivity
    //==============================================================================================

    int off_conn,buff_conn,conn0,conn1,conn2;
    int prev_elements[num_procs]; 
    int *ele_conn,*ele_conn_p;
   
    dummy = settings->getMienFile();
    MPI::File file_conn = MPI::File::Open(MPI::COMM_WORLD,dummy.c_str(),MPI::MODE_RDONLY,MPI::INFO_NULL);
    buff_conn = elementperm[ne+my_rank]*nen;  // buffer size of element connectivity = no.of elements in each processor * no.of element nodes 
    ele_conn = (int *)malloc(buff_conn*sizeof(int));
 
    // Calculate previous elements in each processor

    prev_elements[0]=0;
    for(int i=1;i<=num_procs;i++)
          prev_elements[i] = prev_elements[i-1] + elementperm[ne+i-1];
   
    // Calcuate beginning element index and offset for reading file in each processor   

    if(my_rank==0)
    {
      element_index = 0;
      off_conn = 0;
    } 
    else
    { 
      element_index = prev_elements[my_rank];
      off_conn = prev_elements[my_rank] * nen * sizeof(int);
    }
 
    // Reading element connectivity data

    file_conn.Set_view( off_conn, MPI::INT, MPI::INT, "native" , MPI::INFO_NULL);    
    file_conn.Read( ele_conn, buff_conn, MPI::INT );
    swapBytes((char*)ele_conn, buff_conn, sizeof(int)); 

    // Element permutation using mprm file
    // ( A dummy array of size equal to buffer size of ele_conn in each processor is created to hold permuted element connectivity data

    ele_conn_p = (int *)malloc(buff_conn*sizeof(int));

    // Creating window of dummy array in each processor where  permuted element connectivity data to be received 

    MPI::Win window_2 = MPI::Win::Create(ele_conn_p,buff_conn*4,4,MPI::INFO_NULL,MPI::COMM_WORLD);
    window_2.Fence(0);

    // Sending element connectivity data of every element to destination process in destination position calculated based on element permutation data
  
    j=0;
    int elementperm_value;
   
    for(int i=0; i < ne_pro; i++)
    {
    
       elementperm_value = elementperm[element_index + i]-1;    // Calculating permuted element index of every element

       // Calculating destination rank and position of permuted element index

       for(int k=0;k<num_procs;k++)
       {
           if(elementperm_value>=prev_elements[k] && elementperm_value<prev_elements[k+1])
           {
               dest_rank = k;
               dest_pos = elementperm_value - prev_elements[k];
           }
       }

       // Put element connectivity data from each processor to destination rank in destination position
    
       window_2.Put( &ele_conn[i+j], 3, MPI::INT, dest_rank, 3*dest_pos, 3, MPI::INT);
      
       j = j+2;
    
    } 

    window_2.Fence(0);
    window_2.Free();
  
    // Setting permuted element connectivity data to each element in processor

    j=0;
    for(int i=0; i< ne_pro; i++)         
    {

     // Checking whether all nodes corresponding to element connectivity data of all elements in each processor is matching with node information in each processor
     // If not set flag as 1. ie. Setting that particular node coordinates data is not present in its own processor 

      conn0 = nodeperm[ele_conn_p[i+j]-1]-1;                                                
          
      if( conn0<node_index || conn0>node_index + (buff_xyz/nsd) )
           node[conn0].set_common(1);
                 
      conn1 = nodeperm[ele_conn_p[i+j+1]-1]-1;

      if( conn1<node_index || conn1>node_index + (buff_xyz/nsd) )
            node[conn1].set_common(1);
         
      conn2 = nodeperm[ele_conn_p[i+j+2]-1]-1;

      if( conn2<node_index || conn2>node_index + (buff_xyz/nsd) )
            node[conn2].set_common(1);
                  
      elem[ i ].setConn(0,conn0);
      elem[ i ].setConn(1,conn1);
      elem[ i ].setConn(2,conn2);

      j= j+2;
     
     }
  
    cout << "> File read complete: " << dummy << endl;
  
    file_conn.Close();
  
    delete [] ele_conn;
    delete [] ele_conn_p; 
   
    //==============================================================================================
    // READ THE MRNG FILE
    // Each processor reads parallely the boundary data of elements it contains
    // This file contains the boundry information
    //==============================================================================================

    int off_mrng, buff_mrng,mrng_index; 
    int *ele_mrng,*ele_mrng_p;

    dummy = settings->getMrngFile();   
    MPI::File file_mrng = MPI::File::Open(MPI::COMM_WORLD,dummy.c_str(),MPI::MODE_RDONLY,MPI::INFO_NULL); 
    buff_mrng = elementperm[ne+my_rank]*nen;         // buffer size of element boundary = no.of elements in each processor * no.of element edges 
    ele_mrng = (int *)malloc(buff_mrng*sizeof(int));
  
    // Calcuate beginning element index and offset for reading file in each processor  
 
    if(my_rank==0)
    {
      mrng_index = 0;
      off_mrng = 0;
    } 
    else
    { 
      mrng_index = prev_elements[my_rank];
      off_mrng = prev_elements[my_rank] * nen * sizeof(int);
    }
   
    // Reading element boundary data

    file_mrng.Set_view( off_mrng, MPI::INT, MPI::INT, "native" , MPI::INFO_NULL);    
    file_mrng.Read( ele_mrng, buff_mrng, MPI::INT );
    swapBytes((char*)ele_mrng, buff_mrng, sizeof(int)); 
  
    // Element boundary permutation using mprm file
    // ( A dummy array of size equal to buffer size of ele_mrng in each processor is created to hold permuted element boundary data

    ele_mrng_p = (int *)malloc(buff_mrng*sizeof(int));

    //  Creating window of dummy array in each processor where  permuted element boundary data to be received  

    MPI::Win window_3 = MPI::Win::Create(ele_mrng_p,buff_mrng*4,4,MPI::INFO_NULL,MPI::COMM_WORLD);
    window_3.Fence(0);  
 
   j=0;
   
    for(int i=0; i< ne_pro; i++)
    {
    
       elementperm_value = elementperm[element_index + i]-1;   // Calculating permuted element index of every element

       // Calculating destination rank and position of permuted element index
     
       for(int k=0;k<num_procs;k++)
       {
         if(elementperm_value>=prev_elements[k] && elementperm_value<prev_elements[k+1])
         {
           dest_rank = k;
           dest_pos = elementperm_value - prev_elements[k];
         }
       } 

       // Put element boundary data from each processor to destination rank in destination position
  
       window_3.Put( &ele_mrng[i+j], 3, MPI::INT, dest_rank, 3*dest_pos, 3, MPI::INT);
      
       j = j+2;
    
    } 

    window_3.Fence(0);
    window_3.Free();

    // Setting permuted element boundary data to each element in processor

    j=0;
    for(int i=0; i< ne_pro; i++)         
    {
        elem[  i ].setFG(0, ele_mrng_p[i+j]);
        elem[  i ].setFG(1, ele_mrng_p[i+j+1]);
        elem[  i ].setFG(2, ele_mrng_p[i+j+2]);
  
        j= j+2;
    }
  
    cout << "> File read complete: " << dummy << endl;
    file_mrng.Close(); 
   
    delete [] ele_mrng;
    delete [] ele_mrng_p; 
  

    //================================================================================================================
    // INTERMEDIATE STORAGE OF NODE INFORMATION FOR ELEMENTS THAT SHARE NODE WITH OTHER PROCESSORS
    // Node coordinates data of node that is missing in its own processor is received from other processors where the data is
    //================================================================================================================


    // Find total number of nodes whose information is missing in the current processor

    int sum;
    
    sum = 0;
    for(int i=0;i<nn;i++)
          sum = sum + node[i].get_common(); 
       
    
    // Setting permuted node coordinates data to each node in processor

    j=0;
    for(int i=0; i< nn_pro; i++)         
    {
        node[node_index + i ].setX(xyz_p[i+j]);
        node[ node_index + i ].setY(xyz_p[i+j+1]);
  
        j= j+1;
    } 

    int src_rank,src_pos; 
    double xyz_n[2];

    xyz_n[0] = 0.0;
    xyz_n[1] = 0.0; 

    //  Creating window of node coordinates array in each processor from which node data of certain nodes will be get 

    MPI::Win window_4 = MPI::Win::Create(xyz_p,buff_xyz*8,8,MPI::INFO_NULL,MPI::COMM_WORLD);
    window_4.Fence(0);
   
    // Getting node coordinates data by nodes whose flag is marked as 1. ie.whose node data is missing in its own processor

    for(int i=0;i<nn;i++)
    {

      window_4.Fence(0);
 
      // Proceed only for nodes whose flag is marked as 1   
      
      if( node[i].get_common()==1);
      {
         
         // Calcuate source rank and destination of missing nodes   
                   
         for(int k=0;k<num_procs;k++)
         {
            if( i>=prev_nodes[k] && i<prev_nodes[k+1])
            {
               src_rank = k;
               src_pos = i - prev_nodes[k];
            }
         } 

         // Get node connectivity data to dummy array

         window_4.Fence(0);       
         window_4.Get( &xyz_n,2, MPI::DOUBLE, src_rank, 2*src_pos, 2, MPI::DOUBLE);
         window_4.Fence(0);

         // Set the node connectivity data to that particular node in processor

         node[ i ].setX(xyz_n[0]);
         node[ i ].setY(xyz_n[1]);

       }
    }   
 
    window_4.Fence(0);
    window_4.Free(); 

    delete [] xyz_p;
  
  
 /*   //==============================================================================================
    // READ THE INITIAL FILE OR INITIALISE
    // This file contains initial field distribution
    //==============================================================================================
    dummy = settings->getDataFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    readStream = new char [sizeof(double)];
    file.seekg (0, ios::beg);
    for(int i=0; i<nn; i++)
    {
        file.read (readStream, sizeof(double));
        swapBytes(readStream, 1, sizeof(double));
        node[i].setT(*((double*)readStream));
    }
    cout << "> File read complete: " << dummy << endl;
    file.close();
*/

    // Setting initial temperature value and flag for all nodes in processor

    dummyDouble = settings->getInitT();
    for(int i=0;i<nn;i++)
    {
        node[i].setT(dummyDouble);
        node[i].set_flag(0);
    }
   
    // Creating directory to hold vtk files in post processing

    std::ostringstream ostr; 	// output string stream
    string	dir;		// directory


    cout << "====== Creating processor_"<<my_rank<<" directory=====" << endl;
    ostr << my_rank; 	     //use the string stream just like cout,except the stream prints not to stdout but to a string.
    dir = settings->getWdir();
    dir = dir.append("proc_").append(ostr.str());
    mkdir(dir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    cout << "> tri.cpp read complete: " << endl;


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
// GAUSS QUADRATURE PIONTS AND WEIGHTS ARE SET FOR 7 POINT QUADRATURE FORMULA
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


