//==================================================================================================
// Name        : solver.cpp
// Author      : Raghavan Lakshmanan
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the main functions to solve the fem problem.
//==================================================================================================

#include "solver.h"
#include "postProcessor.h"
#include "mpi.h"
//==================================================================================================
// solverControl
//==================================================================================================
void femSolver::solverControl(inputSettings* argSettings, triMesh* argMesh,int my_rank, int num_procs, int prev, int next)
{
    mesh = argMesh;
    settings = argSettings;

    // Calculate Jacobian for all elements in each processor

    for(int i=0;i< mesh->getNe_pro();i++)      
       calculateJacobian(i);
      
    cout << "> Calculate jacobian completed:"<<"\t"<<my_rank<<endl;

    // Calculate element matrices for all elements in each processor 
       
    for(int i=0;i< mesh->getNe_pro();i++)        
          calculateElementMatrices(i);
    
    cout << "> Calculate element matrices completed:"<<"\t"<<my_rank<<endl;

    // Apply the boundary conditions 
 
    for(int i=0;i<  mesh->getNe_pro();i++) 
          applyBoundaryConditions(i);

    cout << "> Apply boundary conditions completed:"<<"\t"<<my_rank<<endl;

    // Solve the equation system

    explicitSolver();  

    return;
}

//==================================================================================================
// calculateJacobian
//==================================================================================================
void femSolver::calculateJacobian(const int e)
{
  
  int i,j,k;  

  double detJ,inv_detJ;
  double dSdKsiEta[2][3];
  double xy[3][2];
  double J[2][2];
  double J_inv[2][2];

  // Derivatives of shape functions with respect to master element coordinates
  
  for(j=0;j<3;j++)
  {
    dSdKsiEta[0][j] = mesh->getME(0)->getDSdKsi(j);   
    dSdKsiEta[1][j] = mesh->getME(0)->getDSdEta(j);  
  }
  
  // Matrix of node coordinates
  
  for(i=0;i<3;i++)
  {
    xy[i][0] = mesh->getNode(mesh->getElem(e)->getConn(i))->getX();  
    xy[i][1] = mesh->getNode(mesh->getElem(e)->getConn(i))->getY();
  }
 
  // Evaluation of Jacobian matrix - matrix matrix multiplication

  for(i=0;i<2;i++)
  {
    for(j=0;j<2;j++)
    {
      J[i][j] = 0.0;
      for(k=0;k<3;k++)
        J[i][j] = J[i][j] + ( dSdKsiEta[i][k] * xy[k][j]) ;
    }
  }  
   
  // Evaluation of determinant of Jacobian matrix    
   
  detJ = (J[0][0]*J[1][1]) - (J[1][0]*J[0][1]) ;

  // Evaluation of inverse of Jacobian matrix

  inv_detJ=1/detJ;
  J_inv[0][0] = J[1][1]*inv_detJ;
  J_inv[0][1] =-J[0][1]*inv_detJ;
  J_inv[1][0] =-J[1][0]*inv_detJ;
  J_inv[1][1] = J[0][0]*inv_detJ;

   
  // Setting absolute value of Jacobian matrix to element
    
  mesh->getElem(e)->setdetJ(abs(detJ));

  // Setting Jacobian matrix and Inverse of Jacobian matrix to element

  for(i=0;i<2;i++)
  {
    for(j=0;j<2;j++)
    {
        mesh->getElem(e)->setJ(i,j,J[i][j]); 
        mesh->getElem(e)->setJ_inv(i,j,J_inv[i][j]);
    }
  } 

 return;

}

//==================================================================================================
// calculateElementMatrices
//==================================================================================================
void femSolver::calculateElementMatrices(const int e)
{
 
  int i,j,f;
 
  double dKsidx,dKsidy,dEtadx,dEtady,detJ;
  double F[3],K[3][3],M[3][3],M_L[3];
  double sum_all,sum_diag;
  
  //Setting value of k from input file

  double k=settings->getD();

  // Evaluation of element stiffness matrix and element mass matrix

  // Initialisation to zero

  for(i=0;i<3;i++)
  {
    F[i] = 0.0;
    M_L[i]=0.0;
    for(j=0;j<3;j++)
    {
       K[i][j]=0.0;
       M[i][j]=0.0;
    }
  }
  
  // Getting shape function derivatives wrt space for element as scalars

  dKsidx = mesh->getElem(e)->getJ_inv(0,0);
  dKsidy = mesh->getElem(e)->getJ_inv(1,0);
  dEtadx = mesh->getElem(e)->getJ_inv(0,1);
  dEtady = mesh->getElem(e)->getJ_inv(1,1);
  detJ   = mesh->getElem(e)->getdetJ();

  for(i=0;i<3;i++)
  {
    
   for(j=0;j<3;j++)
   {   

    // Implementing Gauss quadrature rule of integration

    for(f=0;f<nGQP;f++)
    {          

     K[i][j] = K[i][j] + ( (  ( ( (mesh->getME(f)->getDSdKsi(i)*dKsidx) + (mesh->getME(f)->getDSdEta(i)*dEtadx) ) * \
                                ( (mesh->getME(f)->getDSdKsi(j)*dKsidx) + (mesh->getME(f)->getDSdEta(j)*dEtadx) ) ) + \
                              ( ( (mesh->getME(f)->getDSdKsi(i)*dKsidy) + (mesh->getME(f)->getDSdEta(i)*dEtady) ) * \
                                ( (mesh->getME(f)->getDSdKsi(j)*dKsidy) + (mesh->getME(f)->getDSdEta(j)*dEtady) ) ) ) * k * detJ * mesh->getME(f)->getWeight() );

     M[i][j] = M[i][j] +( mesh->getME(f)->getS(i) * mesh->getME(f)->getS(j) * detJ * mesh->getME(f)->getWeight() );
      
     }

   }
  
  }

  // Lumping of Mass matrix 
 
  // Sum of all elements and diagonal elements

  sum_all  = 0.0;
  sum_diag = 0.0;
  
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
       sum_all = sum_all + M[i][j];

       // sum of diagonal elements
       if(i==j)
       {
           sum_diag = sum_diag + M[i][j];
        }
    }
  }

  for(i=0;i<3;i++)
     M_L[i] = M[i][i] * sum_all * (1/sum_diag);
       
  // Computing flux for element

  for(i=0;i<3;i++)
  {
    // Implementing gauss quadrature rule of integration
    
    for(f=0;f<nGQP;f++)
    {
      F[i] = F[i] +( mesh->getME(f)->getS(i) * settings->getSource() * detJ * mesh->getME(f)->getWeight() );
     }
  } 
     

  // Setting element stiffness matrix,mass matrix,lumped mass matrix and flux vector to element e

  for(i=0;i<3;i++)
  {
    mesh->getElem(e)->setele_flux(i,F[i]);
    mesh->getElem(e)->setele_lum_mass(i,M_L[i]);
    
    for(j=0;j<3;j++)
    {
      mesh->getElem(e)->setele_mat(i,j,K[i][j]); 
    }
  }

  return;
}

//==================================================================================================
// applyBoundaryConditions
//==================================================================================================
void femSolver::applyBoundaryConditions(const int e)
{

  int i,j,FG,ln1,ln2,BCtype,conn[2];
 
  double T1,T2,L,X0,X1,Y0,Y1,BCvalue,HTC,temp;
  double B[3]={0.0};
    
  for(i=0;i<3;i++)
         mesh->getElem(e)->setB(i,0.0);
   
   
  // Loop over 3 edges of element
    
  for(i=0;i<3;i++)
  {
    
   // Check whether the edge is on boundary
   FG = mesh->getElem(e)->getFG(i);
      
   if(FG!=0)
   {
    // Get Boundary condition type and value of that edge
 
    BCtype  = settings->getBC(FG)->getType();
    BCvalue = settings->getBC(FG)->getValue();
   
    // Get local node values of the edge

    ln1 = edgeNodes[i][0];
    ln2 = edgeNodes[i][1];
 
    // Mapping from local node values of edge numbers to node numbers

    conn[0] = mesh->getElem(e)->getConn(ln1);
    conn[1] = mesh->getElem(e)->getConn(ln2);

    if(BCtype==1)  // Dirichilet conditions
    {
      // Get temperature value at the node

      T1 = mesh->getNode(conn[0])->getT();
      T2 = mesh->getNode(conn[1])->getT();

      // check whether temperature at any node is equal to initial temp.If true (which is mostly for all nodes), set new value to node-(default)
      // If false (which is for corner nodes),that node has been previously set for boundary condition value.Take average of 2 values. 

      B[ln1] = BCvalue;
      B[ln2] = BCvalue;

      if(T1 != settings->getInitT())
          B[ln1] = (BCvalue+T1)/2;  

      if(T2 != settings->getInitT())
          B[ln2] = (BCvalue+T2)/2; 
          
      // Set Boundary vector terms to element

      mesh->getElem(e)->setB(ln1,B[ln1]);
      mesh->getElem(e)->setB(ln2,B[ln2]);
    
      // Set Temperature values to nodes that are in Dirichilet boundary

      mesh->getNode(conn[0])->setT(B[ln1]);
      mesh->getNode(conn[1])->setT(B[ln2]);
      
      // Modifying correspoding rows of element stiffness matrix for nodes whose temperature are known
       
      for(j=0;j<3;j++)                                               
      {
        mesh->getElem(e)->setele_mat(ln1,j,0);
        mesh->getElem(e)->setele_mat(ln2,j,0); 
      }

      mesh->getElem(e)->setele_mat(ln1,ln1,1);
      mesh->getElem(e)->setele_mat(ln2,ln2,1);  

      // Flag the nodes which are on boundary - used in solver part

      mesh->getNode(conn[0])->set_flag(BCtype);
      mesh->getNode(conn[1])->set_flag(BCtype);
     
    }

    else if(BCtype==2) // Neumann conditions
    {

     // Find X and Y coordinates of nodes on boundary

     X0 = mesh->getNode(conn[0])->getX();
     Y0 = mesh->getNode(conn[0])->getY();
     X1 = mesh->getNode(conn[1])->getX();
     Y1 = mesh->getNode(conn[1])->getY();

     // Calculate length of edge

     L = pow( ( (X1-X0)*(X1-X0) + (Y1-Y0) *(Y1-Y0) ) ,0.5); 

     // Get temperature value at the node

     T1 = mesh->getNode(conn[0])->getT();
     T2 = mesh->getNode(conn[1])->getT();

     // check whether it is equal to initial temp.If true,set new value. For false i.e for nodes 
     // whose temp. is computed, flux need not be computed again.

     if(T1 == settings->getInitT())
          B[ln1] = (BCvalue*L)/2;  

     if(T2 == settings->getInitT())
          B[ln2] = (BCvalue*L)/2; 
          
     // Set Boundary vector terms to element

     mesh->getElem(e)->setB(ln1,B[ln1]);
     mesh->getElem(e)->setB(ln2,B[ln2]);

    }

    else if(BCtype==3) // Robin (mixed) conditions
    {

     // Get Heat Transfer Coefficient
     HTC = settings->getBC(FG)->getHTC();

     // Find X and Y coordinates of nodes on boundary

     X0 = mesh->getNode(conn[0])->getX();
     Y0 = mesh->getNode(conn[0])->getY();
     X1 = mesh->getNode(conn[1])->getX();
     Y1 = mesh->getNode(conn[1])->getY();

     // Calculate length of edge

     L = pow( ( (X1-X0)*(X1-X0) + (Y1-Y0) *(Y1-Y0) ) ,0.5); 

     // Get temperature value at the node

     T1 = mesh->getNode(conn[0])->getT();
     T2 = mesh->getNode(conn[1])->getT();

     // check whether it is equal to initial temp.If true,set new value. For false i.e for nodes 
     // whose temp. is computed, flux need not be computed again.

     if(T1 == settings->getInitT())
          B[ln1] = (BCvalue*HTC*L)/2;  

     if(T2 == settings->getInitT())
          B[ln2] = (BCvalue*HTC*L)/2; 
          
     // Set Boundary vector terms to element

     mesh->getElem(e)->setB(ln1,B[ln1]);
     mesh->getElem(e)->setB(ln2,B[ln2]);

    // Make the coresponding changes in element level K matrix 
    // Getting appropriate elements in element stiffness matrix,storing it in temp variable,modifying it and setting it back	
    
     temp = mesh->getElem(e)->getele_mat(ln1,ln1);
     temp = temp - ( (-HTC*L)/3);
     mesh->getElem(e)->setele_mat(ln1,ln1,temp);

     temp = mesh->getElem(e)->getele_mat(ln1,ln2);
     temp = temp - ( (-HTC*L)/6);  
     mesh->getElem(e)->setele_mat(ln1,ln2,temp); 

     temp = mesh->getElem(e)->getele_mat(ln2,ln1);
     temp = temp - ( (-HTC*L)/6);
     mesh->getElem(e)->setele_mat(ln2,ln1,temp);

     temp = mesh->getElem(e)->getele_mat(ln2,ln2);
     temp = temp - ( (-HTC*L)/3);
     mesh->getElem(e)->setele_mat(ln2,ln2,temp);      

    }
   
    else
    {
       cout<<">Warning! Unknown boundary conditions!"<<endl;
    }


   } // end of if(FG!=0 ) loop
 
  } // end of for loop
  
   
  return;

}

//==================================================================================================
// explicitSolver
//==================================================================================================
void femSolver::explicitSolver()
{

  int n;
  int conn[3];

  n = mesh->getNn();

  int flag[n];
  int flag_all[n];

  double RHS_e[3];

  double M[n], M_all[n];
  double RHS[n],RHS_all[n];
  
  double time = 0.0;
  double dt = settings->getDt();

  for(int i=0;i<n;i++)
  {
     flag[i]   = 0;
     flag_all[i] = 0;  
  }
     
  // Get flag information of nodes which are in dirichilet boundary of mesh
 
  for(int i=0;i<mesh->getNn();i++)
     flag[i] = mesh->getNode(i)->get_flag();

  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&flag, &flag_all,n,MPI::INT,MPI::SUM);
 
  postProcessor* postP = new postProcessor;

  // Time loop

  for(int t=0;t<settings->getNIter();t++)
  {

   // Write solution at certain time steps
   
   if(t==0)	postP->postProcessorControl(settings, mesh, t, time);
   if(t%settings->getDwf()==0)	postP->postProcessorControl(settings, mesh, t, time);

   // Initialise node level variables at each time step

   for(int i=0;i<n;i++)
   {
     M[i]   = 0.0;
     RHS[i] = 0.0; 
     M_all[i] = 0.0;
     RHS_all[i]=0.0; 
   }

   // Assembling lumped mass matrix M and RHS

   // Loop through all elements

   for(int e=0;e< mesh->getNe_pro();e++) 
   {
    
     // Initialise element level variables at each time step

     for(int i=0;i<3;i++)
        RHS_e[i] = 0.0;

     // Construct element level RHS
     // Access the connectivity of element e

     for(int i=0;i<3;i++)
        conn[i] = mesh->getElem(e)->getConn(i);   
    
     // K[3][3] * T[3]

     for(int i=0;i<3;i++)
     {
       for(int j=0;j<3;j++)
       {
          RHS_e[i] = RHS_e[i] + ( mesh->getElem(e)->getele_mat(i,j) * mesh->getNode(conn[j])->getT() );
       }
     }  
   
     
     // dt * (F + B - K*T)

     for(int i=0;i<3;i++)
     {
        RHS_e[i] = settings->getDt() * (mesh->getElem(e)->getele_flux(i) + mesh->getElem(e)->getB(i) - RHS_e[i]);
     }

     // M[3][3]*T[3] + dt*(F + B - K*T)

     for(int i=0;i<3;i++)
     {
        RHS_e[i] = RHS_e[i] + (mesh->getElem(e)->getele_lum_mass(i) * mesh->getNode(conn[i])->getT() );
     }

     for(int i=0;i<3;i++)
     { 
       M[conn[i]]   = M[conn[i]] +  mesh->getElem(e)->getele_lum_mass(i) ;
       RHS[conn[i]] = RHS[conn[i]] + RHS_e[i];
     } 

    
   }  // end of element loop 

   // Communicating RHS and M across nodes which are shared by processors

   MPI::COMM_WORLD.Barrier();

   MPI::COMM_WORLD.Allreduce(&M, &M_all,n,MPI::DOUBLE,MPI::SUM);
   MPI::COMM_WORLD.Allreduce(&RHS, &RHS_all,n,MPI::DOUBLE,MPI::SUM);
  
     
   // Setting calculated temperature to the nodes which are not on the dirichlet boundary

   for(int i=0;i<n;i++)
   {
     if(flag_all[i]==0)
     {
      mesh->getNode(i)->setT(RHS_all[i]/M_all[i]);
     }
   }  

   // Increase time by dt

   time += dt;
        
 } // end of time loop

/*  if(MPI::COMM_WORLD.Get_rank()==0)
  {

  for(int i=mesh->getnode_index();i< (mesh->getnode_index() + mesh->getNn_pro() );i++)
         cout<<" > temp at node"<<"\t"<<i<<"\t"<<mesh->getNode(i)->getT()<<endl;
}*/

 return;

}

       

