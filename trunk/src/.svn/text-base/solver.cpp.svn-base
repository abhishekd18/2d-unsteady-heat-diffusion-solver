//==================================================================================================
// Name        : solver.cpp
// Authors     : Abhishek Y. Deshmukh, Raghavan Lakshmanan, Mohsin Ali Chaudry
// Version     : 1.3
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the main functions to solve the fem problem.
//==================================================================================================

#include "solver.h"
#include "postProcessor.h"

//==================================================================================================
// solverControl
//==================================================================================================
void femSolver::solverControl(inputSettings* argSettings, triMesh* argMesh)
{
    mesh = argMesh;
    settings = argSettings;

    ///Calculate Jacobian for all elements
    for(int e=0;e<mesh->getNe();e++)
    	femSolver::calculateJacobian(e);

    ///Calculate element matrices for all elements
    for(int e=0;e<mesh->getNe();e++)
	femSolver::calculateElementMatrices(e);

    ///Apply the boundary conditions 
    for(int e=0;e<mesh->getNe();e++)
	femSolver::applyBoundaryConditions(e);

    ///Solve the equation system 
    femSolver::explicitSolver();

    return;
}

//==================================================================================================
// calculateJacobian
//==================================================================================================
void femSolver::calculateJacobian(const int e)
{
    int* conn = new int[3]();
    double* X = new double[3]();
    double* Y = new double[3]();

    ///Access the coordinates of the nodes of the element 'e'
    for(int i=0;i<3;i++){
	conn[i] = mesh->getElem(e)->getConn(i);
	X[i] = mesh->getNode(conn[i])->getX();
	Y[i] = mesh->getNode(conn[i])->getY();
    }

    ///Access the derivatives of shape functions
    ///dSdXi and dSdEta are constant for the reference element so it is sufficient to access values
    ///for 0th Master element.
    double* dSdXi = new double[3]();
    double* dSdEta = new double[3]();
    for(int i=0;i<3;i++){
	dSdXi[i] = mesh->getME(0)->getDSdKsi(i);
	dSdEta[i] = mesh->getME(0)->getDSdEta(i);
    }

    ///Matrix-Matrix multiplication for calculation of Jacobian Matrix
    ///		     _	        _
    ///		    | J[0]  J[1] |
    ///		J = | 		 |
    ///		    |_J[2]  J[3]_|		
		    
    double* J = new double[4]();
    for(int j=0;j<3;j++){
	J[0] = J[0] + dSdXi[j]*X[j];
	J[1] = J[1] + dSdXi[j]*Y[j];
	J[2] = J[2] + dSdEta[j]*X[j];
	J[3] = J[3] + dSdEta[j]*Y[j];
    }
    
    ///Set the Jacobian matrix for each element
    mesh->getElem(e)->setJ(J);

    ///Calculate determinant of J
    double det_J = J[0]*J[3] - J[1]*J[2];
    mesh->getElem(e)->setDetJ(fabs(det_J));

    ///Invert Jacobian Matrix by Adjoint method
    double* J_inv = new double[4]();
    J_inv[0] = J[3]/det_J;
    J_inv[1] = -J[1]/det_J;
    J_inv[2] = -J[2]/det_J;
    J_inv[3] = J[0]/det_J;

    ///Set the Inverse of Jacobian matrix for each element
    mesh->getElem(e)->setJinv(J_inv);

    delete[] conn;
    delete[] X;
    delete[] Y;
    delete[] dSdXi;
    delete[] dSdEta;
    delete[] J;
    delete[] J_inv;
    return;
}

//==================================================================================================
// calculateElementMatrices
//==================================================================================================
void femSolver::calculateElementMatrices(const int e)
{
    double K[9] = {0.0};
    double M[9] = {0.0};
    double dSdX[3] = {0.0};
    double dSdY[3] = {0.0};
    double F[3] = {0.0};
    double S[3] = {0.0};

    ///Access the derivatives of shape functions
    ///dSdXi and dSdEta are constant for the reference element so it is sufficient to access values
    ///for 0th Master element.
    double dSdXi[3] = {0.0};
    double dSdEta[3] = {0.0};

    ///Get Diffusion coefficient
    double k = settings->getD();

    ///Access the coordinates of the nodes of the element 'e'
    double X[3], Y[3];
    int conn[3];
    for(int i=0;i<3;i++){
	conn[i] = mesh->getElem(e)->getConn(i);
	X[i] = mesh->getNode(conn[i])->getX();
	Y[i] = mesh->getNode(conn[i])->getY();
    }

    ///Get absolute value of determinant of Jacobian
    double det_J_abs = mesh->getElem(e)->getDetJ();

    /// Calculate the right hand side source term
    double source = settings->getSource()/(settings->getRho() * settings->getCp());

    double weight;
    ///Gauss Quadrature rule for integration
    for(int g=0;g<nGQP;g++){

	/// Access weight for current gauss point
	weight = mesh->getME(g)->getWeight();

	///Access Shape function values
    	for(int i=0;i<3;i++)
		S[i] = mesh->getME(g)->getS(i);

	///Access derivatives of shape functions
    	for(int i=0;i<3;i++){
		dSdXi[i] = mesh->getME(g)->getDSdKsi(i);
		dSdEta[i] = mesh->getME(g)->getDSdEta(i);
	}
     
        ///Evaluate dSdX and dSdY at Gauss Quadrature Points. 
        for(int i=0;i<3;i++){
		dSdX[i] = mesh->getElem(e)->getJinv()[0]*dSdXi[i] + mesh->getElem(e)->getJinv()[1]*dSdEta[i];
		dSdY[i] = mesh->getElem(e)->getJinv()[2]*dSdXi[i] + mesh->getElem(e)->getJinv()[3]*dSdEta[i];
        }
 
        ///Construct the element matrix
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			K[3*i+j] = K[3*i+j] + \
				  k*(dSdX[i]*dSdX[j] + dSdY[i]*dSdY[j])*det_J_abs*weight;

		}
	}

	///Construct the mass matrix
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			M[3*i+j] = M[3*i+j] + S[i]*S[j]*det_J_abs*weight;
		}
	}
	
	///Construct the element level Source vector
	for(int i=0;i<3;i++){
		if(Y[i]>=0.0004)
			F[i] = F[i] + S[i]*source*det_J_abs*weight;
	}

    }//GQR loop end

    ///Set the element matrix in element level structure
    mesh->getElem(e)->setK(K);

    ///Lumping of Mass matrix for all elements

    ///Initialize sums
    double sum_all = 0.0;
    double sum_diag = 0.0;

    ///Lumped mass matrix is diagonal matrix so only diagonal is stored
    double M_l[3] = {0.0};

    ///Sum of all elements
    for(int i=0;i<9;i++)
	sum_all = sum_all + M[i];

    //Sum of diagonal elements
    for(int i=0;i<3;i++)
	sum_diag = sum_diag + M[3*i+i];

    ///Lumped local mass matrix (only diagonal elements)
    for(int i=0;i<3;i++)
	M_l[i] = M[3*i+i]*sum_all/sum_diag;

    ///Set the mass matrix in element level structure
    mesh->getElem(e)->setM(M_l);

    ///Set the element matrix in element level structure
    mesh->getElem(e)->setF(F);

    return;
}

//==================================================================================================
// applyBoundaryConditions
//==================================================================================================
void femSolver::applyBoundaryConditions(const int e)
{ 
    int BCType, FG, n1, n2;
    int conn[2];
    double BCValue, HTC, T1, T2, X0, Y0, X1, Y1, L;
    double B[3]={0.0};
    double K[9]={0.0};
    mesh->getElem(e)->setB(B);

    // Material properties
    double rho = settings->getRho();
    double cp = settings->getCp();

    ///Loop over Faces of element 'e'
    for(int i=0;i<3;i++){
	
	for(int j=0;j<9;j++)	K[j] = 0.0;

	FG = mesh->getElem(e)->getFG(i);

	if(FG!=0){

	    BCType = settings->getBC(FG)->getType();
	    BCValue = settings->getBC(FG)->getValue();

	    ///Access the local indices for nodes of the current Face
	    n1 = edgeNodes[i][0];
	    n2 = edgeNodes[i][1];

    	    ///Access the connectivity of nodes of the Face
	    conn[0] = mesh->getElem(e)->getConn(n1);
	    conn[1] = mesh->getElem(e)->getConn(n2);

	    if(BCType==1){	///Dirichlet

   		///Set the temperature values of the nodes of the element on boundaries

		///Get the temperature at the node. Default is Intial value.
		T1 = mesh->getNode(conn[0])->getT();
		T2 = mesh->getNode(conn[1])->getT();

		///Check if Temperature at some node is already set to some value. This can
		///be the case with nodes which lie on multiple boundaries. If this is the 
		///case, take average of the both values
/*		if(T1 != settings->getInitT()){
			B[n1] = (BCValue+T1)/2.0;
			B[n2] = BCValue;
			mesh->getNode(conn[0])->setT((BCValue+T1)/2.0);
			mesh->getNode(conn[1])->setT(BCValue);
		}else if(T2 != settings->getInitT()){
			B[n1] = BCValue;
			B[n2] = (BCValue+T2)/2.0;
			mesh->getNode(conn[0])->setT(BCValue);
			mesh->getNode(conn[1])->setT((BCValue+T2)/2.0);
		}else{*/
			B[n1] = BCValue;
			B[n2] = BCValue;
			mesh->getNode(conn[0])->setT(BCValue);
			mesh->getNode(conn[1])->setT(BCValue);
//		}
		
		///Make the corresponding entries in element level K matrix = 1 and rest equate to zero
		///in the same row
		K[3*n1+n1] = 1.0;
		K[3*n2+n2] = 1.0;

		///Copy row correspoding to unknown temperature to modified K
		for(int m=0;m<3;m++)
			K[3*(3-(n1+n2))+m] = mesh->getElem(e)->getK()[3*(3-(n1+n2))+m];

		///Set B
		mesh->getElem(e)->setB(B);
		///Set K
		mesh->getElem(e)->setK(K);

	    	///Flag the nodes on Dirichlet boundary and assign BC type
	    	mesh->getNode(conn[0])->setBC_type(BCType);
	    	mesh->getNode(conn[1])->setBC_type(BCType);

	    }else if(BCType==2){	///Neumann

		///Get the X and Y coordinates of the nodes on boundary
		X0 = mesh->getNode(conn[0])->getX();
		Y0 = mesh->getNode(conn[0])->getY();
		X1 = mesh->getNode(conn[1])->getX();
		Y1 = mesh->getNode(conn[1])->getY();
		
		///Calculate length of the edge
		L = pow((X0 - X1)*(X0 - X1) + (Y0 - Y1)*(Y0 -Y1), 0.5);

   		///Set the flux values at the nodes of the element on boundaries
		///Get the temperature at the node. Default is Intial value.
		T1 = mesh->getNode(conn[0])->getT();
		T2 = mesh->getNode(conn[1])->getT();

		///Check if Temperature at some node is already set to some value. If this
		///is the case, we don't apply any flux value to that node, because we 
		///already know the value of temperature and we don't solve for it. 
//		if(T1 == settings->getInitT())
			B[n1] = BCValue*L/(2.0*rho*cp);
//		if(T2 == settings->getInitT())
			B[n2] = BCValue*L/(2.0*rho*cp);


		mesh->getElem(e)->setB(B);

	    }else if(BCType==3){	///Robin (Mixed)
		
		///Get Heat Transfer Coefficient
	    	HTC = settings->getBC(FG)->getHTC();

		///Get the X and Y coordinates of the nodes on boundary
		X0 = mesh->getNode(conn[0])->getX();
		Y0 = mesh->getNode(conn[0])->getY();
		X1 = mesh->getNode(conn[1])->getX();
		Y1 = mesh->getNode(conn[1])->getY();
		
		///Calculate length of the edge
		L = pow((X0 - X1)*(X0 - X1) + (Y0 - Y1)*(Y0 -Y1), 0.5);

   		///Set the flux values at the nodes of the element on boundaries
		///Get the temperature at the node. Default is Intial value.
		T1 = mesh->getNode(conn[0])->getT();
		T2 = mesh->getNode(conn[1])->getT();

		///Check if Temperature at some node is already set to some value. If this
		///is the case, we don't apply any flux value to that node, because we 
		///already know the value of temperature and we don't solve for it. 
//		if(T1 == settings->getInitT())
			B[n1] = BCValue*HTC*L/(2.0*rho*cp);
//		if(T2 == settings->getInitT())
			B[n2] = BCValue*HTC*L/(2.0*rho*cp);				
		

		mesh->getElem(e)->setB(B);

		std::memcpy(K,mesh->getElem(e)->getK(),sizeof(K));

		///Make the coresponding changes in element level K matrix 
		K[3*n1+n1] = K[3*n1+n1] - (-HTC)*L/(3.0*rho*cp);
		K[3*n2+n2] = K[3*n2+n2] - (-HTC)*L/(3.0*rho*cp);
		K[3*n1+n2] = K[3*n1+n2] - (-HTC)*L/(6.0*rho*cp);
		K[3*n2+n1] = K[3*n2+n1] - (-HTC)*L/(6.0*rho*cp);


		mesh->getElem(e)->setK(K);
		//cout<<"HTC = "<<HTC<<"Tinf = "<<BCValue<<endl;
	    }else{
		cout<<">Warning! Unknown boundary conditions!"<<endl;
	    }

        }//End of if(FG!=0)

    }///Loop over faces end

    return;
}


//==================================================================================================
// explicitSolver
//==================================================================================================
void femSolver::explicitSolver()
{
    ///Element level variables
    int conn[3];
    double RHS_e[3];
    double M_l[3];

    ///Node level variables
    int nn = mesh->getNn();
    double* M = new double [nn]();
    double* RHS = new double [nn]();

    double time = 0.0;
    double dt = settings->getDt();
    postProcessor*  postP = new postProcessor;

    ///Time loop start	
    for(int t=0;t<=settings->getNIter();t++){

	///Write solution at certain time steps
	if(t%settings->getDwf()==0)	postP->postProcessorControl(settings, mesh, t, time);

   	///Initialize node level variables
	for(int node=0;node<nn;node++){
		M[node] = 0.0;
		RHS[node] = 0.0;
	}

	///Loop through all elements
	for(int e=0;e<mesh->getNe();e++){

		for(int i=0;i<3;i++)
   			RHS_e[i] = 0.0;	

		///Construct element level Right Hand Side
	    	///Access the connectivity of the element 'e'
	    	for(int i=0;i<3;i++)
			conn[i] = mesh->getElem(e)->getConn(i);

	    	///K[i][j]*T[j]
	    	for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
		    		RHS_e[i] = RHS_e[i] + mesh->getElem(e)->getK()[3*i+j] * mesh->getNode(conn[j])->getT();
			}	
	    	}

		///dt*(F + B - K*T)
		for(int i=0;i<3;i++){
			RHS_e[i] = dt*(mesh->getElem(e)->getF()[i] + mesh->getElem(e)->getB()[i] - RHS_e[i]);
	    	}

		///Access lumped mass matrix 
		std::memcpy(M_l,mesh->getElem(e)->getM(),3*sizeof(double));
		
		///M[3][3]*T[3] + dt(F + B - K*T)
		for(int i=0;i<3;i++)
			RHS_e[i] = RHS_e[i] + M_l[i]*mesh->getNode(conn[i])->getT();
	
		///Assemble global Diagonal Mass matrix and RHS
		for(int i=0;i<3;i++){
			M[conn[i]] = M[conn[i]] + M_l[i];
			RHS[conn[i]] = RHS[conn[i]] + RHS_e[i];
		}
		

	}///element loop end

	///Loop through all nodes, calculate and set the temperature (Also check if it reached steady state)
	double rate, max_rate = 0.0, T_prev, T_curr, T_max = 0;
	for(int node=0;node<nn;node++){	
		// Get previous Temperature of node
		T_prev = mesh->getNode(node)->getT();		
	
		///Set the calculated temperature to the nodes which are not on the Dirichlet Boundary
		if(mesh->getNode(node)->getBC_type()!=1)
			mesh->getNode(node)->setT(RHS[node]/M[node]);

		// Get current Temperature of node
		T_curr = mesh->getNode(node)->getT();

		// Calculate the square difference between previous and current temperature and add it to sum
		rate = fabs((T_curr - T_prev)/dt);
		if(rate>max_rate)	max_rate = rate;
		if(T_curr>T_max)	T_max = T_curr; 
	}
	

	if(max_rate<0.001){
		cout<<">> Solution reached Steady state! \n"<<endl;
		cout<<"> Maximum temperature in the domain: "<<T_max<<" K\ttime = "<<time<<" s\n"<<endl;
		break;
	}

	///Increase time by dt	
	time += dt;

    }///Time loop end

    delete postP;
    delete[] M;
    delete[] RHS;
    return;
}
