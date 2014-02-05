//==================================================================================================
// Name        : tri.h
// Author      : Raghavan Lakshmanan
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : Class definitions for triangular elements. There are different classes for storing
//               variables at node and element level.
//==================================================================================================

#ifndef TRI_H_
#define TRI_H_

#include "settings.h"
#include "mpi.h"

//==================================================================================================
// class triMasterElement == > LINEAR TRIANGULAR MASTER ELEMENT
//==================================================================================================
// In the mesh, there are elements with different sizes and shapes. Instead of integrating the 
// equations on each element, we use a reference element, which simplifies the integration.
// Note that the reference element is same for all the elements in the mesh.
//==================================================================================================
class triMasterElement
{
    private:
        /// PRIVATE VARIABLES
        double point[2];    // ksi and eta for each GQ point
        double weight;      // weight of each GQ point
        double S[3];        // Shape functions
        double dSdKsi[3];   // ksi derivatives of shape functions
        double dSdEta[3];   // eta derivatives of shape functions
        
    protected:

    public:
        /// GETTERS
        double getPoint(int i)  {return point[i];};
        double getWeight()      {return weight;};
        double getS(int i)      {return S[i];};
        double getDSdKsi(int i) {return dSdKsi[i];};    
        double getDSdEta(int i) {return dSdEta[i];};    

        /// PUBLIC INTERFACE METHODS
        void setupGaussQuadrature();
        void evaluateShapeFunctions();
    
};


//==================================================================================================
// class trNode ==> NODE LEVEL DATA STRUCTURE
//==================================================================================================
// This data structure holds any kind of variable that are defined at node points.
// These variables can be x,y coordinates, temperature values or any other value that is necessary
// during the solution stage. 
//==================================================================================================
class triNode
{
    private:
        /// PRIVATE VARIABLES
        double x;         // x-coordinate
        double y;         // y-ccordinate
        double T;         // Temperature
        int flag;         // Flag nodes which are in dirichilet boundary - set BCtype
        int common_node;  // Flag nodes which are shared by other processors
        // maybe some more variables can be necessary for the solver ?

    protected:

    public:
        /// DEFAULT CONSTRUCTOR
        triNode(){x=0.0f; y=0.0f; T=0.0f;flag=0,common_node=0;};

        /// DESTRUCTOR
        ~triNode(){};

        /// SETTERS
        void setX      (double value) {x = value;};
        void setY      (double value) {y = value;};
        void setT      (double value) {T = value;};
        void set_flag  (int value)    {flag=value;}; 
        void set_common(int value)    {common_node=value;}; 


        /// GETTERS
        double getX() {return x;};
        double getY() {return y;};
        double getT() {return T;};
        int    get_flag() {return flag;};
        int    get_common() {return common_node;};
};


//==================================================================================================
// class triElement ==> ELEMENT LEVEL DATA STRUCTURE
//==================================================================================================
// This data structure holds the variables determined for an element.
// An element have a connectivity (the node numbers of the element) and face groups which indicates
// which face has which boundary condition type. Other than these you will need to define element
// mass stiffness matrices or any other matrices that will be necessary during the solution stage.
//==================================================================================================
class triElement
{
    private:
        /// PRIVATE VARIABLES
        int conn[nen];
        int FG[nef];
        // and some more variables that you will need during the solution stage...
        double J[2][2];
        double detJ;
        double J_inv[2][2];
        double ele_mat[3][3];
        double ele_mass[3][3];
        double ele_lum_mass[3];
        double ele_flux[3];
        double B[3];

    protected:

    public:
        /// DEFAULT CONSTRUCTOR
        triElement(){conn[0]=0;conn[1]=0;conn[2]=0;};

        /// DESTRUCTOR
        ~triElement(){};

        /// SETTERS
        void setConn    (int i, int value) {conn[i] = value;};
        void setFG      (int i, int value) {FG[i] = value;};
        void setJ       (int i, int j, double value) {J[i][j]=value;};
        void setdetJ    (double value) {detJ=value;};
        void setJ_inv   (int i, int j, double value) {J_inv[i][j]=value;}; 
        void setele_mat (int i, int j, double value) {ele_mat[i][j] = value;};
        void setele_mass(int i, int j, double value) {ele_mass[i][j] = value;};
        void setele_lum_mass(int i, double value) {ele_lum_mass[i] = value;};
        void setele_flux(int i, double value) {ele_flux[i] = value;};
        void setB      (int i, double value) {B[i]=value;};

        /// GETTERS
        int     getConn (int index) {return conn[index];};
        int     getFG   (int index) {return FG[index];};
        double  getJ    (int i,int j) {return J[i][j];};
        double  getdetJ () {return detJ;};        
        double  getJ_inv(int i,int j) {return J_inv[i][j];};
        double  getele_mat (int i,int j) {return ele_mat[i][j];};
        double  getele_mass (int i,int j) {return ele_mass[i][j];};
        double  getele_lum_mass (int i) {return ele_lum_mass[i];};
        double  getele_flux (int i) {return ele_flux[i];};
        double  getB    (int i) {return B[i];};

};


//==================================================================================================
// class triMesh ==> MESH DATA STRUCTURE
//==================================================================================================
// This class is used to keep the pointers to node and element data structures. It is more
// convinient to pass the pointer to this class rather than passin the pointers to both element and
// node arrays during function calls. So it is just for simplification.
//==================================================================================================
class triMesh
{
    private:
        /// PRIVATE VARIABLES
        int ne;                     // total number of elements
        int nn;                     // total number of nodes
        int ne_pro;                 // number of elements per processor
        int nn_pro;                 // number of nodes per processor
        int element_index;          // starting element index in each processor
        int node_index;             // starting node number in each processor
        triNode*            node;   // pointer for node level data structure
        triElement*         elem;   // pointer for element level data structure
        triMasterElement*   ME;     // pointer for reference element
        int my_rank;                // rank of the processor
        int num_procs;              // total no.of processors
        int prev;                   // rank of previous processor
        int next;                   // rank of next processor
       

        /// PRIVATE METHODS
        void swapBytes(char*, int, int);
        
    protected:

    public:
        /// DEFAULT CONSTRUCTOR
        triMesh(){};

        /// DESTRUCTOR       
        ~triMesh()
        {
            delete[] node;
            delete[] elem;
            delete[] ME;
        };

       
        /// GETTERS
        int                 getNe()             {return ne;};
        int                 getNn()             {return nn;};
        int                 getNe_pro()         {return ne_pro;};
        int                 getNn_pro()         {return nn_pro;};
        int                 getnode_index()     {return node_index;};
        int                 getelem_index()     {return element_index;};

        triNode*            getNode (int index) {return &node[index];};
        triElement*         getElem (int index) {return &elem[index];};
        triMasterElement*   getME   (int index) {return &ME[index];};

        /// PUBLIC INTERFACE METHOD
        void readMeshFiles(inputSettings*, int, int, int, int);
  
};


#endif /* TRI_H_ */
