//==================================================================================================
// Name        : solver.h
// Author      : Raghavan Lakshmanan
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the solver class.
//==================================================================================================

#ifndef SOLVER_H_
#define SOLVER_H_

#include "settings.h"
#include "tri.h"

//==================================================================================================
// Solver class
//==================================================================================================
class femSolver
{
    private:
        /// PRIVATE VARIABLES
        inputSettings*  settings;   // a local pointer for the settings
        triMesh*        mesh;       // a local pointer for the mesh
        int my_rank;
        int num_procs;
        int prev;
        int next; 

        /// PRIVATE METHODS
        void calculateJacobian(const int);
        void calculateElementMatrices(const int);
        void applyBoundaryConditions(const int);
        void explicitSolver();
    

    protected:

    public:
        /// DEFAULT CONSTRUCTOR
        femSolver(){};

        /// DESTRUCTOR
        ~femSolver(){};

        /// INTERFACE FUNCTION
        void solverControl(inputSettings*, triMesh*, int, int, int, int);
    
};

#endif /* SOLVER_H_ */


