//==================================================================================================
// Name        : solver.h
// Author      : A. Emre Ongut
// Version     : 1.3
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the solver class.
//==================================================================================================

#ifndef SOLVER_H_
#define SOLVER_H_

#include "settings.h"
#include "tri.h"

/*!
 * \brief This class defines the solver control and solver member functions
 */

class femSolver
{
    private:
        /// PRIVATE VARIABLES
        inputSettings*  settings;   // a local pointer for the settings
        triMesh*        mesh;       // a local pointer for the mesh

        /// PRIVATE METHODS
        void calculateJacobian(const int);
        void calculateElementMatrices(const int);
        void applyBoundaryConditions(const int);
        void globalAssembly();
        void explicitSolver();

    protected:

    public:
        /// DEFAULT CONSTRUCTOR
        femSolver(){};

        /// DESTRUCTOR
        ~femSolver(){};

        /// INTERFACE FUNCTION
        void solverControl(inputSettings*, triMesh*);

};

#endif /* SOLVER_H_ */


