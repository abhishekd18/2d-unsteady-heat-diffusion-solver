//==================================================================================================
// Name        : decomposePar.h
// Author      : Abhishek Y. Deshmukh
// Version     : 1.3
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the solver class.
//==================================================================================================

#ifndef DECOMPOSEPAR_H_
#define DECOMPOSEPAR_H_

#include "settings.h"
#include "tri.h"

/*!
 * \brief This class defines the decomposition methods.
 */

class decomposePar
{
    private:
        /// PRIVATE VARIABLES
	triMesh* 	meshL;	    // a local pointer for the mesh
	int* 		nel;

        /// PRIVATE METHODS
	void generateProcB(const int, inputSettings*);

    protected:

    public:
        /// DEFAULT CONSTRUCTOR
        decomposePar(){};

        /// DESTRUCTOR
        ~decomposePar(){};

        /// INTERFACE FUNCTION
        void decompose(inputSettings*, triMesh*);
};

#endif /* SOLVER_H_ */


