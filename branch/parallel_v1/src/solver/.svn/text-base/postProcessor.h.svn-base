//==================================================================================================
// Name        : postProcessor.h
// Author      : 
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : Defines the postProcessor class which contains the post processing routines and
//               variables.
//==================================================================================================

#ifndef POSTPROCESSOR_H_
#define POSTPROCESSOR_H_

#include <sstream>

#include "solver.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkPropAssembly.h>
#include <vtkSmartPointer.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>

/*! \brief 
 * This class contains post processing routines and
 * variables.
 */

class postProcessor
{
    private:
        /// PRIVATE VARIABLES
        inputSettings*  settings;   // a local pointer for the settings
        triMesh*        mesh;       // a local pointer for the mesh
        double          minT;       // min value of the Temperature field
        double          maxT;       // max value of the Temperature field

        /// PRIVATE METHODS
        void evaluateLimits();
        void vtkVisualization(int ts, double time);
        // Here you can include your own postProcessing routine which creates the legacy VTK file
        // without using the VTK library.

    protected:

    public:
        /// DEFAULT CONSTRUCTOR
        postProcessor(){};

        /// DESTRUCTOR
        ~postProcessor(){};

        /// PUBLIC INTERFACE METHOD
        void postProcessorControl(inputSettings*, triMesh*, int, double);
};


#endif /* POSTPROCESSOR_H_ */
