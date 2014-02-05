//==================================================================================================
// Name        : postProcessor.cpp
// Author      : Raghavan Lakshmanan
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the necessary routines for post processing the results.
//==================================================================================================

#include "postProcessor.h"

//==================================================================================================
// preProcessorControl
// Public interface for post process control. It calls other post processing routines.
//==================================================================================================
void postProcessor::postProcessorControl(inputSettings* argSettings, triMesh* argMesh, int ts, double time)
{
        cout << endl << "===== Post-processing =====" << endl;

        mesh = argMesh;
        settings = argSettings;
        
        evaluateLimits();
        vtkVisualization(ts,time);
        
        return;
}

//==================================================================================================
// Evaluates the maximum and minimum temperatures in the field
//==================================================================================================
void postProcessor::evaluateLimits()
{
        int nn = mesh->getNn();
        double T;

        //Let's get the largest and smallest numbers possible.
        minT = std::numeric_limits<double>::max();
        maxT = std::numeric_limits<double>::min();

        //Find max and min.
       
      
        for(int i=0; i<nn; i++)
        {
        T = mesh->getNode(i)->getT();
                if(T < minT)
                        minT = T;
                if(T > maxT)
                        maxT = T;
        }
        cout << "Tmin" << minT << endl;
        cout << "Tmax" << maxT << endl;

        return;
}

//==================================================================================================
// Main visualization function
//==================================================================================================
void postProcessor::vtkVisualization(int ts,double time)
{
        int nn = mesh->getNn();
        int ne = mesh->getNe();
        int ne_pro = mesh->getNe_pro(); // Get no.of elements in each processor
        string dummy;

        // VTK Double Array
        vtkDoubleArray* pcoords = vtkDoubleArray::New();
        pcoords->SetNumberOfComponents(nen);
        pcoords->SetNumberOfTuples(nn);

        // vtkDoubleArray type pcoords is filled with the data in meshPoints.
        for (int i=0; i<nn; i++)
          pcoords->SetTuple3(i,mesh->getNode(i)->getX(),mesh->getNode(i)->getY(),0.0f);

        // vtkPoints type outputPoints is filled with the data in pcoords.
        vtkPoints* outputPoints = vtkPoints::New();
        outputPoints->SetData(pcoords);

        // Connectivity is written to vtkCellArray type outputCells
        vtkCellArray* connectivity = vtkCellArray::New();
        for(int i=0; i<ne_pro; i++)
        {
            connectivity->InsertNextCell(nen);
            for(int j=0; j<nen; j++)
                connectivity->InsertCellPoint(mesh->getElem(i)->getConn(j));
        }

        // Scalar property
        vtkDoubleArray* pressure = vtkDoubleArray::New();
        pressure->SetName("Temperature");
       
        for(int i=0; i<nn; i++) 
            pressure->InsertNextValue(mesh->getNode(i)->getT());

        // Previously collected data which are outputPoints, outputCells, scalarProperty, are written to
        // vtkPolyData type polydata var.
        vtkPolyData* polydata = vtkPolyData::New();
        polydata->SetPoints(outputPoints);
        polydata->SetPolys(connectivity);
        polydata->GetPointData()->SetScalars(pressure);

        // Whatever collected in polydata above is written to the vtk file below.
        // vtkDataSetWriter is for leagacy VTK format, vtkXMLDataSetWriter is for VTK XML format.
        vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
        // vtkXMLDataSetWriter *writer = vtkXMLDataSetWriter::New();
   
        dummy = settings->getWdir();                                         
   	std::ostringstream myrank_str; 	// output string stream              
	myrank_str << MPI::COMM_WORLD.Get_rank();                              
	dummy.append("proc_").append(myrank_str.str()).append("/");           

        dummy.append(settings->getTitle());
        dummy.append(".");
    	ostringstream ss; ss << ts;
        dummy.append(ss.str());
        dummy.append(".vtk");
        //cout<<dummy<<endl;
        
       /// Add time to vtkDataSet
    	vtkDoubleArray *t = vtkDoubleArray::New();
    	t->SetName("TIME");
    	t->SetNumberOfTuples(1);
    	t->SetTuple1(0, time);
    	polydata->GetFieldData()->AddArray(t);


        writer->SetFileName(dummy.c_str());
        writer->SetInput(polydata);
        writer->Write();

 /*       // In the below section, data collected in polydata var is rendered.
        // Create the mapper and set the appropriate scalar range. (default is (0,1)
        vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
        mapper->SetInput(polydata);
        mapper->SetScalarRange(minT, maxT);

        // Create an actor.
        vtkActor* actor = vtkActor::New();
        actor->SetMapper(mapper);

        // Scalar bar
        vtkSmartPointer<vtkScalarBarActor> scalarBar = 
            vtkSmartPointer<vtkScalarBarActor>::New();
        scalarBar->SetLookupTable(mapper->GetLookupTable());
        scalarBar->SetTitle("Temperature");
        scalarBar->SetNumberOfLabels(10);
 
        // Create a lookup table to share between the mapper and the scalarbar
        vtkSmartPointer<vtkLookupTable> hueLut =
                vtkSmartPointer<vtkLookupTable>::New();
        hueLut->SetTableRange (0, 1);
        hueLut->SetHueRange (0.7, 0.0);
        hueLut->SetSaturationRange (1, 1);
        hueLut->SetValueRange (1, 1);
        hueLut->Build();

        // llokup table shared between mapper and scalarbar
        mapper->SetLookupTable( hueLut );
        scalarBar->SetLookupTable( hueLut );

        // Create the rendering objects.
        vtkRenderer* renderer = vtkRenderer::New();
        renderer->GradientBackgroundOn();
        renderer->SetBackground(1,1,1);
        renderer->SetBackground2(0,0,0);

        vtkRenderWindow* renWin = vtkRenderWindow::New();
        renWin->AddRenderer(renderer);
        renWin->SetSize(1000,600);

        cout << "> Rendering window is initiated." << endl;
        vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
        iren->SetRenderWindow(renWin);

        vtkSmartPointer<vtkAxesActor> axes = 
            vtkSmartPointer<vtkAxesActor>::New();
    
        vtkSmartPointer<vtkOrientationMarkerWidget> widget = 
                vtkSmartPointer<vtkOrientationMarkerWidget>::New();
        widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
        widget->SetOrientationMarker( axes );
        widget->SetInteractor( iren );
        widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
        widget->SetEnabled( 1 );
        widget->InteractiveOff();

        renderer->AddActor(actor);
        renderer->AddActor2D(scalarBar);
 
        renderer->ResetCamera();
        renWin->Render();
 
        // Begin mouse interaction
        iren->Start();*/

        // Unnecessary data is destroyed.
        pcoords->Delete();
        outputPoints->Delete();
        connectivity->Delete();
        pressure->Delete();
        polydata->Delete();
/*        mapper->Delete();
        actor->Delete();
        renderer->Delete();
        renWin->Delete();
        iren->Delete();*/

        return;
}
