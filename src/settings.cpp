//==================================================================================================
// Name        : settings.cpp
// Author      : 
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : This file contains the functions which are used for reading in the settings file.
//==================================================================================================
#include "settings.h"

//==================================================================================================
// void inputSettings::inputSettings() 
// Default constructor for an inputSetting object
//==================================================================================================
inputSettings::inputSettings()
{
    // Default values for the input parameters
    title = "Default Title";
    minfFile = "minf";
    mxyzFile = "mxyz";
    mienFile = "mien";
    mrngFile = "mrng";
    dataFile = "data";
    initT = 0.0;
    D = 1.0;
    source = 0.0;
    nIter = 1;
    dt = 1.0;
    dwf = 1;
    BC[0].BCType = 0;
    BC[0].BCValue = 0;
    BC[0].HTC = 0;
}

//==================================================================================================
// void inputSettings::readSettingsFile()
// Reads the settings file.
//==================================================================================================
void inputSettings::readSettingsFile()
{
    string lineString;
    string dummyString;
    char    dummyChar;
    
    ifstream inputFile;
    inputFile.open("settings.in",ios::in);
    if (inputFile.is_open()==false)
    {
        cout << "Unable to open input file! Aborting... " << endl;
        exit(0);
    }
    cout.precision(7);
    cout << scientific;
    cout << endl << "====== Settings ======" << endl;
    
    while (!inputFile.eof())
    {
        // Get a line and store in lineString
        getline(inputFile, lineString, '\n');
        // If the first character of the line is not a '#'
        if (lineString.c_str()[0] != '#')
        {
            istringstream iss(lineString);
            iss >> dummyString;
            if(dummyString == "title")
                iss >> title;
            else if(dummyString == "minf")
                iss >> minfFile;
            else if(dummyString == "mxyz")
                iss >> mxyzFile;
            else if(dummyString == "mien")
                iss >> mienFile;
            else if(dummyString == "mrng")
                iss >> mrngFile;
            else if(dummyString == "data")
                iss >> dataFile;
            else if(dummyString == "init")
                iss >> initT;
            else if(dummyString == "D")
                iss >> D;
            else if(dummyString == "S")
                iss >> source;
            else if(dummyString == "iter")
                iss >> nIter;
            else if(dummyString == "dt")
                iss >> dt;
            else if(dummyString == "dwf")
                iss >> dwf;
            else if(dummyString == "fg1")
            {
                iss >> BC[1].BCType;
                iss >> BC[1].BCValue;

		if(BC[1].BCType == 3)
		iss >> BC[1].HTC;
            }
            else if(dummyString == "fg2")
            {
                iss >> BC[2].BCType;
                iss >> BC[2].BCValue;

		if(BC[2].BCType == 3)
		iss >> BC[2].HTC;
            }
            else if(dummyString == "fg3")
            {
                iss >> BC[3].BCType;
                iss >> BC[3].BCValue;

		if(BC[3].BCType == 3)
		iss >> BC[3].HTC;
            }
            else if(dummyString == "fg4")
            {
                iss >> BC[4].BCType;
                iss >> BC[4].BCValue;

		if(BC[4].BCType == 3)
		iss >> BC[4].HTC;
            }
            else if(dummyString == "fg5")
            {
                iss >> BC[5].BCType;
                iss >> BC[5].BCValue;

		if(BC[5].BCType == 3)
		iss >> BC[5].HTC;
            }
            else if(dummyString == "fg6")
            {
                iss >> BC[6].BCType;
                iss >> BC[6].BCValue;
		iss >> BC[6].HTC;

		if(BC[6].BCType == 3)
		iss >> BC[6].HTC;
            }
            else
            {
                cout << endl << "Unknown keyword in the settings file : " << dummyString;
                cout << endl << "Aborting...";
                exit(0);    
            }

        }
        
    }

    // Report the settings read from the file.

    //cout << fixed;
    cout << "Title of the simualation                : " << title << endl;
    cout << "Name of the minf file                   : " << minfFile  << endl;
    cout << "Name of the mxyz file                   : " << mxyzFile  << endl;
    cout << "Name of the mien file                   : " << mienFile  << endl;
    cout << "Name of the mrng file                   : " << mrngFile  << endl;
    cout << "Name of the initial distribution file   : " << dataFile  << endl;
    cout << "Initial value of the dependent variable : " << initT  << endl;
    cout << "Diffusion coefficient                   : " << D     << endl;
    cout << "Source term                             : " << source     << endl;
    cout << "Number of maximum time steps            : " << nIter  << endl;
    cout << "Time step size                          : " << dt    << endl;
    cout << "Data Writing Frequency                  : " << dwf    << endl;
    cout << "Type and value of BC on FG1             : " << BC[1].BCType << " " << BC[1].BCValue << endl;
    cout << "Type and value of BC on FG2             : " << BC[2].BCType << " " << BC[2].BCValue << endl;
    cout << "Type and value of BC on FG3             : " << BC[3].BCType << " " << BC[3].BCValue << endl;
    cout << "Type and value of BC on FG4             : " << BC[4].BCType << " " << BC[4].BCValue << endl;
    cout << "Type and value of BC on FG5             : " << BC[5].BCType << " " << BC[5].BCValue << " " << BC[5].HTC << endl;
    cout << "Type and value of BC on FG6             : " << BC[6].BCType << " " << BC[6].BCValue << " " << BC[6].HTC << endl;
    cout << endl << endl;

    inputFile.close();
    
    return;
}


