//=================================================================================================
// File Name   : constants.h
// Author      : 
// Version     : 1.0
// Copyright   : See the copyright notice in the README file.
// Description : You can define global constants here.
//=================================================================================================

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <limits>
#include <time.h>

using namespace std;

const double PI = 3.14159265;   /// Pi constant

const int nGQP = 7; /// Number of Gauss quadrature points

const int nsd = 2;  /// number of space dimensions
const int nen = 3;  /// number of element nodes
const int nef = 3;  /// number of element faces

///Mapping from the edge numbers to node numbers.
const int edgeNodes[3][2] = {{0,1},{1,2},{2,0}};

#endif /* CONSTANTS_H_ */
