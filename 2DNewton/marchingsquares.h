/*
 * marchingsquares.h: Header file for marchingsquares.cpp
 * Author: Patrick Butler
 * Created: 09/01/2010
 * Last Modified: 09/03/2010
 * 
 */

#ifndef _MARCHINGSQUARES_

#define _MARCHINGSQUARES_

#include <assert.h>
#include <iostream>
//#include <fstream>
using namespace std;
#include <math.h>

#include "Vector.h"
#include "polyfuncs.h"

int generateGrid (int incHorz, int incVert, FloatArrArr &gird);

int getGridVales (float upperBound, float lowerBound, 
				  float rightBound, float leftBound, int incHorz, int incVert, 
				  FloatArr fcoeff, V2iArr fdegs, FloatArrArr &grid);

int printGrid (FloatArrArr grid);

int findFirstSquare (float upperBound, float lowerBound, 
				     float rightBound, float leftBound, 
					 int incHorz, int incVert, FloatArrArr grid, V2f &startSquare);

int nextSquare(float upperBound, float lowerBound, 
			   float rightBound, float leftBound, 
			   int incHorz, int incVert, FloatArrArr grid, 
			   V2f &nextSquare, V2fArr &endPoints);

int generateEndPoints (float upperBound, float lowerBound, 
				       float rightBound, float leftBound, int incHorz, int incVert,
				       FloatArrArr grid, V2fArr &endPoints);

int marchingSquares (FloatArr fcoeff, V2iArr fdegs, float upperBound,
				     float lowerBound, float rightBound, float leftBound,
				     int incHorz, int incVert, V2fArr &endPoints);

//#ifndef WIN32
#include "marchingsquares.cpp"
//#endif

#endif
