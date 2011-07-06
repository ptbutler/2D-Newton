/*
 * 2dnewton.h: Header file for 2dnewton.cpp
 * Author: Patrick Butler
 * Created: 02/10/2010
 * Last Modified: 3/20/2010
 * 
 */

#ifndef _2DNEWTON_
#define _2DNEWTON_

#include <assert.h>
#include <iostream>
//using namespace std;
#include <math.h>
#include <float.h>
#include <list>

//#include "Vector.h"
#include "polyfuncs.h"
#include "1dnewton.h"


void homogenize2D (FloatArr fcoeff, V2iArr fdegs, FloatArr &Fcoeff, V3iArr &Fdegs);

void BinomCoeff(int n, IntArr &bcoeff);
void polyToPower(float a, float b, int power, FloatArr &fcoeff);
void reduceCurve (FloatArr fcoeff, V2iArr fdegs, V2f dir, V2f point, FloatArr &gcoeff);
void findTangent (FloatArr fcoeff, V2iArr fdegs, V2f sp, V2f &tan);
bool curveInter (FloatArr fcoeff, V2iArr fdegs, FloatArr gcoeff, V2iArr gdegs,
		V2f startPoint, float length, V2f &endPoint);


void genSturm (FloatArr fcoeff, FloatArrArr &sturm);
int evalSturm (FloatArrArr sturm, float a);
int numOfRoots (V2f interval, FloatArrArr sturm);
void removeMultipleRoots(FloatArr fcoeff, FloatArr &gcoeff);
int numOfRoots (FloatArr fcoeff, V2f interval);
void getRootIntervals (V2f interval, FloatArrArr sturm, V2fArr &bd);
void getRootIntervals (FloatArr fcoeff, V2f interval, V2fArr &bd);
void getStartPoint (FloatArr fcoeff, V2iArr fdegs, V2f guess, V2f &start);
void findMaxInt(FloatArr p1dcoeff, V2f maxInt);
void getSpecialPoints (FloatArr fcoeff, V2iArr fdegs, V2f guess, V2f dir, V2fArr &specials);

void findCloseSingularityResultant (FloatArr fcoeff, V2iArr fdegs, V2f startPoint, V2f &singularity);
bool findCloseSingularityIntersect (FloatArr fcoeff, V2iArr fdegs, V2f startPoint, V2f &singularity);
void singularityToOrigin(FloatArr fcoeff, V2iArr fdegs, V2f singularity, FloatArr &hcoeff, V2iArr &hdegs);
void shiftStart(V2f start, V2f singularity, V2f &hstart);
void inverseToOrigin(V2fArr hpoints, V2f singularity, V2fArr &points);
void quadraticTransform(FloatArr fcoeff, V2iArr fdegs, FloatArr &gcoeff, V2iArr &gdegs);
void quadraticTransform (V2f start, V2f &gstart);
void inverseQuadraticTransform(V2fArr gpoints, V2fArr &fpoints);
void pointsThroughSingularity (FloatArr gcoeff, V2iArr gdegs, V2f gstart, float length, V2fArr &gpoints);
void jumpThroughSingularity (FloatArr fcoeff, V2iArr fdegs, V2f start, V2f singularity, float length, V2fArr &points);

bool isNearOrPastSingularity(FloatArr fcoeff, V2iArr fdegs, V2f lastPoint, V2f thisPoint, V2f tangent, bool &past);

void findNextPointSing(FloatArr fcoeff, V2iArr fdegs, V2f lastPoint, V2f startPoint, float length, int i, V2f &oldtangent, V2f &tangent, V2f &endPoint);
void findNextPoint(FloatArr fcoeff, V2iArr fdegs, V2f lastPoint, V2f startPoint, float length, int i, V2f &oldtangent, V2f &tangent, V2f &endPoint);
int TwoDNewton (FloatArr fcoeff, V2iArr fdegs, V2f guessPoint,
	       float length, int numPoints, V2fArrArr &endPoints);

//#ifndef WIN32
#include "2DNewton.cpp"
//#endif

#endif
