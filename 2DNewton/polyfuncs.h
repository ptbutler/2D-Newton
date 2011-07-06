/*
 * polyfunc.h: Header file for polyfunc.cpp
 * Author: Patrick Butler
 * Created: 09/02/2010
 * Last Modified: 09/02/2010
 * 
 */

#ifndef _POLYFUNCS_
#define _POLYFUNCS_

#include <assert.h>
#include <iostream>
//using namespace std;
#include <math.h>

#include "Vector.h"

float polyEval (FloatArr fcoeff, float a);
float evalBi (FloatArr fcoeff, V2iArr fdegs, V2f point);
float evalTri (FloatArr fcoeff, V3iArr fdegs, V3f point);

void addPoly (FloatArr fcoeff, FloatArr gcoeff, FloatArr &hcoeff);
void subPoly (FloatArr fcoeff, FloatArr gcoeff, FloatArr &hcoeff);
void mulPolyByCoeff (FloatArr &fcoeff, float factor);
void mulPoly(FloatArr fcoeff, FloatArr gcoeff, FloatArr &hcoeff);
void dividePoly (FloatArr fcoeff, FloatArr gcoeff, FloatArr &q, FloatArr &r);
void derive (FloatArr fcoeff, FloatArr &fpcoeff);
void polyGCD(FloatArr fcoeff, FloatArr gcoeff, FloatArr &gcd);

void polyAxisInteresect (FloatArr Fcoeff, V3iArr Fdegs, int axis, float k, FloatArr &fcoeff, V2iArr &fdegs);
void polyPlaneIntersect (FloatArr Fcoeff, V3iArr Fdegs, V4f plane, FloatArr &fcoeff, V2iArr &fdegs);

void pDerive2D (FloatArr fcoeff, V2iArr fdegs, int var, FloatArr &fpcoeff, V2iArr &fpdegs);
void pDerive3D (FloatArr Fcoeff, V3iArr Fdegs, int var, FloatArr &Fpcoeff, V3iArr &Fpdegs);

void mulPolyToPoly(FloatArr xcoeff, FloatArr ycoeff, FloatArr &mulcoeff, V2iArr &muldegs);

//#ifndef WIN32
#include "polyfuncs.cpp"
//#endif

#endif
