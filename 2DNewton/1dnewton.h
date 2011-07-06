/*
 * 1dnewton.h: Header file for 1dnewton.cpp
 * Author: Patrick Butler
 * Created: 02/03/2010
 * Last Modified: 3/20/2010
 * 
 */

#ifndef _1DNEWTON_
#define _1DNEWTON_

//#include <assert.h>
//#include <iostream>
//using namespace std;
//#include <math.h>

#include "Vector.h"
#include "polyfuncs.h"


bool OneDNewton (FloatArr fcoeff, float guess, V2f interval, float &root);

//#ifndef WIN32
#include "1dnewton.cpp"
//#endif

#endif
