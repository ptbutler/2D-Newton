/*
 * 1dnewton.cpp: Finds the closest root of a polynomial given a point
 * Author: Patrick Butler
 * Created: 01/28/2010
 * Last Modified: 03/20/2010
 *
 */

#include "1dnewton.h"



/******************************************************************************
 * Finds the root of a polynomial f(x) using the 1-D Newton method, which is
 * an iterative process. If it converges, the function returns true, otherwise, 
 * it returns false.
 *
 * ->fcoeff: array of coefficients of the function in which the roots are to be found
 * ->guess: place at which to start the process, a good guess close to a root
 * ->interval: the interval in which the root is known to be
 * <-root: the root
 * <-: the answer is valid or not
 ******************************************************************************/

bool OneDNewton (FloatArr fcoeff, float guess, V2f interval, float &root)
{
  float x = guess;
  float xOld = x;
  float relErr = 1;
  float relErrOld = 1;
  float relEps = 0.000001;
  float eps =    0.000001;   // chosen as floating point is understood to have ~6 digits of precision
  int maxIter = 1000;
  int iter = 0;
  FloatArr fpcoeff;

  float fx;
  fx = polyEval(fcoeff, x);
  derive(fcoeff, fpcoeff);


  while (fabs(fx) > eps &&      // continue until sufficiently accurate
	 relErr > relEps &&   
	 iter < maxIter)           // or taking too long
  {                                              
    float fpx = polyEval(fpcoeff, x);

    // special case, lines parallel to x-axis are bad
    if (fabs(fpx) <= eps)
    {
#ifdef _VERBOSEMODE_
      cout<<"Error: derivative of 0!"<<endl;
#endif
      root = x;
      return false;
    }

    // keep a record of the previous x
    xOld = x;
    x = x - (fx / fpx); 

    // did we get better by a diminishing amount?
    relErrOld = relErr;
    relErr = fabs((xOld - x) / xOld);

//     // if not, we are diverging
//     if(relErr > relErrOld)
//     {
//       root = x;
//       return false;
//     }

    // if we jump out of the given bounds, we are diverging
    if (x < interval[0] || x > interval[1])
    {
#ifdef _VERBOSEMODE_
      cout<<"Error: jumped out of bounds!"<<endl;
#endif
      return false;
    }

    fx = polyEval(fcoeff, x);
    iter++;
  }
  root = x;

  if (iter == maxIter) // if we took too long
  {
#ifdef _VERBOSEMODE_
    cout<<"Error: took too long!"<<endl;
#endif
    return false;      // answer might not be correct
  }
  else
    return true;
}
