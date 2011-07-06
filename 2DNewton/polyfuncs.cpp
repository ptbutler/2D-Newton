/*
 * marchingsquares.cpp: Alternative method of generating a point set
 * Author: Patrick Butler
 * Created: 09/01/2010
 * Last Modified: 09/03/2010
 *
 */


#include "polyfuncs.h"


/******************************************************************************
 * Evaluates the polynomial function f(x) described by for a supplied a 
 * using Horner's scheme
 *
 * ->fcoeff: array of coefficients of f(x)
 * ->a: 
 * <-: f(a)
 ******************************************************************************/

float polyEval (FloatArr fcoeff, float a)
{
  float answer = 0;
  
  for(int i = fcoeff.getn() - 1; i >= 0; i--)
    answer = answer * a + fcoeff[i];

  return answer;
}


/************************************************************************
 * Evalutes a bivariate polynomial f(x,y) at a point (a,b)
 * 
 * ->Fcoeff: the coefficients of f(x,y)
 * ->Fdegs: the degrees of f(x,y)
 * ->point: the point (a,b)
 * <-: the value of f(x,y) at (a,b)
 ************************************************************************/

float evalBi (FloatArr fcoeff, V2iArr fdegs, V2f point)
{
  float answer = 0;
  for (int i = 0; i < fcoeff.getn(); i++)
    answer += fcoeff[i] * pow(point[0], fdegs[i][0]) 
                        * pow(point[1], fdegs[i][1]);

  return answer;
}


/************************************************************************
 * Evalutes a trivariate polynomial f(x,y,z) at a point (a,b,c)
 * 
 * ->Fcoeff: the coefficients of f(x,y,z)
 * ->Fdegs: the degrees of f(x,y,z)
 * ->point: the point (a,b,c)
 * <-: the value of f(x,y,z) at (a,b,c)
 ************************************************************************/

float evalTri (FloatArr fcoeff, V3iArr fdegs, V3f point)
{
  float answer = 0;
  for (int i = 0; i < fcoeff.getn(); i++)
    answer += fcoeff[i] * pow(point[0], fdegs[i][0]) 
                        * pow(point[1], fdegs[i][1]) 
                        * pow(point[2], fdegs[i][2]);

  return answer;
}


/************************************************************************
 * Adds two polynomials f(x) & g(x) together to produce another 
 * polynomial h(x)
 *
 * ->fcoeff: the coefficients of f(x)
 * ->gcoeff: the coefficients of g(x)
 * <-hcoeff: the coefficients of h(x)
 *
 * Assumption: all coefficients are of the terms in increasing sequential
 * order
 ************************************************************************/

void addPoly (FloatArr fcoeff, FloatArr gcoeff, FloatArr &hcoeff)
{
  if (fcoeff.getn() > gcoeff.getn())
    hcoeff.allocate(fcoeff.getn());
  else
    hcoeff.allocate(gcoeff.getn());

  
  for (int i = 0; i < hcoeff.getn(); i++)
  {
    if (i < fcoeff.getn() && i < gcoeff.getn())
      hcoeff[i] = fcoeff[i] + gcoeff[i];
    else if (i < fcoeff.getn())
      hcoeff[i] = fcoeff[i];
    else
      hcoeff[i] = gcoeff[i];
  }
}


/************************************************************************
 * Subtracts the polynomial g(x) from f(x) to produce another 
 * polynomial h(x)
 *
 * ->fcoeff: the coefficients of f(x)
 * ->gcoeff: the coefficients of g(x)
 * <-hcoeff: the coefficients of h(x)
 *
 * Assumption: all coefficients are of the terms in increasing sequential
 * order
 ************************************************************************/

void subPoly (FloatArr fcoeff, FloatArr gcoeff, FloatArr &hcoeff)
{
  if (fcoeff.getn() > gcoeff.getn())
    hcoeff.allocate(fcoeff.getn());
  else
    hcoeff.allocate(gcoeff.getn());

  
  for (int i = 0; i < hcoeff.getn(); i++)
  {
    if (i < fcoeff.getn() && i < gcoeff.getn())
      hcoeff[i] = fcoeff[i] - gcoeff[i];
    else if (i < fcoeff.getn())
      hcoeff[i] = fcoeff[i];
    else
      hcoeff[i] = -gcoeff[i];
  }
}


/************************************************************************
 * Multiplys two univariate polynomials f(x) & g(x) to get h(x)
 * 
 * ->fcoeff: the coefficients of f(x)
 * ->gcoeff: the coefficients of g(x)
 * <-hcoeff: the coefficients of h(x)
 *
 * Assumption: all coefficients are of the terms in increasing sequential
 * order
 ************************************************************************/

void mulPoly (FloatArr fcoeff, FloatArr gcoeff, FloatArr &hcoeff)
{
  hcoeff.allocate(fcoeff.getn() + gcoeff.getn() - 1);
  int i;
  for (i = 0; i < hcoeff.getn(); i++)
    hcoeff[i] = 0;

  for (i = 0; i < fcoeff.getn(); i++)
    for (int j = 0; j < gcoeff.getn(); j++)
      hcoeff[i+j] += fcoeff[i] * gcoeff[j];
}


/************************************************************************
 * Multiplys a polynomial f(x) by a factor 
 * 
 * <-fcoeff: the coefficients of f(x)
 * ->factor: the factor to multiply by
 *
 * Assumption: all coefficients are of the terms in increasing sequential
 * order
 ************************************************************************/

void mulPolyByCoeff (FloatArr &fcoeff, float factor)
{
  for (int i = 0; i < fcoeff.getn(); i++)
    fcoeff[i] = fcoeff[i] * factor;
}


/************************************************************************
 * Divides one unviariate polynomial f(x) by another univariate
 * polynomial g(x), giving quotient q(x) and remainder r(x)
 * 
 * ->fcoeff: the coefficients of f(x)
 * ->gcoeff: the coefficients of g(x)
 * <-q: the coefficients of q(x)
 * <-r: the coefficients of r(x)
 ************************************************************************/

void dividePoly (FloatArr fcoeff, FloatArr gcoeff, FloatArr &q, FloatArr &r)
{
  int n = fcoeff.getn() - gcoeff.getn() + 1;
  q.allocate(n);
  
  int i;
  for (i = 0; i < n; i++)
    q[i] = 0;

  for (i = n-1; i >= 0; i--)
  {
    q[i] = fcoeff[i + gcoeff.getn() - 1] / gcoeff[gcoeff.getn() - 1];
    FloatArr scoeff;
    FloatArr dcoeff;
    dcoeff.allocate(i+1);
    int j;
    for (j = 0; j < dcoeff.getn(); j++)
      dcoeff[j] = 0;

    dcoeff[i] = q[i];
    mulPoly(gcoeff, dcoeff, scoeff);

    FloatArr tmpcoeff;
    subPoly(fcoeff, scoeff, tmpcoeff);
    
    fcoeff.allocate(fcoeff.getn() - 1);

    for (j = 0; j < fcoeff.getn(); j++)
      fcoeff[j] = tmpcoeff[j];
  }
  r = fcoeff;

  for (i = r.getn(); i > 0; i--)
  {
    if(r[i-1] != 0)
      break;
    
    r.shrink(r.getn()-1);
  }
}


/******************************************************************************
 * Derives the polynomial f(x) and returns the derivative f'(x)
 *
 * ->fcoeff: array of coefficients of f(x)
 * <-fpcoeff: array of coefficients of f'(x)
 *
 * Assumption: coefficients are in increasing order, and all are there
 ******************************************************************************/

void derive (FloatArr fcoeff, FloatArr &fpcoeff)
{
  if (fcoeff.getn() < 2)
  {
    fpcoeff.allocate(1);
    fpcoeff[0] = 0;
	return;
  }

  fpcoeff.allocate(fcoeff.getn()-1);

  for (int i = 0; i < fpcoeff.getn(); i++)
    fpcoeff[i] = fcoeff[i+1] * (i + 1);
}


/******************************************************************************
 * Find the greatest common divisor of f(x) and g(x) using Euler's method
 *
 * ->fcoeff: array of coefficients of f(x)
 * ->gcoeff: array of coefficients of g(x)
 * <-gcd: array of coefficients of the gcd(x)
 *
 * Assumption: coefficients are in increasing order, and all are there
 ******************************************************************************/

void polyGCD(FloatArr fcoeff, FloatArr gcoeff, FloatArr &gcd)
{
	while(true)
	{
		FloatArr q, rem;
		dividePoly(fcoeff, gcoeff, q, rem);
	
		int zeros = 0;
		for (int i = 0; i < rem.getn(); i++)
			if (rem[i] == 0)
				zeros++;

		if (zeros == rem.getn())
		{
			gcd = gcoeff;
			return;
		}
		fcoeff = gcoeff;
		gcoeff = rem;
	}
}


/************************************************************************
 * Generates a bivariate polynomial f from a trivariate polynomial F and 
 * an axis at a given value k
 *
 * ->Fcoeff: the coefficients of F
 * ->Fdegs: the powers of the variables of F
 * ->axis: the axis
 * ->k: the value
 * <-fcoeff: the coefficients of f
 * <-fdegs: the powers of the variables of f
 *
 ************************************************************************/

void polyAxisIntersect (FloatArr Fcoeff, V3iArr Fdegs, int axis, float k, FloatArr &fcoeff, V2iArr &fdegs)
{
	int altAxis1;
	int altAxis2;
	// set up the variables of the new polynomial
	if(axis == 0)
	{
		altAxis1 = 1;
		altAxis2 = 2;
	}
	else if(axis == 1)
	{
		altAxis1 = 0;
		altAxis2 = 2;
	}
	else
	{
		altAxis1 = 0;
		altAxis2 = 1;
	}

	int terms = Fcoeff.getn();
	
	fdegs.allocate(terms);
	fcoeff.allocate(terms);
	// initialize the values
	for(int i = 0; i < terms; i++)
		fcoeff[i] = 0;
	
	int count = 0;
	// evaluate the polynomial at the value of the axis
	for(int i = 0; i < Fcoeff.getn(); i++)
	{
		float coeff  = Fcoeff[i] * pow(k, Fdegs[i][axis]);
		V2i degs;
		degs[0] = Fdegs[i][altAxis1];
		degs[1] = Fdegs[i][altAxis2];
		bool same = false;
		for(int j = 0; j < count; j++)
		{
			if(fdegs[j] == degs)  // combine like terms
			{
				fcoeff[j] += coeff;
				same = true;
				break;
			}
		}
		if(!same)
		{
			fcoeff[count] = coeff;
			fdegs[count] = degs;
			count++;
		}
	}
	// shrink to distinct terms
	fcoeff.shrink(count);
	fdegs.shrink(count);

}


/************************************************************************
 * Generates a bivariate polynomial f from a trivariate polynomial
 * F(x) and a plane
 *
 * ->Fcoeff: the coefficients of F
 * ->Fdegs: the powers of the variables of F
 * ->plane: the plane, in the form of: 
 *          plane[0]x + plane[1]y + plane[2]z + plane[3] = 0
 * <-fcoeff: the coefficients of f
 * <-fdegs: the powers of the variables of f
 *
 ************************************************************************/

void polyPlaneIntersect (FloatArr Fcoeff, V3iArr Fdegs, V4f plane, FloatArr &fcoeff, V2iArr &fdegs)
{
	// solve plane for z, y or x, must be non-zero coefficient
	FloatArr gcoeff;
	V2iArr gdegs;
	gcoeff.allocate(3);
	gdegs.allocate(3);
	if(plane[2] != 0)
	{
		gcoeff[0] = plane[0]/plane[2];
		gdegs[0][0] = 1;
		gdegs[0][1] = 0;
		gcoeff[1] = plane[1]/plane[2];
		gdegs[0][0] = 0;
		gdegs[0][1] = 1;
		gcoeff[3] = plane[3]/plane[2];
		gdegs[0][0] = 0;
		gdegs[0][1] = 0;
	}
	else if(plane[1] != 0)
	{
		gcoeff[0] = plane[0]/plane[1];
		gdegs[0][0] = 1;
		gdegs[0][1] = 0;
		gcoeff[2] = plane[2]/plane[1];
		gdegs[0][0] = 0;
		gdegs[0][1] = 1;
		gcoeff[3] = plane[3]/plane[1];
		gdegs[0][0] = 0;
		gdegs[0][1] = 0;
	}
	else if(plane[0] != 0)
	{
		gcoeff[1] = plane[1]/plane[0];
		gdegs[0][0] = 1;
		gdegs[0][1] = 0;
		gcoeff[2] = plane[2]/plane[0];
		gdegs[0][0] = 0;
		gdegs[0][1] = 1;
		gcoeff[3] = plane[3]/plane[0];
		gdegs[0][0] = 0;
		gdegs[0][1] = 0;
	}
	else
	{
#ifdef _VERBOSEMODE_
		cout<<"Error: This is not a plane!"<<endl;
#endif
		return;
	}

	// substitute plane for z in F


	// reduce F to f
}


/************************************************************************
* Finds the partial derivative fp(x,y) of the bivariate polynomial 
* f(x,y) with respect to the given variable var
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->var: the variable to take the partial derivative with respect to
* <-fpcoeff: the coefficients of fp(x,y)
* <-fpdegs: the degrees of fp(x,y)
************************************************************************/

void pDerive2D (FloatArr fcoeff, V2iArr fdegs, int var, FloatArr &fpcoeff, V2iArr &fpdegs)
{
	if (var > 1 || var < 0)
		return;

	int var2;
	if (var == 0)
		var2 = 1;
	else
		var2 = 0;

	int diff = 0;
	int i;
	for (i = 0; i < fdegs.getn(); i++)
		if(fdegs[i][var] == 0)
			diff++;

	fpcoeff.allocate(fcoeff.getn() - diff);
	fpdegs.allocate(fdegs.getn() - diff);

	int j = 0;
	for (i = 0; i < fcoeff.getn(); i++)
	{
		if (fdegs[i][var] != 0)
		{
			fpcoeff[j] = fcoeff[i] * fdegs[i][var];
			fpdegs[j][var] = fdegs[i][var] - 1;
			fpdegs[j][var2] = fdegs[i][var2];
			j++;
		}
	}
}


/************************************************************************
* Finds the partial derivative fp(x,y,z) of the trivariate polynomial 
* f(x,y,z) with respect to the given variable var
* 
* ->fcoeff: the coefficients of f(x,y,z)
* ->fdegs: the degrees of f(x,y,z)
* ->var: the variable to take the partial derivative with respect to
* <-fpcoeff: the coefficients of fp(x,y,z)
* <-fpdegs: the degrees of fp(x,y,z)
************************************************************************/

void pDerive3D (FloatArr Fcoeff, V3iArr Fdegs, int var, FloatArr &Fpcoeff, V3iArr &Fpdegs)
{
	if (var > 2 || var < 0)
		return;

	int otherVar1, otherVar2;
	if (var == 0)
	{
		otherVar1 = 1;
		otherVar2 = 2;
	}
	else if (var == 1)
	{
		otherVar1 = 0;
		otherVar2 = 2;
	}
	else
	{
		otherVar1 = 0;
		otherVar2 = 1;
	}

	int diff = 0;
	int i;
	for (i = 0; i < Fdegs.getn(); i++)
		if (Fdegs[i][var] == 0)
			diff++;

	Fpcoeff.allocate(Fcoeff.getn() - diff);
	Fpdegs.allocate(Fdegs.getn() - diff);

	int j = 0;
	for (i = 0; i < Fcoeff.getn(); i++)
	{
		if (Fdegs[i][var] != 0)
		{
			Fpcoeff[j] = Fcoeff[i] * Fdegs[i][var];
			Fpdegs[j][var] = Fdegs[i][var] - 1;
			Fpdegs[j][otherVar1] = Fdegs[i][otherVar1];
			Fpdegs[j][otherVar2] = Fdegs[i][otherVar2];
			j++;
		}
	}
}


void mulPolyToPoly(FloatArr xcoeff, FloatArr ycoeff, FloatArr &mulcoeff, V2iArr &muldegs)
{
	mulcoeff.allocate(xcoeff.getn()*ycoeff.getn());
	muldegs.allocate(xcoeff.getn()*ycoeff.getn());
	int place = 0;
	for(int i = 0; i < xcoeff.getn(); i++)
	{
		for(int j = 0; j < ycoeff.getn(); j++)
		{
			bool found = false;
			for(int k = 0; k < place; k++)
			{
				if(muldegs[k][0] == i && muldegs[k][1] == j)
				{
					mulcoeff[k] += xcoeff[place] * ycoeff[place];
					found = true;
					break;
				}
			}
			if(!found)
			{
				mulcoeff[place] = xcoeff[i] * ycoeff[j];
				muldegs[place][0] = i;
				muldegs[place][1] = j;
				place++;
			}
		}
	}
	mulcoeff.shrink(place);
	muldegs.shrink(place);
}