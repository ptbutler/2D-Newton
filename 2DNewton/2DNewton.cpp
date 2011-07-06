/*
* 2dnewton.cpp: Runs 2D Newton method
* Author: Patrick Butler
* Created: 02/10/2010
* Last Modified: 8/19/2010
* 
*/

#include "2dnewton.h"

/************************************************************************
* Takes a polynomial f(x,y) in 2D space and returns a homogenized 
* polynomial F(x,y,w) in projective 2D space
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees for f(x,y)
* <-Fcoeff: the coefficients of F(x,y,w)
* <-Fdegs: the degrees for F(x,y,w)
************************************************************************/

void homogenize2D (FloatArr fcoeff, V2iArr fdegs, FloatArr &Fcoeff, V3iArr &Fdegs)
{
	Fcoeff = fcoeff;
	int highest = 0;

	for(int i = 0; i < fdegs.getn(); i++)
		if(fdegs[i][0] + fdegs[i][1] > highest)
			highest = fdegs[i][0] + fdegs[i][1];

	Fdegs.allocate(fdegs.getn());

	for(int i = 0; i < Fdegs.getn(); i++)
	{
		Fdegs[i][0] = fdegs[i][0];
		Fdegs[i][1] = fdegs[i][1];
		Fdegs[i][2] = highest - (fdegs[i][0] + fdegs[i][1]);
	}
}


/************************************************************************
* Finds the row n of Pascal's triangle in order to get all the binomial
* coefficients of a binomial raised to the nth power
* 
* ->n: the row of Pascal's triangle that is needed
* <-coeff: the binomial coefficients
************************************************************************/

void BinomCoeff (int n, IntArr &bcoeff)
{
	bcoeff.allocate(n + 1);
	bcoeff[0] = 1;

	for (int i = 1; i <= n; i++)
	{
		bcoeff[i] = 1;
		for (int j = i - 1; j > 0; j--)
			bcoeff[j] += bcoeff[j - 1];
	}
}


/************************************************************************
* Finds a polynomial g(x) produced from a binomial raised to a
* a power
* 
* ->a: the first coefficient of the binomial
* ->b: the second coefficient of the binomial
* ->power: the power the binomial is raised to
* <-coeff: the coefficients of g(x)
*
* Assumption: all coefficients are of the terms in increasing sequential
* order
************************************************************************/

void polyToPower (float a, float b, int power, FloatArr &fcoeff)
{
	int size = power + 1;
	fcoeff.allocate(size);

	FloatArr apower, bpower;
	apower.allocate(size);
	bpower.allocate(size);

	apower[0] = 1;
	bpower[0] = 1;

	// Set up arrays of all the terms used to find the polynomial
	int i;
	for (i = 1; i < size; i++)
	{
		apower[i] = apower[i-1]*a;
		bpower[i] = bpower[i-1]*b;
	}
	IntArr bcoeff;
	BinomCoeff(power, bcoeff);

	// Actually multiply all the terms together to get the polynomial
	for (i = 0; i < size; i++)
		fcoeff[i] = apower[i] * bpower[power - i] * bcoeff[i];
}


/************************************************************************
* Reduces an algebraic curve f(x,y) to a univariate polynomial g(t)
* using a parametric line L(t).
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->dir: direction of the line
* ->point: a point on the line
* <-gcoeff: the univariate curve
************************************************************************/

void reduceCurve (FloatArr fcoeff, V2iArr fdegs, V2f dir, V2f point, FloatArr &gcoeff)
{
	gcoeff.allocate(0);
	int i;
	for (i = 0; i < fcoeff.getn(); i++)
	{
		FloatArr hcoeff;

		if (fdegs[i][0] > 0 && fdegs[i][1] > 0)
		{
			FloatArr xcoeff, ycoeff;
			polyToPower(dir[0], point[0], fdegs[i][0], xcoeff);
			polyToPower(dir[1], point[1], fdegs[i][1], ycoeff);
			mulPoly(xcoeff, ycoeff, hcoeff);
		}
		else if (fdegs[i][0] > 0)
			polyToPower(dir[0], point[0], fdegs[i][0], hcoeff);
		else if (fdegs[i][1] > 0)
			polyToPower(dir[1], point[1], fdegs[i][1], hcoeff);
		else
		{
			hcoeff.allocate(1);
			hcoeff[0] = 1;
		}

		// get the coefficient multiplied out
		mulPolyByCoeff(hcoeff, fcoeff[i]);

		// add poly to poly    
		FloatArr tmpcoeff;
		addPoly(hcoeff, gcoeff, tmpcoeff);
		gcoeff = tmpcoeff;
	}

	for (i = gcoeff.getn(); i > 0; i--)
	{
		if(gcoeff[i-1] != 0)
			break;

		gcoeff.shrink(gcoeff.getn()-1);
	}
}


/************************************************************************
* Finds the tangent to an algebraic curve f(x,y) at a given point (a,b)
*  
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->sp: the point (a,b)
* <-tan: the tangent at f(a,b)
************************************************************************/

void findTangent (FloatArr fcoeff, V2iArr fdegs, V2f sp, V2f &tan)
{
	FloatArr fpxcoeff, fpycoeff;
	V2iArr fpxdegs, fpydegs;

	pDerive2D(fcoeff, fdegs, 0, fpxcoeff, fpxdegs);
	pDerive2D(fcoeff, fdegs, 1, fpycoeff, fpydegs);

	V2f normal;
	normal[0] = evalBi(fpxcoeff, fpxdegs, sp);
	normal[1] = evalBi(fpycoeff, fpydegs, sp);
	tan[0] = -normal[1];
	tan[1] = normal[0];  
}


/************************************************************************
* Iterative method to find the intersection of two algebraic curves 
* f(x,y) & g(x,y) utilizing tangents given a starting point (a,b) on
* f(x,y)
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->gcoeff: the coefficients of g(x,y)
* ->gdegs: the degrees of g(x,y)
* ->startPoint: the point (a,b)
* <-endPoint: the intersection of the two curves 
************************************************************************/

bool curveInter (FloatArr fcoeff, V2iArr fdegs, FloatArr gcoeff, V2iArr gdegs,
	V2f startPoint, float length, V2f &endPoint)
{
	float eps = 0.0001;

	while (true)
	{
		V2f tangent;
		findTangent(fcoeff, fdegs, startPoint, tangent);

		FloatArr hcoeff;
		reduceCurve(gcoeff, gdegs, tangent, startPoint, hcoeff);
		FloatArr p1dcoeff;
		removeMultipleRoots(hcoeff, p1dcoeff);

#ifdef _VERBOSEMODE_		
		// debugging code
		cout<<"p1dcoeff"<<endl;
		for(int i = 0; i < p1dcoeff.getn(); i++)
			cout<<i<<": "<<p1dcoeff[i]<<endl;
		cout<<"startpoint[0]: "<<startPoint[0]<<endl;
		cout<<"startpoint[1]: "<<startPoint[1]<<endl;
#endif

		// pass poly to 1dnewton
		V2f interval;
		interval[0] = -length;
		interval[1] = 1000 * length;
		float guess = 0.01;
		float root = 0;

		bool success;
		success = OneDNewton(p1dcoeff, guess, interval, root);
		if (!success)
		{
			cout<<"***** Error with OneDNewton! *** Curve Intersection *****"<<endl;
#ifdef _VERBOSEMODE_
			for(int i = 0; i < p1dcoeff.getn(); i++)
				cout<<"p1dcoeff["<<i<<"]: "<<p1dcoeff[i]<<endl;
			cout<<"guess: "<<guess<<endl;
			cout<<"interval[0]: "<<interval[0]<<", interval[1]:"<<interval[1]<<endl;
			cout<<"root: "<<root<<endl;
			cout<<"***** End Error Message *********************************"<<endl;
#endif			
			success = OneDNewton(p1dcoeff, interval[1], interval, root);
			if(!success)
			{
				success = OneDNewton(p1dcoeff, (interval[0]+interval[1])/2, interval, root);
				if(!success)
					return false;
			}
				
		}

		endPoint[0] = tangent[0] * root + startPoint[0];
		endPoint[1] = tangent[1] * root + startPoint[1];

		// if endpoint is close enough to being on algebraic curve
		float eval = evalBi(fcoeff, fdegs, endPoint);
		if (fabs(eval) < eps)
			return true;                              // we are done

		startPoint = endPoint;
		FloatArr tcoeff = fcoeff;
		V2iArr tdegs = fdegs;
		fcoeff = gcoeff;
		fdegs = gdegs;
		gcoeff = tcoeff;
		gdegs = tdegs;
	}
}


/************************************************************************
* Compute the Sturm sequence of a polynomial.
* 
* ->fcoeff: the coefficients of the polynomial
* <-sturm: the strum sequence
************************************************************************/

void genSturm (FloatArr fcoeff, FloatArrArr &sturm)
{
	sturm.allocate(2);
	sturm[0] = fcoeff;
	FloatArr fprime;
	derive(fcoeff, fprime);
	sturm[1] = fprime;

	while ((sturm[sturm.getn()-1]).getn() > 1)
	{
		FloatArr q, rem;
		dividePoly(sturm[sturm.getn()-2], sturm[sturm.getn()-1], q, rem);

		int zeros = 0;
		for (int i = 0; i < rem.getn(); i++)
			if (rem[i] == 0)
				zeros++;

		if (zeros == rem.getn())
			break;

		mulPolyByCoeff(rem, -1);
		FloatArrArr tmp;
		tmp.allocate(sturm.getn()+1);
		for (int i = 0; i < sturm.getn(); i++)
			tmp[i] = sturm[i];

		tmp[tmp.getn()-1] = rem;
		sturm = tmp;
	}
}


/************************************************************************
* Evaluate a Sturm sequence at a specified value a, and count the sign changes.
* 
* ->sturm: the Sturm sequence
* ->a: a
* <-: the number of sign changes
************************************************************************/

int evalSturm (FloatArrArr sturm, float a)
{
	int count = 0;
	bool pos;
	bool neg;
	float value;

	FloatArr tmp = sturm[0];
	value = polyEval(tmp, a);

	if (value > 0)
	{
		pos = true;
		neg = false;
	}
	else if (value < 0)
	{
		pos = false;
		neg = true;
	}
	else
	{
		pos = false;
		neg = false;
	}
	for (int i = 1; i < sturm.getn(); i++)
	{
		value = polyEval(sturm[i], a);
		if (value > 0 && neg == true)
		{
			neg = false;
			pos = true;
			count++;
		}
		else if (value < 0 && pos == true)
		{
			neg = true;
			pos = false;
			count++;
		}
	}
	return count;
}


/************************************************************************
* Given a Sturm sequence and bounding interval (a,b)
* find the number of roots contained within
* 
* ->interval: the candidate interval (a,b)
* ->sturm:the Sturm Sequence
* <-: the number of roots in this interval
* 
************************************************************************/

int numOfRoots (V2f interval, FloatArrArr sturm)
{
	return evalSturm(sturm, interval[0]) - evalSturm(sturm, interval[1]);
}


/************************************************************************
* Given a univariate polynomial f(x), remove any multiple roots and 
* create g(x)
* 
* ->fcoeff: the coefficients of f(x)
* <-gcoeff: the coefficients of g(x)
* 
************************************************************************/

void removeMultipleRoots(FloatArr fcoeff, FloatArr &gcoeff)
{
	FloatArr fpx;
	derive(fcoeff, fpx);
	FloatArr gcd;
	polyGCD(fcoeff, fpx, gcd); // find gcd of f and f'
	FloatArr r;
	dividePoly(fcoeff, gcd, gcoeff, r); // g = f / gcd, r should come back 0
}

/************************************************************************
* Given a univariate polynomial f(x) and bounding interval (a,b)
* find the number of roots contained within
* 
* ->fcoeff: the coefficients of f(x)
* ->interval: the candidate interval (a,b)
* <-: the number of roots in this interval
* 
************************************************************************/

int numOfRoots (FloatArr fcoeff, V2f interval)
{
	//FloatArr gcoeff;
	//removeMultipleRoots(fcoeff, gcoeff);
	FloatArrArr sturm;
	genSturm(fcoeff, sturm);
	return evalSturm(sturm, interval[0]) - evalSturm(sturm, interval[1]);
}


/************************************************************************
* Given a univariate polynomial f(x), a bounding interval (a,b), and a
* Sturm sequence, find the bounding intervals for all roots within (a,b).
* 
* ->fcoeff: the coefficients of f(x)
* ->interval: (a,b)
* ->sturm: the Sturm sequence
* <-bd: ith root lies in the interval (bd[i][0], bd[i][1]) 
************************************************************************/

void getRootIntervals (V2f interval, FloatArrArr sturm, V2fArr &bd)
{
	int numRoots = numOfRoots(interval, sturm);
	float eps = 0.0001;
	if (numRoots == 1)
	{
		bd.allocate(1);
		bd[0] = interval;
	}
	else if (numRoots == 0)
		bd.allocate(0);
	else if (fabs(interval[0] - interval[1]) < eps)
	{
#ifdef _VERBOSEMODE_
		cout<<"**** Error: Interval is too small *** Interval finding ****"<<endl;
		cout<<"numRoots: "<<numRoots<<endl;
		cout<<"interval[0]: "<<interval[0]<<endl;
		cout<<"interval[1]: "<<interval[1]<<endl;
		cout<<"**** End Error message ******************************************"<<endl;
#endif
		bd.allocate(0);
	}
	else if (numRoots > 1)
	{
		float half = (interval[1] + interval[0]) / 2;
		float val = polyEval(sturm[0], half);
		if(fabs(val) < eps)
			half += eps;

		V2fArr half1, half2;
		V2f int1, int2;
		int1[0] = interval[0];
		int1[1] = half;
		int2[0] = half;
		int2[1] = interval[1];
		getRootIntervals(int1, sturm, half1);
		getRootIntervals(int2, sturm, half2);

		bd.allocate(half1.getn() + half2.getn());
		int i;
		for(i = 0; i < half1.getn(); i++)
			bd[i] = half1[i];
		int j;
		for(j = 0; j < half2.getn(); j++)
			bd[i+j] = half2[j];
	}
	else
	{
		cout<<"**** Error: Sturm sequence value of < 0 *** Interval finding ****"<<endl;
#ifdef _VERBOSEMODE_
		cout<<"numRoots: "<<numRoots<<endl;
		cout<<"interval[0]: "<<interval[0]<<endl;
		cout<<"interval[1]: "<<interval[1]<<endl;
		cout<<"**** End Error message ******************************************"<<endl;
#endif
		bd.allocate(0);
	}
}


/************************************************************************
* Given a univariate polynomial f(x) and a bounding interval (a,b),
* find the bounding intervals for all roots within (a,b).
* 
* ->fcoeff: the coefficients of f(x)
* ->interval: (a,b)
* <-bd: ith root lies in the interval (bd[i][0], bd[i][1]) 
************************************************************************/

void getRootIntervals (FloatArr fcoeff, V2f interval, V2fArr &bd)
{
	//FloatArr gcoeff;
	//removeMultipleRoots(fcoeff, gcoeff);
	FloatArrArr sturm;
	genSturm(fcoeff, sturm);
	getRootIntervals(interval, sturm, bd);
}


/************************************************************************
* Finds a point (a,b) on the algebraic curve f(x,y) to start at
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->guess: the guess of a starting point
* <-start: the point (a,b)
************************************************************************/

void getStartPoint (FloatArr fcoeff, V2iArr fdegs, V2f guess, V2f &start)
{
	V2f dir;
	dir[0] = 1;
	dir[1] = 0;
	dir.normalize();

	FloatArr p1dcoeff;
	reduceCurve(fcoeff, fdegs, dir, guess, p1dcoeff);

	V2f maxInt;
	maxInt[0] = -10000; //change to what min could be
	maxInt[1] = 10000;  //change to what max could be

	V2fArr intervals;
	getRootIntervals(p1dcoeff, maxInt, intervals);

	float half = (intervals[0][0] + intervals[0][1]) / 2;
	float root;
	V2f newtonInt;
	newtonInt[0] = intervals[0][0];
	newtonInt[1] = intervals[0][1];

	bool success;
	while (true)
	{
		success = OneDNewton(p1dcoeff, half, newtonInt, root);
		if (success)
			break;

		half = (newtonInt[0] + newtonInt[1]) / 2;
		V2f tmpInt;
		tmpInt[0] = newtonInt[0];
		tmpInt[1] = half;

		if (numOfRoots(p1dcoeff, tmpInt) > 0)
			newtonInt[1] = half;
		else
			newtonInt[0] = half;
	}
	start[0] = dir[0] * root + guess[0];
	start[1] = dir[1] * root + guess[1];
}


/************************************************************************
* Finds the maximum bounds of the real roots of the curve p1d
* 
* ->p1dcoeff: the coefficients of p1d(x,y)
* <-maxInt: the bounds
************************************************************************/

void findMaxInt(FloatArr p1dcoeff, V2f maxInt)
{
	float max = fabs(p1dcoeff[1]);
	for(int i = 2; i < p1dcoeff.getn()-1; i++)
		if(fabs(p1dcoeff[i]) > max)
			max = fabs(p1dcoeff[i]);

	max = 1 + max / fabs(p1dcoeff[p1dcoeff.getn()-1]);
	maxInt[0] = -max;
	maxInt[1] = max;
}


/************************************************************************
* Finds points on the algebraic curve f(x,y) to start at
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->guess: the guess of a starting point
* <-specials: the points
************************************************************************/

void getSpecialPoints (FloatArr fcoeff, V2iArr fdegs, V2f guess, V2f dir, V2fArr &specials)
{
	//float eps = 0.1;
	int i;
	dir.normalize();

	FloatArr p1dcoeff;
	reduceCurve(fcoeff, fdegs, dir, guess, p1dcoeff);
	FloatArr gcoeff;
	removeMultipleRoots(p1dcoeff, gcoeff);

	V2f maxInt;
	maxInt[0] = -10000; //change to what min could be
	maxInt[1] = 12500;  //change to what max could be
	findMaxInt(p1dcoeff, maxInt);

#ifdef _VERBOSEMODE_
	cout<<"p1dcoeff to find starting points"<<endl;
	for(int i = 0; i < p1dcoeff.getn(); i++) {
		cout<<i<<": ("<<p1dcoeff[i]<<")"<<endl;
	}
	cout<<"gcoeff, multiple roots removed(?)"<<endl;
	for(int i = 0; i < gcoeff.getn(); i++) {
		cout<<i<<": ("<<gcoeff[i]<<")"<<endl;
	}
	cout<<"interval: ("<<maxInt[0]<<", "<<maxInt[1]<<")"<<endl;
#endif

	V2fArr intervals;
	//getRootIntervals(p1dcoeff, maxInt, intervals);
	getRootIntervals(gcoeff, maxInt, intervals);
	specials.allocate(intervals.getn());

	for (i = 0; i < intervals.getn(); i++)
	{
		float half = (intervals[i][0] + intervals[i][1]) / 2;
		float root;
		V2f newtonInt;
		newtonInt[0] = intervals[i][0];
		newtonInt[1] = intervals[i][1];

		bool success;
		V2f point;
		int infloop = 1000;
		int count = 0;
		while (true && newtonInt[0] < newtonInt[1] && count < infloop)
		{
			//success = OneDNewton(p1dcoeff, half, newtonInt, root);
			success = OneDNewton(gcoeff, half, newtonInt, root);
			point[0] = dir[0] * root + guess[0];
			point[1] = dir[1] * root + guess[1];
			if (success) //&& fabs(evalBi(fcoeff, fdegs, point)) < eps)
				break;
			else
			{
				//success = OneDNewton(p1dcoeff, newtonInt[0], newtonInt, root);
				success = OneDNewton(gcoeff, newtonInt[0], newtonInt, root);
				point[0] = dir[0] * root + guess[0];
				point[1] = dir[1] * root + guess[1];
				if (success)
					break;
				else
				{
					//success = OneDNewton(p1dcoeff, newtonInt[1], newtonInt, root);
					success = OneDNewton(gcoeff, newtonInt[1], newtonInt, root);
					point[0] = dir[0] * root + guess[0];
					point[1] = dir[1] * root + guess[1];
					if (success)
						break;
				}
			}
			

			half = (newtonInt[0] + newtonInt[1]) / 2;
			V2f tmpInt;
			tmpInt[0] = newtonInt[0];
			tmpInt[1] = half;

			if (numOfRoots(p1dcoeff, tmpInt) > 0)
				newtonInt[1] = half;
			else
				newtonInt[0] = half;
			count++;
		}
		specials[i] = point;
	}
}


/************************************************************************
* Finds the position of a singularity close to a point on the algebraic 
* curve f. This is found using resultants.
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->startPoint: point close to singularity
* ->singularity: the singularity on f
* <-hcoeff: the coefficients of h(x,y)
* <-hdegs: the degrees of h(x,y)
************************************************************************/

void findCloseSingularityResultant (FloatArr fcoeff, V2iArr fdegs, V2f startPoint, V2f &singularity)
{
	// Get the partial derivatives of f
	FloatArr fxcoeff, fycoeff;
	V2iArr fxdegs, fydegs;
	pDerive2D(fcoeff, fdegs, 0, fxcoeff, fxdegs);
	pDerive2D(fcoeff, fdegs, 1, fycoeff, fydegs);

	// Build the resultant with the partials


	// Find the root of resultant close to startPoint
}


/************************************************************************
* Finds the position of a singularity close to a point on the algebraic 
* curve f. This is found using Newton's curve intersection.
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->startPoint: point close to singularity
* ->singularity: the singularity on f
* <-hcoeff: the coefficients of h(x,y)
* <-hdegs: the degrees of h(x,y)
************************************************************************/

bool findCloseSingularityIntersect (FloatArr fcoeff, V2iArr fdegs, V2f startPoint, V2f &singularity)
{
	// Get the partial derivatives of f
	FloatArr fxcoeff, fycoeff;
	V2iArr fxdegs, fydegs;
	pDerive2D(fcoeff, fdegs, 0, fxcoeff, fxdegs);
	pDerive2D(fcoeff, fdegs, 1, fycoeff, fydegs);

	int length = 1;
	// Use these to find the intersection
	bool found = curveInter(fxcoeff, fxdegs, fycoeff, fydegs, startPoint, length, singularity);
	if(found)
		cout<<"Singularity found at ("<<singularity[0]<<", "<<singularity[1]<<")"<<endl;
	return found;
}


/************************************************************************
* Moves a known singularity of an algebraic curve f to the origin by 
* shifting the curve to make a new curve h
* 
* ->fcoeff: the coefficients of f(x,y)
* ->fdegs: the degrees of f(x,y)
* ->singularity: the singularity on f
* <-hcoeff: the coefficients of h(x,y)
* <-hdegs: the degrees of h(x,y)
************************************************************************/

void singularityToOrigin(FloatArr fcoeff, V2iArr fdegs, V2f singularity, FloatArr &hcoeff, V2iArr &hdegs)
{
	hcoeff.allocate(fcoeff.getn());
	hdegs.allocate(fcoeff.getn());
	int place = 0;
	// for each term in f, expand it, then add it to h
	for(int i = 0; i < fcoeff.getn(); i++)
	{
		FloatArr xcoeff, ycoeff, mulcoeff;
		V2iArr muldegs;
		IntArr xbincoeff, ybincoeff;

		// expand this term in x
		BinomCoeff(fdegs[i][0], xbincoeff);
		xcoeff.allocate(xbincoeff.getn());
		for(int j = 0; j < xbincoeff.getn(); j++)
		{
			xcoeff[j] = xbincoeff[j] * pow(singularity[0], xbincoeff.getn()-(j+1));
			float temp = xcoeff[j];
		}

		// expand this term in y
		BinomCoeff(fdegs[i][1], ybincoeff);
		ycoeff.allocate(ybincoeff.getn());
		for(int j = 0; j < ybincoeff.getn(); j++)
		{
			ycoeff[j] = ybincoeff[j] * pow(singularity[1], ybincoeff.getn()-(j+1));
			float temp = ycoeff[j];
			float temp2 = singularity[1];
		}

		// multiply these together
		mulPolyToPoly(xcoeff, ycoeff, mulcoeff, muldegs);

		// multiply by the original coefficient
		for(int j = 0; j < mulcoeff.getn(); j++)
		{
			mulcoeff[j] *= fcoeff[i];
			float tmp = mulcoeff[j];
		}

		// add to new curve h
		for(int j = 0; j < mulcoeff.getn(); j++)
		{
			bool found = false;
			for(int k = 0; k < place; k++)
			{
				if(muldegs[j][0] == hdegs[k][0] && muldegs[j][1] == hdegs[k][1])
				{
					hcoeff[k] += mulcoeff[j];
					found = true;
					break;
				}
			}
			if(!found)
			{
				hcoeff[place] = mulcoeff[j];
				hdegs[place] = muldegs[j];
				place++;
				if(place == hcoeff.getn())
				{
					hcoeff.grow(10);
					hdegs.grow(10);
				}
			}
		}
#ifdef _VERBOSEMODE_
		// debugging code
		cout<<"f["<<i<<"]: "<<fcoeff[i]<<" "<<fdegs[i][0]<<" "<<fdegs[i][1]<<endl;
		cout<<"h:"<<endl;
		for(int i = 0; i < place; i++)
			cout<<i<<": "<<hcoeff[i]<<" "<<hdegs[i][0]<<" "<<hdegs[i][1]<<endl;
#endif
	}
	hcoeff.shrink(place);
	hdegs.shrink(place);

	// Do I need to remove any zero entries, and possibly and entries that are very close to zero?
	float max = 0;
	for(int i = 0; i < hcoeff.getn(); i++)
		if(fabs(hcoeff[i]) > max)
			max = fabs(hcoeff[i]);

	for(int i = 0; i < hcoeff.getn(); i++)
		if(fabs(hcoeff[i] * 1000000) < max)
			hcoeff[i] = 0;
	int i;
	for(i = 0; i < hcoeff.getn(); i++)
	{
		if(hcoeff[i] == 0)
		{
			bool replaced = false;
			for(int j = hcoeff.getn()-1; j > i; j--)
			{
				if(hcoeff[j] != 0)
				{
					hcoeff[i] = hcoeff[j];
					hcoeff[j] = 0;
					hdegs[i] = hdegs[j];
					replaced = true;
					break;
				}
			}
			if(!replaced)
				break;
		}
	}
	
	hcoeff.shrink(i);
	hdegs.shrink(i);
#ifdef _VERBOSEMODE_
	// debugging code
	cout<<"h:"<<endl;
	for(int i = 0; i < hcoeff.getn(); i++)
		cout<<i<<": "<<hcoeff[i]<<" "<<hdegs[i][0]<<" "<<hdegs[i][1]<<endl;
#endif
}


/************************************************************************
* Moves the point close to the singularity toward the origin the same 
* distance the singularity moved
* 
* ->start: point to be shifted
* ->singularity: the singularity on f
* <-hstart: the shift point
************************************************************************/

void shiftStart(V2f start, V2f singularity, V2f &hstart)
{
	hstart[0] = start[0] + singularity[0];
	hstart[1] = start[1] + singularity[1];
}


/************************************************************************
* Moves the points from the shifted curve h back to the original curve 
* f
* 
* ->hpoints: points to be shifted
* ->singularity: the singularity on f
* <-points: the shifted points
************************************************************************/

void inverseToOrigin(V2fArr hpoints, V2f singularity, V2fArr &points)
{
	points.allocate(hpoints.getn());
	for(int i = 0; i < points.getn(); i++)
	{
		points[i][0] = hpoints[i][0] + singularity[0];
		points[i][1] = hpoints[i][1] + singularity[1];
	}
}


/************************************************************************
 * Preforms a quadratic transformation on the algebraic curve f to 
 * generate the new algebraic curve g. The transformation used is:
 *    x = x1    y = x1 * y1
 * 
 * ->fcoeff: the coefficients of f(x,y)
 * ->fdegs: the degrees of f(x,y)
 * <-gcoeff: the coefficients of g(x,y)
 * <-gdegs: the degrees of g(x,y)
 ************************************************************************/

void quadraticTransform (FloatArr fcoeff, V2iArr fdegs, FloatArr &gcoeff, V2iArr &gdegs)
{
	gcoeff.allocate(fcoeff.getn());
	gdegs.allocate(fdegs.getn());

	for(int i = 0; i < gcoeff.getn(); i++)
	{
		gcoeff[i] = fcoeff[i];
		gdegs[i][0] = fdegs[i][0] + fdegs[i][1];
		gdegs[i][1] = fdegs[i][1];
	}
}


/************************************************************************
 * Preforms a quadratic transformation on the point start to generate the 
 * new point gstart. The transformation used is:
 *    x = x1    y = x1 * y1
 * 
 * ->start: the original point
 * <-gstart: the transformed point
 ************************************************************************/

void quadraticTransform (V2f start, V2f &gstart)
{
	gstart[0] = start[0];
	if(gstart[0] != 0)
		gstart[1] = start[1] / gstart[0];
	else if(start[1] == 0)
		gstart[1] = 0;
	else if(start[1] > 0)
		gstart[1] = FLT_MAX;
	else
		gstart[1] = -FLT_MAX;
}


/************************************************************************
 * Preforms an inverse of the quadratic transformation on the points 
 * gpoints to generate the new points fpoints in the original space. The 
 * transformation being reveresed is:
 *    x = x1    y = x1 * y1
 * 
 * ->gpoints: the original point
 * <-fpoints: the transformed point
 ************************************************************************/

void inverseQuadraticTransform(V2fArr gpoints, V2fArr &fpoints)
{
	fpoints.allocate(gpoints.getn());

	for(int i = 0; i < gpoints.getn(); i++)
	{
		fpoints[i][0] = gpoints[i][0];
		fpoints[i][1] = gpoints[i][0] * gpoints[i][1];
	}
}


/************************************************************************
 * Finds the points from gstart to get through the singularity in the
 * the original curve while on the transformed curve g by going past the 
 * x-axis
 * 
 * ->gcoeff: the coefficients of g(x,y)
 * ->gdegs: the degrees of g(x,y)
 * ->gstart: the point before the x-axis
 * ->length: distance to jump between each point
 * <-gpoints: the points moving past the x-axis
 ************************************************************************/

void pointsThroughSingularity (FloatArr gcoeff, V2iArr gdegs, V2f gstart, float length, V2fArr &gpoints)
{
	V2f tangent, oldtangent, gprevious, glast, endPoint;
	// need to continue as before

	int j = 1;
	float epsilon = 0.001;
	if(gstart[0] < 0)
	{
		oldtangent[0] = 1;
		oldtangent[1] = 0;
	}
	else
	{
		oldtangent[0] = -1;
		oldtangent[1] = 0;
	}
	gprevious = gstart;
	glast = gstart;
	while(true)
	{
		findNextPointSing(gcoeff, gdegs, gprevious, glast, length, j, oldtangent, tangent, endPoint);
		float tempa = endPoint[0];
		float tempb = endPoint[1];
		gpoints.append(endPoint);
		//oldtangent = tangent;
		gprevious = glast;
		glast = endPoint;

		// get far enough past the line x = 0
		if(endPoint[0] > (0 + epsilon) && gstart[0] < 0)
			break;
		else if(endPoint[0] < (0 + epsilon) && gstart[0] > 0)
			break;
	}
}


/************************************************************************
 * Handles moving through a singularity by transforming the curve f, 
 * moving past the singularity, and then transforming back
 * 
 * ->fcoeff: the coefficients of f(x,y)
 * ->fdegs: the degrees of f(x,y)
 * ->start: point before the singularity
 * ->singularity: the singularitys
 * ->length: distance to jump between each point 
			 (passed on to pointsThroughSingularity)
 * <-points: the points moving past the singularity
 ************************************************************************/

void jumpThroughSingularity (FloatArr fcoeff, V2iArr fdegs, V2f start, V2f singularity, float length, V2fArr &points)
{
	// shift singularity to origin
	FloatArr hcoeff;
	V2iArr hdegs;
	V2f hstart;
	singularityToOrigin(fcoeff, fdegs, singularity, hcoeff, hdegs);
	shiftStart(start, singularity, hstart);

	// perform quadratic transformation x = xi, y = xi*yi
	FloatArr gcoeff;
	V2iArr gdegs;
	V2f gstart;
	quadraticTransform(hcoeff, hdegs, gcoeff, gdegs);
	quadraticTransform(hstart, gstart);

	// move past representation of singularity
	V2fArr gpoints;
	pointsThroughSingularity(gcoeff, gdegs, gstart, length, gpoints);

	// convert these points back using inverse of quadratic transformation
	V2fArr hpoints;
	inverseQuadraticTransform(gpoints, hpoints);

	// convert these points back to original using inverse of move
	inverseToOrigin(hpoints, singularity, points);
}


/************************************************************************
 * Checks to see if a singularity on  the curve f is close, or has been 
 * jumped. If a singularity is detected, its location is determined.
 * 
 * ->fcoeff: the coefficients of f(x,y)
 * ->fdegs: the degrees of f(x,y)
 * ->lastPoint: last point, before the singularity
 * ->thisPoint: current point, unknown compared to singularity
 * ->tangent: direction we are heading
 * <-singularity: the singularity if found
 * <-: if we are near or just past a singularity true is returned
 *
 * Assumption: After testing, +/- 0.00001 is the tolerance of the 
 * quadratic transformation. The singularity must be located within 
 * these bounds.
 ************************************************************************/

bool isNearOrPastSingularity(FloatArr fcoeff, V2iArr fdegs, V2f lastPoint, V2f thisPoint, V2f tangent, bool &past)
{
	float singError = 0.1;
	if (fabs(tangent[0]) < singError && fabs(tangent[1]) < singError)
		return true;
	else if (fabs(tangent[0]) < singError*10 && fabs(tangent[1]) < singError*10)
	{
		V2f lastTangent;
		findTangent(fcoeff, fdegs, lastPoint, lastTangent);
		bool xtrue = false, ytrue = false;
		if((lastTangent[0] > 0 && tangent[0] < 0) || (lastTangent[0] < 0 && tangent[0] > 0))
		{
			xtrue = true;
			past = true;
		}
		else
		{
			FloatArr fpxcoeff;
			V2iArr fpxdegs;
			pDerive2D(fcoeff, fdegs, 0, fpxcoeff, fpxdegs);

			FloatArr fp2xcoeff;
			V2iArr fp2xdegs;
			pDerive2D(fcoeff, fdegs, 0, fp2xcoeff, fp2xdegs);

			if(lastTangent[0] > 0 && evalBi(fp2xcoeff, fp2xdegs, lastPoint) < 0 && tangent[0] > 0 && evalBi(fp2xcoeff, fp2xdegs, thisPoint) > 0)
			{
				xtrue = true;
				past = true;
			}
			else if(lastTangent[0] < 0 && evalBi(fp2xcoeff, fp2xdegs, lastPoint) > 0 && tangent[0] < 0 && evalBi(fp2xcoeff, fp2xdegs, thisPoint) < 0)
			{
				xtrue = true;
				past = true;
			}
		}

		if((lastTangent[1] > 0 && tangent[1] < 0) || (lastTangent[1] < 0 && tangent[1] > 0))
		{
			ytrue = true;
			past = true;
		}
		else
		{
			FloatArr fpycoeff;
			V2iArr fpydegs;
			pDerive2D(fcoeff, fdegs, 1, fpycoeff, fpydegs);

			FloatArr fp2ycoeff;
			V2iArr fp2ydegs;
			pDerive2D(fcoeff, fdegs, 1, fp2ycoeff, fp2ydegs);

			if(lastTangent[1] > 0 && evalBi(fp2ycoeff, fp2ydegs, lastPoint) < 0 && tangent[1] > 0 && evalBi(fp2ycoeff, fp2ydegs, thisPoint) > 0)
			{
				ytrue = true;
				past = true;
			}
			else if(lastTangent[1] < 0 && evalBi(fp2ycoeff, fp2ydegs, lastPoint) > 0 && tangent[1] < 0 && evalBi(fp2ycoeff, fp2ydegs, thisPoint) < 0)
			{
				ytrue = true;
				past = true;
			}
		}
		return xtrue && ytrue;
	}
	return false;
}


void findNextPointSing(FloatArr fcoeff, V2iArr fdegs, V2f lastPoint, V2f startPoint, float length, int i, V2f &oldtangent, V2f &tangent, V2f &endPoint)
{
	FloatArr circoeff;
	V2iArr cirdegs;
	V2f point;

	circoeff.allocate(5);
	cirdegs.allocate(5);
	cirdegs[0][0] = 2;
	cirdegs[0][1] = 0;
	cirdegs[1][0] = 1;
	cirdegs[1][1] = 0;
	cirdegs[2][0] = 0;
	cirdegs[2][1] = 1;
	cirdegs[3][0] = 0;
	cirdegs[3][1] = 2;
	cirdegs[4][0] = 0;
	cirdegs[4][1] = 0;

	circoeff[0] = 1;
	circoeff[1] = -2 * startPoint[0];
	circoeff[2] = -2 * startPoint[1];
	circoeff[3] = 1;
	circoeff[4] = startPoint[0] * startPoint[0] + startPoint[1] * startPoint[1] - (length * length);   

	findTangent(fcoeff, fdegs, startPoint, tangent);	
	tangent.normalize();

	if (i == 0)
	{
		if((startPoint[0] > 0 && tangent[0] > 0) || (startPoint[0] < 0 && tangent[0] < 0)) {
			tangent[0] = -tangent[0];
			tangent[1] = -tangent[1];
		}
		oldtangent = tangent;
	}
	if (tangent.dot(oldtangent) < 0)
	{
		tangent[0] = -tangent[0];
		tangent[1] = -tangent[1];
	}
	

	point[0] = tangent[0] * length + startPoint[0];
	point[1] = tangent[1] * length + startPoint[1];

	curveInter(circoeff, cirdegs, fcoeff, fdegs, point, length, endPoint);
}


void findNextPoint(FloatArr fcoeff, V2iArr fdegs, V2f lastPoint, V2f startPoint, float length, int i, V2f &oldtangent, V2f &tangent, V2f &endPoint)
{
	FloatArr circoeff;
	V2iArr cirdegs;
	V2f point;

	circoeff.allocate(5);
	cirdegs.allocate(5);
	cirdegs[0][0] = 2;
	cirdegs[0][1] = 0;
	cirdegs[1][0] = 1;
	cirdegs[1][1] = 0;
	cirdegs[2][0] = 0;
	cirdegs[2][1] = 1;
	cirdegs[3][0] = 0;
	cirdegs[3][1] = 2;
	cirdegs[4][0] = 0;
	cirdegs[4][1] = 0;

	circoeff[0] = 1;
	circoeff[1] = -2 * startPoint[0];
	circoeff[2] = -2 * startPoint[1];
	circoeff[3] = 1;
	circoeff[4] = startPoint[0] * startPoint[0] + startPoint[1] * startPoint[1] - (length * length);   

	findTangent(fcoeff, fdegs, startPoint, tangent);

	// need this to not necessarily be arbitrary
	//float singError = 0.1;
	// fabs(tangent[0]) < singError && fabs(tangent[1]) < singError
	bool past = false;
	if (isNearOrPastSingularity(fcoeff, fdegs, lastPoint, startPoint, tangent, past))
	//if (fabs(tangent[0]) < singError && fabs(tangent[1]) < singError)
	{
		// get past the singularity
		V2fArr points;
		V2f singularity;

#ifdef _VERBOSEMODE_
		cout<<"Singularity nearby: ("<<startPoint[0]<<", "<<startPoint[1]<<")"<<endl;
#endif
		bool singFound = false;
		if(past)
			singFound = findCloseSingularityIntersect (fcoeff, fdegs, lastPoint, singularity);
		else
			singFound = findCloseSingularityIntersect (fcoeff, fdegs, startPoint, singularity);

		if(singFound)
		{
			if(past)
				jumpThroughSingularity (fcoeff, fdegs, lastPoint, singularity, length, points);
			else
				jumpThroughSingularity (fcoeff, fdegs, startPoint, singularity, length, points);

			// return the next point
			endPoint = points[points.getn()-1];
			if(past)
			{
				tangent[0] = endPoint[0] - lastPoint[0];
				tangent[1] = endPoint[1] - lastPoint[1];
			}
			else
			{
				tangent[0] = endPoint[0] - startPoint[0];
				tangent[1] = endPoint[1] - startPoint[1];
			}
			tangent.normalize();
			return;
		}
	}

	tangent.normalize();
	if (i == 0)
		oldtangent = tangent;
	if (tangent.dot(oldtangent) < 0)
	{
		tangent[0] = -tangent[0];
		tangent[1] = -tangent[1];
	}

	point[0] = tangent[0] * length + startPoint[0];
	point[1] = tangent[1] * length + startPoint[1];

	curveInter(circoeff, cirdegs, fcoeff, fdegs, point, length, endPoint);
}


/************************************************************************
* Iterative method to find the next point along a bivariate algebraic 
* curve f(x,y)=0 using tangents given a guessing point (a,b) near the 
* curve and the length to move along the curve
* 
* ->fcoeff: the coefficients of f(x,y)=0
* ->fdegs: the degrees of f(x,y)=0
* ->guessPoint: (a,b)
* ->length: length along the tangent
* ->numPoint: number of points to find along the curve
* <-endPoint: the next point along the curve
* <-: returns the number of successfully found points
************************************************************************/

int TwoDNewton (FloatArr fcoeff, V2iArr fdegs, V2f guessPoint,
	float length, int numPoints, V2fArrArr &endPoints)
{
	int successes = 0;
	int totalcount = 0;

	V2f startPoint;
	V2fArr specials;
	V2f tangent;
	V2f endPoint;
	int i;

	// erroring out here, need to investigate
	V2f dir;
	dir[0] = 1;
	dir[1] = 0;
	//getSpecialPoints(fcoeff, fdegs, guessPoint, dir, specials);
	// just hard coding the points to start at
	specials.allocate(3);
	specials[0][0] = 0;
	specials[0][1] = -1;
	specials[1][0] = 0;
	specials[1][1] = 1;
	specials[2][0] = 1;
	specials[2][1] = 0.8104655;

#ifdef _VERBOSEMODE_
	for(i = 0; i < specials.getn(); i++)
		cout<<"specials["<<i<<"][0]: "<<specials[i][0]<<endl
			<<"specials["<<i<<"][1]: "<<specials[i][1]<<endl;
#endif

	endPoints.allocate(specials.getn());
	list<int> unvisited;
	for(i = 0; i < specials.getn(); i++)
		unvisited.push_back(i);

	i = 0;
	endPoints[i].allocate(10);
	while(unvisited.size() > 0)
	{
		int current = unvisited.front();
		startPoint = specials[current];
		if(endPoints[i].getn() == successes)
			endPoints[i].grow(10);
		endPoints[i][successes] = startPoint;
		successes++;
		V2f firsttangent;
		V2f oldtangent;
		V2f lastPoint = startPoint;
		V2f thisPoint = startPoint;

		bool flipped = false;
		int count = 0;
		while(true)
		{
			if(flipped && count == 0)
				findNextPoint(fcoeff, fdegs, lastPoint, thisPoint, length, 1, oldtangent, tangent, endPoint);
			else
				findNextPoint(fcoeff, fdegs, lastPoint, thisPoint, length, count, oldtangent, tangent, endPoint);

			if(count == 0)
				firsttangent = tangent;
			if(endPoints[i].getn() == successes)
				endPoints[i].grow(10);
			endPoints[i][successes] = endPoint;
			lastPoint = thisPoint;
			thisPoint = endPoint;
			oldtangent = tangent;
			count++;
			successes++;

			list<int>::iterator iter;
			int prevSize = unvisited.size();
			bool quit = false;
			for(iter = unvisited.begin(); iter != unvisited.end(); iter++)
			{
				if(*iter == current && count < 2)
					continue;  
				else if(*iter == current && count > 200 && !flipped)
				{
					flipped = true;
					lastPoint = startPoint;
					thisPoint = startPoint;
					oldtangent[0] = -firsttangent[0];
					oldtangent[1] = -firsttangent[1];
					count = 0;
					break;
				}
				else if(*iter == current && count > 200 && flipped)
				{
					unvisited.erase(iter);
					if(endPoints[i].getn() > successes)
						endPoints[i].shrink(successes);
					i++;
					totalcount += successes;
					successes = 0;
					quit = true;
					break;
				}
				if(endPoint.dist(specials[*iter]) < length)
				{
					if(*iter == current)
					{
						unvisited.erase(iter);
						if(endPoints[i].getn() > successes)
							endPoints[i].shrink(successes);
						i++;
						totalcount += successes;
						successes = 0;
						quit = true;
						break;
					}
					list<int>::iterator iter2;
					iter2 = iter;
					iter--;
					unvisited.erase(iter2);
					endPoints.shrink(endPoints.getn()-1);
				}
			}
			int newSize = unvisited.size();
			if(newSize - prevSize  > 1)
				cout<<"**** Error: Step size is too large! ********"<<endl;
			if(quit)
				break;
		}
	}

#ifdef _VERBOSEMODE_
	for(int j = 0; j < endPoints.getn(); j++){
		cout<<"Contour "<<j<<":"<<endl;
		for(int k = 0; k < endPoints[j].getn(); k++) {
			cout<<"("<<endPoints[j][k][0]<<", "<<endPoints[j][k][1]<<")"<<endl;
		}
	}
#endif

	return totalcount;
}

