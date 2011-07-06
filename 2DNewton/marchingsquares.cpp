/*
 * marchingsquares.cpp: Alternative method of generating a point set
 * Author: Patrick Butler
 * Created: 09/01/2010
 * Last Modified: 09/03/2010
 *
 */

#include "marchingsquares.h"

int generateGrid (int incHorz, int incVert, FloatArrArr &grid)
{
	grid.allocate(incHorz);
	for(int i = 0; i < incHorz; i++)
		grid[i].allocate(incVert);
	return 0;
}

int getGridValues (float upperBound, float lowerBound, 
				  float rightBound, float leftBound, int incHorz, int incVert,
				  FloatArr fcoeff, V2iArr fdegs, FloatArrArr &grid)
{
	float xdiff = (rightBound - leftBound) / incHorz;
	float ydiff = (upperBound - lowerBound) / incVert;
	V2f point;
	for(int i = 0; i < incHorz; i++)
	{
		point[0] = xdiff * i + leftBound;
		for(int j = 0; j < incVert; j++)
		{
			point[1] = ydiff * j + lowerBound;
			grid[i][j] = evalBi(fcoeff, fdegs, point);
		}
	}

	return 0;
}

int printGrid(FloatArrArr grid)
{
	cout<<"Grid:"<<endl;

	for(int j = grid[0].getn()-1; j >= 0; j--)
	{
		for(int i = 0; i < grid.getn(); i++)
		{
			cout<<grid[i][j]<<" ";
		}
		cout<<endl;
	}
	return 0;
}

int findFirstSquare (float upperBound, float lowerBound, 
				     float rightBound, float leftBound, 
					 int incHorz, int incVert, FloatArrArr grid, V2f &startSquare)
{
	// Start in the middle
	V2f initPos;
	bool found = false;
	initPos[0] = incHorz / 2;
	initPos[1] = incVert / 2;

	// Check each direction, moved towards closest to 0
	// Continue until we find a place that crosses 0
	while(!found)
	{
		int initVal = grid[initPos[0]][initPos[1]];
		int leftVal, rightVal, topVal, bottomVal;
		leftVal = grid[initPos[0]-1][initPos[1]];
		rightVal = grid[initPos[0]+1][initPos[1]];
		topVal = grid[initPos[0]][initPos[1]-1];
		bottomVal = grid[initPos[0]][initPos[1]+1];

		// Return the top-left corner of the square
		if((leftVal > 0 && initVal < 0 ) || (leftVal < 0 && initVal > 0))
		{
			initPos[0]--; // Return left position
			startSquare = initPos;
			found = true;
		}
		else if((rightVal > 0 && initVal < 0 ) || (rightVal < 0 && initVal > 0))
		{
			startSquare = initPos; // Return init position
			found = true;
		}
		else if((topVal > 0 && initVal < 0 ) || (topVal < 0 && initVal > 0))
		{
			initPos[1]++; // Return top position
			startSquare = initPos;
			found = true;
		}
		else if((bottomVal > 0 && initVal < 0 ) || (bottomVal < 0 && initVal > 0))
		{
			startSquare = initPos; // Return init position
			found = true;
		}
		else
		{
			// Using the absolute values. Which is closer to zero
			if(abs(leftVal) < abs(rightVal) && abs(leftVal) < abs(topVal) && abs(leftVal) < abs(bottomVal))
				initPos[0]--; // Move Left
			else if(abs(rightVal) < abs(leftVal) && abs(rightVal) < abs(topVal) && abs(rightVal) < abs(bottomVal))
				initPos[0]++; // Move right
			else if(abs(topVal) < abs(leftVal) && abs(topVal) < abs(rightVal) && abs(topVal) < abs(bottomVal))
				initPos[1]--; // Move top
			else
				initPos[1]++; // Move bottom
		
			// Is biggest jump towards zero better?
			/*
			int deltaLeft = abs(initVal - leftVal);
			int deltaRight = abs(initVal - rightVal);
			int deltaTop = abs(initVal - topVal);
			int deltaBottom = abs(initVal - bottomVal);
			*/
		}
	}

	return 0;
}

int nextSquare(float upperBound, float lowerBound, 
			   float rightBound, float leftBound, 
			   int incHorz, int incVert, FloatArrArr grid, 
			   V2f &nextSquare, V2fArr &endPoints)
{
	// Based on this square, return the next square
	// Which direction do we pick? (Is there 4 or 8 choices?)
		// Look at where the zero - points are, they will dictate
	//Return the top-left corner of the next square

	return 0;
}

int generateEndPoints(float upperBound, float lowerBound, 
					  float rightBound, float leftBound, int incHorz, int incVert,
				      FloatArrArr grid, V2fArr &endPoints)
{
	float xdiff = (rightBound - leftBound) / incHorz;
	float ydiff = (upperBound - lowerBound) / incVert;
	int count = 0;
	for(int i = 0; i < grid.getn(); i++)
	{
		for(int j = 0; j < grid[i].getn()-1; j++)
		{
			float dist;
			V2f point;
			if(grid[i][j] < 0 && grid[i][j+1] > 0)
			{
				dist = -grid[i][j]/(grid[i][j+1] - grid[i][j]);
				point[0] = xdiff * i + leftBound;
				point[1] = ydiff * (j + dist) + lowerBound;
				endPoints.append(point);
				count++;
			}
			if(grid[i][j] > 0 && grid[i][j+1] < 0)
			{
				dist = -grid[i][j+1]/(grid[i][j] - grid[i][j+1]);
				point[0] = xdiff * i + leftBound;
				point[1] = ydiff * (j+1 - dist) + lowerBound;
				endPoints.append(point);
				count++;
			}
		}
	}

	for(int j = 0; j < grid.getn(); j++)
	{
		for(int i = 0; i < grid[0].getn()-1; i++)
		{
			float dist;
			V2f point;
			if(grid[i][j] < 0 && grid[i+1][j] > 0)
			{
				dist = -grid[i][j]/(grid[i+1][j] - grid[i][j]);
				point[0] = xdiff * (i + dist) + leftBound;
				point[1] = ydiff * j + lowerBound;
				endPoints.append(point);
				count++;
			}
			if(grid[i][j] > 0 && grid[i+1][j] < 0)
			{
				dist = -grid[i+1][j]/(grid[i][j] - grid[i+1][j]);
				point[0] = xdiff * (i+1 - dist) + leftBound;
				point[1] = ydiff * j + lowerBound;
				endPoints.append(point);
				count++;
			}
		}
	}

	return count;
}

int marchingSquares (FloatArr fcoeff, V2iArr fdegs, float upperBound,
				   float lowerBound, float rightBound, float leftBound,
				   int incHorz, int incVert, V2fArr &endPoints)
{
	// Set up a 2-D grid from upperBound to lowerBound and rightBound to leftBound
	// with the increments specificed by incHorz and incVert.
	FloatArrArr grid;
	generateGrid(incHorz, incVert, grid);
	getGridValues(upperBound, lowerBound, rightBound, leftBound, 
				  incHorz, incVert, fcoeff, fdegs, grid);
	printGrid(grid);

	int count = generateEndPoints(upperBound, lowerBound, rightBound, leftBound, 
				      incHorz, incVert, grid, endPoints);

	// The middle steps ... ???

	// endPoints is the output array


	return count;
}
