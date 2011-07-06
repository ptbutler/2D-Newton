/*
 * 2ddriver.cpp: Driver for running 2D Newton on an algebraic curve
 * Author: Patrick Butler
 * Created: 02/17/2010
 * Last Modified: 8/19/2010
 * 
 */

//#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
//#include <math.h>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define _VERBOSEMODE_

#include "Vector.h"
#include "1dnewton.h"
#include "2dnewton.h"
#include "marchingsquares.h"


#define TIMER_MS 25
V2fArrArr endPoints;
//V2fArrArr testPoints;
V3fArrArrArr scanPoints;
string comment = "";
int pointSet = 0;
bool scan;

void cleanup() 
{
	// Need to make sure there are no memory leaks
}

void handleKeypress(unsigned char key, int x, int y) 
{
	switch (key) 
	{
	case 27: // Escape key
		exit(0);
		break;
		// add left and right keys to flip through scans
	case '[':
		pointSet++;
		if(pointSet >= scanPoints.getn())
			pointSet = 0;
		break;
	case ']':
		pointSet--;
		if(pointSet < 0)
			pointSet = scanPoints.getn() - 1;
		break;
	}
}

void handleSpecial(int key, int x, int y)
{
	switch (key) 
	{
	case GLUT_LEFT_BUTTON:
		pointSet++;
		if(pointSet >= scanPoints.getn())
			pointSet = 0;
		break;
	case GLUT_RIGHT_BUTTON:
		pointSet--;
		if(pointSet < 0)
			pointSet = scanPoints.getn() - 1;
		break;
	}
}

void initRendering() 
{
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);
  glClearColor(0.65f, 0.65f, 0.65f, 1.0f);
  glPointSize(5);
/*
  testPoints.allocate(3);
  testPoints[0].allocate(1);
  testPoints[1].allocate(2);
  testPoints[2].allocate(3);
  testPoints[0][0][0] = -1;
  testPoints[0][0][1] = 0;
  testPoints[1][0][0] = -1;
  testPoints[1][0][1] = 0;
  testPoints[2][0][0] = -1;
  testPoints[2][0][1] = 0;
  testPoints[2][1][0] = 0;
  testPoints[2][1][1] = 0;
  testPoints[1][1][0] = 0;
  testPoints[1][1][1] = 0;
  testPoints[2][2][0] = 1;
  testPoints[2][2][1] = 0; */
}

void handleResize(int w, int h) 
{
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-5, 5, -5, 5);
}

void drawScene() 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_POINTS);
	if(scan)
		for(int j = 0; j < scanPoints[pointSet].getn(); j++) {
			for(int k = 0; k < scanPoints[pointSet][j].getn(); k++) {
				glVertex2f(scanPoints[pointSet][j][k][0], scanPoints[pointSet][j][k][1]);
			}
		}
	else
		for(int i = 0; i < endPoints.getn(); i++)
			for(int j = 0; j < endPoints[i].getn(); j++)
				glVertex2f(endPoints[i][j][0], endPoints[i][j][1]);

	glEnd();
/*
	glColor3f(0.0f, 1.0f, 0.0f);
	glBegin(GL_POINTS);
	if(scan)
		for(int j = 0; j < testPoints[pointSet].getn(); j++)
			glVertex2f(testPoints[pointSet][j][0], testPoints[pointSet][j][1]);
	glEnd();
*/

	glutSwapBuffers();
}

// Called every TIMER_MS milliseconds
void update(int value) 
{
  glutTimerFunc(TIMER_MS, update, 0);
  glutPostRedisplay();
}

int TestMarchingSquares(FloatArr fcoeff, V2iArr fdegs)
{
	float upperBound, lowerBound, rightBound, leftBound;
	upperBound = 5;
	lowerBound = -5;
	rightBound = 5;
	leftBound = -5;
	int incHorz = 101;
	int incVert = 101;
	return marchingSquares(fcoeff, fdegs, upperBound, lowerBound, rightBound, 
				  leftBound, incHorz, incVert, endPoints[0]);
}

int readAlgebraicCurve()
{
	string filename;
	cout<<"Please enter the name of the input file: "<<endl;
	cin>>filename;

	ifstream input;
	input.open(filename.c_str());

	int k = 0;
	while (input.good())
	{
		char line[256];
		char c;
		input.getline(line, 256);
		int terms = atoi(line);

		if (terms <= 0)
			break;

		cout<<"Curve "<<k<<":"<<endl;

		FloatArr fcoeff;
		fcoeff.allocate(terms);
		V2iArr fdegs;
		fdegs.allocate(terms);

		int i;
		for (i = 0; i < terms; i++)
		{
			input.get(line, 256, ' ');
			fcoeff[i] = atof(line);
			input.get(c);
			input.get(line, 256, ' ');
			fdegs[i][0] = atoi(line);
			input.get(c);
			input.getline(line, 256);
			fdegs[i][1] = atoi(line);
		}

		V2f guessPoint;
		input.get(line, 256, ' ');
		guessPoint[0] = atof(line);
		input.get(c);
		input.get(line, 256, ' ');
		guessPoint[1] = atof(line);
		input.get(c);

		float length;
		input.get(line, 256, ' ');
		length = atof(line);
		input.get(c);

		int numPoints;
		input.getline(line, 256);
		numPoints = atoi(line);

		int successes = TwoDNewton(fcoeff, fdegs, guessPoint, length, numPoints, endPoints);
		//int successes = TestMarchingSquares(fcoeff, fdegs);
#ifdef _VERBOSEMODE_
		for (i = 0; i < endPoints.getn(); i++)
			for(int j = 0; j < endPoints[i].getn(); j++) {
				cout<<"Next Point ("<<i<<"): ("<<endPoints[i][j][0]<<", "<<endPoints[i][j][1]<<")"<<endl;
			}
#endif
		cout<<"Curve "<<k<<" finished."<<endl;
		k++;   
	}

	return 0;
}

bool prayAndSpray(FloatArr fcoeff, V2iArr fdegs) {
	V2f guessPoint, dir;
	guessPoint[0] = -1000;
	dir[0] = 0;
	for(int i = -5; i < 6; i++) {
		V2fArr specials;
		guessPoint[1] = i;
		dir[1] = i;
		getSpecialPoints(fcoeff, fdegs, guessPoint, dir, specials);
		endPoints[0].append(endPoints[0], specials);
	}
	return true;
}

int shootAlgebraicCurve() {
	string filename;
	cout<<"Please enter the name of the input file: "<<endl;
	cin>>filename;

	ifstream input;
	input.open(filename.c_str());

	int k = 0;
	while (input.good())
	{
		char line[256];
		char c;
		input.getline(line, 256);
		int terms = atoi(line);

		if (terms <= 0)
			break;

		cout<<"Curve "<<k<<":"<<endl;

		FloatArr fcoeff;
		fcoeff.allocate(terms);
		V2iArr fdegs;
		fdegs.allocate(terms);

		int i;
		for (i = 0; i < terms; i++)
		{
			input.get(line, 256, ' ');
			fcoeff[i] = atof(line);
			input.get(c);
			input.get(line, 256, ' ');
			fdegs[i][0] = atoi(line);
			input.get(c);
			input.getline(line, 256);
			fdegs[i][1] = atoi(line);
		}

		V2f guessPoint;
		input.get(line, 256, ' ');
		guessPoint[0] = atof(line);
		input.get(c);
		input.get(line, 256, ' ');
		guessPoint[1] = atof(line);
		input.get(c);

		float length;
		input.get(line, 256, ' ');
		length = atof(line);
		input.get(c);

		int numPoints;
		input.getline(line, 256);
		numPoints = atoi(line);


		int successes = TestMarchingSquares(fcoeff, fdegs);
		//int successes = TwoDNewton(fcoeff, fdegs, guessPoint, length, numPoints, endPoints);
		//int successes = TestMarchingSquares(fcoeff, fdegs);
#ifdef _VERBOSEMODE_
		for (i = 0; i < endPoints.getn(); i++)
			cout<<"Next Point ("<<i<<"): ("<<endPoints[i][0]<<", "<<endPoints[i][1]<<")"<<endl;
#endif
		cout<<"Curve "<<k<<" finished."<<endl;
		k++;   
	}

	return 0;
}

int readSurfAxisInter()
{
	char line[256];
	string filename;
	cout<<"Please enter the name of the input file: "<<endl;
	cin>>filename;

	ifstream input;
	input.open(filename.c_str());

	int k = 0;
	while (input.good())
	{
		char c;
		input.getline(line, 256);
		int terms = atoi(line);

		if (terms <= 0)
			break;

		cout<<"Curve "<<k<<":"<<endl;

		FloatArr Fcoeff;
		FloatArr fcoeff;
		Fcoeff.allocate(terms);
		V3iArr Fdegs;
		V2iArr fdegs;
		Fdegs.allocate(terms);

		int i;
		for (i = 0; i < terms; i++)
		{
			input.get(line, 256, ' ');
			Fcoeff[i] = atof(line);
			input.get(c);
			input.get(line, 256, ' ');
			Fdegs[i][0] = atoi(line);
			input.get(c);
			input.get(line, 256, ' ');
			Fdegs[i][1] = atoi(line);
			input.get(c);
			input.getline(line, 256);
			Fdegs[i][2] = atoi(line);
		}

		int axis;
		input.get(line, 256, ' ');
		axis = atoi(line);
		input.get(c);

		float k;
		input.get(line, 256, ' ');
		k = atof(line);
		input.get(c);

		V2f guessPoint;
		input.get(line, 256, ' ');
		guessPoint[0] = atof(line);
		input.get(c);
		input.get(line, 256, ' ');
		guessPoint[1] = atof(line);
		input.get(c);

		float length;
		input.get(line, 256, ' ');
		length = atof(line);
		input.get(c);

		int numPoints;
		input.getline(line, 256);
		numPoints = atoi(line);

		polyAxisIntersect(Fcoeff, Fdegs, axis, k, fcoeff, fdegs);

		int successes = TwoDNewton(fcoeff, fdegs, guessPoint, length, numPoints, endPoints);
		//int successes = TestMarchingSquares(fcoeff, fdegs);
		for (i = 0; i < endPoints.getn(); i++)
			cout<<"Next Point ("<<i<<"): ("<<endPoints[i][0]<<", "<<endPoints[i][1]<<")"<<endl;

		cout<<"Curve "<<k<<" finished."<<endl;
		k++;   
	}

	return 0;
}

int readSurfPlaneInter()
{
	string filename;
	cout<<"Please enter the name of the input file: "<<endl;
	cin>>filename;

	ifstream input;
	input.open(filename.c_str());

	int k = 0;
	while (input.good())
	{
		char line[256];
		char c;
		input.getline(line, 256);
		int terms = atoi(line);

		if (terms <= 0)
			break;

		cout<<"Curve "<<k<<":"<<endl;

		FloatArr Fcoeff;
		FloatArr fcoeff;
		Fcoeff.allocate(terms);
		V3iArr Fdegs;
		V2iArr fdegs;
		Fdegs.allocate(terms);

		int i;
		for (i = 0; i < terms; i++)
		{
			input.get(line, 256, ' ');
			Fcoeff[i] = atof(line);
			input.get(c);
			input.get(line, 256, ' ');
			Fdegs[i][0] = atoi(line);
			input.get(c);
			input.get(line, 256, ' ');
			Fdegs[i][1] = atoi(line);
			input.get(c);
			input.getline(line, 256);
			Fdegs[i][2] = atoi(line);
		}

		V4f plane;
		input.get(line, 256, ' ');
		plane[0] = atof(line);
		input.get(c);
		input.get(line, 256, ' ');
		plane[1] = atof(line);
		input.get(c);
		input.get(line, 256, ' ');
		plane[2] = atof(line);
		input.get(c);
		input.get(line, 256, ' ');
		plane[3] = atof(line);
		input.get(c);

		V2f guessPoint;
		input.get(line, 256, ' ');
		guessPoint[0] = atof(line);
		input.get(c);
		input.get(line, 256, ' ');
		guessPoint[1] = atof(line);
		input.get(c);

		float length;
		input.get(line, 256, ' ');
		length = atof(line);
		input.get(c);

		int numPoints;
		input.getline(line, 256);
		numPoints = atoi(line);

		polyPlaneIntersect(Fcoeff, Fdegs, plane, fcoeff, fdegs);

		int successes = TwoDNewton(fcoeff, fdegs, guessPoint, length, numPoints, endPoints);
		//int successes = TestMarchingSquares(fcoeff, fdegs);
		for (i = 0; i < endPoints.getn(); i++)
			cout<<"Next Point ("<<i<<"): ("<<endPoints[i][0]<<", "<<endPoints[i][1]<<")"<<endl;

		cout<<"Curve "<<k<<" finished."<<endl;
		k++;   
	}

	return 0;
}

int printSlicesToFile() {
	string filename = "";
	cout<<"Please enter the name of the output file: "<<endl;
	cin>>filename;

	ofstream output;
	output.open(filename.c_str());

	// print out the slices in Johnstone's file format

	// possibly print out a comment: '[ comment ]'
	output<<comment<<endl;

	// #points #contours #slices separated by spaces
	int totalpoints = 0;
	int contours = 0;
	int slices = 0;
	slices = scanPoints.getn();
	for(int i = 0; i < scanPoints.getn(); i++) {
		contours += scanPoints[i].getn();
		for(int j = 0; j < scanPoints[i].getn(); j++) {
			totalpoints += scanPoints[i][j].getn();
		}
	}
	output<<totalpoints<<" "<<contours<<" "<<slices<<endl;

	// for each point
	for(int i = 0; i < scanPoints.getn(); i++) {
		for(int j = 0; j < scanPoints[i].getn(); j++) {
			for(int k = 0; k < scanPoints[i][j].getn(); k++) {
				// point in x y z separated by spaces with newline
				output<<scanPoints[i][j][k][0]<<" "<<scanPoints[i][j][k][1]<<" "<<scanPoints[i][j][k][2]<<endl;
			}
		}
	}

	// for each slice
	int curpoint = 0;
	for(int i = 0; i < scanPoints.getn(); i++) {
		// 'S #contours'
		output<<"S "<<scanPoints[i].getn()<<endl;
		// for each contour
		for(int j = 0; j < scanPoints[i].getn(); j++) {
			// 'C #points listOfPoints (where listOfPoints can be 'start#-end#')
			output<<"  C "<<scanPoints[i][j].getn()<<" "<<curpoint<<"-"<<curpoint+scanPoints[i][j].getn()-1<<endl;
			curpoint += scanPoints[i][j].getn();
		}
	}

	output.close();
	return 0;
}

int scanSurface()
{
	char line[256];
	// Read in the input file
	string filename = "";
	cout<<"Please enter the name of the input file: "<<endl;
	cin>>filename;

	ifstream input;
	input.open(filename.c_str());

	input.getline(line, 256);
	if(line[0] == '[') {
		// read in comment
		comment += line;
		comment += "\n";
		while(input.good()) {
			input.getline(line, 256);
			if(line[0] == ']')
				break;
			comment += line;
			comment += "\n";
		}
	}
	else {
		input.close();
		input.clear();
		input.open(filename.c_str());
	}

	FloatArr Fcoeff;
	V3iArr Fdegs;
	V2f guessPoint;
	float length;
	int numPoints;

	while (input.good())
	{
		char c;
		input.getline(line, 256);
		int terms = atoi(line);

		if (terms <= 0)
			break;

		Fcoeff.allocate(terms);
		Fdegs.allocate(terms);

		int i;
		for (i = 0; i < terms; i++)
		{
			input.get(line, 256, ' ');
			Fcoeff[i] = atof(line);
			input.get(c);
			input.get(line, 256, ' ');
			Fdegs[i][0] = atoi(line);
			input.get(c);
			input.get(line, 256, ' ');
			Fdegs[i][1] = atoi(line);
			input.get(c);
			input.getline(line, 256);
			Fdegs[i][2] = atoi(line);
		}
		
		input.get(line, 256, ' ');
		guessPoint[0] = atof(line);
		input.get(c);
		input.get(line, 256, ' ');
		guessPoint[1] = atof(line);
		input.get(c);
		
		input.get(line, 256, ' ');
		length = atof(line);
		input.get(c);

		input.getline(line, 256);
		numPoints = atoi(line);   
	}

	// ask for axis
	int axis = 0;
	cout<<"Please enter the axis to scan along (x = 0, y = 1, x = 2): "<<endl;
	cin>>axis;

	// ask for range
	V2f range;
	cout<<"Please enter the range to scan: "<<endl;
	cin>>range[0];
	cin>>range[1];

	// ask for # of increments
	int increments = 0;
	cout<<"Please enter the # of increments: "<<endl;
	cin>>increments;

	scanPoints.allocate(increments);
	float k = range[0];
	// for each increment
	for(int i = 0; i < increments; i++)
	{
		// get the curve
		FloatArr fcoeff;
		V2iArr fdegs;
		polyAxisIntersect(Fcoeff, Fdegs, axis, k, fcoeff, fdegs);

		string temp = "  interesected with plane ";
		if(axis == 0)
			temp += "x = ";
		else if(axis == 1)
			temp += "y = ";
		else
			temp += "z = ";
		char ktemp[32];
		sprintf(ktemp, "%f", k);
		temp += ktemp;
		comment += temp;
		comment += "\n]\n";

		// render the curve
		V2fArrArr points;
		int successes = TwoDNewton(fcoeff, fdegs, guessPoint, length, numPoints, points);
		scanPoints[i].allocate(points.getn());
		cout<<"Scan "<<i<<":"<<endl;
		for(int j = 0; j < points.getn(); j++){
			scanPoints[i][j].allocate(points[i].getn());
			for(int m = 0; m < points[j].getn(); m++) {
				if(axis == 0) {
					scanPoints[i][j][m][0] = k;
					scanPoints[i][j][m][1] = points[j][m][0];
					scanPoints[i][j][m][2] = points[j][m][1];
				}
				else if(axis == 1) {
					scanPoints[i][j][m][0] = points[j][m][0];
					scanPoints[i][j][m][1] = k;
					scanPoints[i][j][m][2] = points[j][m][1];
				}
				else {
					scanPoints[i][j][m][0] = points[j][m][0];
					scanPoints[i][j][m][1] = points[j][m][1];
					scanPoints[i][j][m][2] = k;
				}
				
				cout<<j<<": ("<<points[j][m][0]<<", "<<points[j][m][1]<<")"<<endl;
			}
		}
		// display important info (# of components, singularities, etc.)

		// set up next increments
		if(increments > 1)
			k += (range[1] - range[0]) / (increments - 1);
	}

	if(scanPoints.getn() > 0)
		pointSet = 0;

	input.close();

	printSlicesToFile();

	return 0;
}


int main(int argc, char** argv)
{
	int choice;
	cout<<"Select an Option:"<<endl
		<<"(1) Read an algebraic curve"<<endl
		<<"(2) Read a surface-axis intersection"<<endl
		<<"(3) Scan a surface along an axis"<<endl
		<<"(4) Use ray-gun shooting on algebraic curve"<<endl;
	cin>>choice;

	if(choice == 1) {
		scan = false;
		readAlgebraicCurve();
	}
	else if(choice == 2) {
		scan = false;
		readSurfAxisInter();
	}
	else if(choice == 3) {
		scan = true;
		scanSurface();
	}
	else if(choice == 4) {
		scan = false;
		shootAlgebraicCurve();
	}

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500, 500);

	glutCreateWindow("Algebraic Curve");
	initRendering();

	glutDisplayFunc(drawScene);
	glutKeyboardFunc(handleKeypress);
	glutSpecialFunc(handleSpecial);
	glutReshapeFunc(handleResize);
	glutTimerFunc(TIMER_MS, update, 0);
	atexit(cleanup);

	glutMainLoop();

	return 0;
}
