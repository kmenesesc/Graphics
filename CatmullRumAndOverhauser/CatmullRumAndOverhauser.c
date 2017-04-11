/*
 *  CatmullRumAndOverhauser.c
 *
 *	Adaptation of "ConnectDots.c" from software accompanying the book
 *		3D Computer Graphics: A Mathematical Introduction with OpenGL,
 *		by S. Buss, Cambridge University Press, 2003.
 *
 * Usage:  
 *   Left click to place a control point.  
 *		Maximum number of control points allowed is currently set at 64.
 *	 Press "f" to remove the first control point
 *	 Press "l" to remove the last control point.
 *   Press "0" for Catmull-Rum Interpolation
 *   Press "1" for Overhauser interpolation with chord length parametrization
 *   Press "2" for Overhauser interpolation with centripedal parametrization
 *	 Press escape to exit.
 */

#include "CatmullRumAndOverhauser.h"
#include <stdlib.h>
#include <GL/glut.h>
#include <stdio.h>
#include <math.h>

#define MaxNumPts 64
float PointArray[MaxNumPts][2];
int NumPts = 0;

// Window size in pixels
int WindowHeight;
int WindowWidth;

int sublines = 50;//Controls the resolution of the curves
float derivs[MaxNumPts][2]; //Holds the derivatives of each point
float timeStance[MaxNumPts];//Holds the time spent at each subcurve
float interVel[MaxNumPts][2];//Holds the intermediate velocity for the Overhauser interpolation
float vel[MaxNumPts][2]; //Holds the velocity at each point for the Overhauser Interpolation

int type = 0;

float distPts(int x, int y) {
	return sqrtf(pow(PointArray[x][0] - PointArray[y][0], 2.0f) + pow(PointArray[x][1] - PointArray[y][1], 2.0));
}

void myKeyboardFunc (unsigned char key, int x, int y)
{
	switch (key) {
	case 'f':
		removeFirstPoint();
		glutPostRedisplay();
		break;
	case 'l':
		removeLastPoint();
		glutPostRedisplay();
		break;
	case '0':
		type = 0;
		glutPostRedisplay();
		break;
	case '1':
		type = 1;
		glutPostRedisplay();
		break;
	case '2':
		type = 2;
		glutPostRedisplay();
		break;
	case 27:			// Escape key
		exit(0);
		break;
	}
}

void removeFirstPoint() {
	int i;
	if ( NumPts>0 ) {
		// Remove the first point, slide the rest down
		NumPts--;
		for ( i=0; i<NumPts; i++ ) {
			PointArray[i][0] = PointArray[i+1][0];
			PointArray[i][1] = PointArray[i+1][1];
		}
	}
}

// Left button presses place a control point.
void myMouseFunc( int button, int state, int x, int y ) {
	if ( button==GLUT_LEFT_BUTTON && state==GLUT_DOWN ) {
		float xPos = ((float)x)/((float)(WindowWidth-1));
		float yPos = ((float)y)/((float)(WindowHeight-1));

		yPos = 1.0f-yPos;			// Flip value since y position is from top row.

		if (!(xPos == PointArray[NumPts-1][0] && yPos == PointArray[NumPts-1][1])) {
			addNewPoint(xPos, yPos);
		}
		glutPostRedisplay();
	}
}

// Add a new point to the end of the list.  
// Remove the first point in the list if too many points.
void removeLastPoint() {
	if ( NumPts>0 ) {
		NumPts--;
	}
}

// Add a new point to the end of the list.  
// Remove the first point in the list if too many points.
void addNewPoint( float x, float y ) {
	if ( NumPts>=MaxNumPts ) {
		removeFirstPoint();
	}
	PointArray[NumPts][0] = x;
	PointArray[NumPts][1] = y; 
	NumPts++;
}


/*Catmull-rum Curve*/
void drawCatmullRum() {
	//Control Points of each subcurve
	float cntlpts[4][3];
	//Find the derivatives of each p_i 
	//Set the first derivative at the first point equal to zero
	derivs[0][0] = 0;
	derivs[0][1] = 0;
	//Set the first derivative of the points in between
	for (int i = 1; i < NumPts - 1; i++) {
		derivs[i][0] = ((float)1 / (float)2)*(PointArray[i + 1][0] - PointArray[i - 1][0]);
		derivs[i][1] = ((float)1 / (float)2)*(PointArray[i + 1][1] - PointArray[i - 1][1]);
	}
	//Set the last derivative at the last point equal to zero
	derivs[NumPts - 1][0] = 0;
	derivs[NumPts - 1][1] = 0;

	//Draw each subcurve
	for (int i = 0; i < NumPts - 1; i++) {
		//Set the interpolation points of each subcurve
		cntlpts[0][0] = PointArray[i][0];
		cntlpts[1][0] = PointArray[i][0] + ((float)1 / (float)3)*derivs[i][0];
		cntlpts[2][0] = PointArray[i + 1][0] - ((float)1 / (float)3)*derivs[i + 1][0];
		cntlpts[3][0] = PointArray[i + 1][0];

		cntlpts[0][1] = PointArray[i][1];
		cntlpts[1][1] = PointArray[i][1] + ((float)1 / (float)3)*derivs[i][1];
		cntlpts[2][1] = PointArray[i + 1][1] - ((float)1 / (float)3)*derivs[i + 1][1];
		cntlpts[3][1] = PointArray[i + 1][1];

		cntlpts[0][2] = 0;
		cntlpts[1][2] = 0;
		cntlpts[2][2] = 0;
		cntlpts[3][2] = 0;

		//Generate the Bezier Curve
		glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, &cntlpts[0][0]);

		//Draw the curve on the screen
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j <= sublines; j++) {
			glEvalCoord1f((float)j / (float)sublines);
		}
		glEnd();
	}
}

/*Overhauser chord parametrization draw curve*/
void drawOverChord() {
	float cntlpts[4][3]; //Control points of each subcurve

						 //calculate the u_i's
	timeStance[0] = 0;
	for (int i = 1; i < NumPts; i++) {
		timeStance[i] = timeStance[i - 1] + distPts(i, i - 1);
	}

	//calculate the intermediate velocity between poitns p_i and p_i+1
	for (int i = 0; i < NumPts; i++) {
		interVel[i][0] = (PointArray[i + 1][0] - PointArray[i][0]) / (timeStance[i + 1] - timeStance[i]);
		interVel[i][1] = (PointArray[i + 1][1] - PointArray[i][1]) / (timeStance[i + 1] - timeStance[i]);
	}

	//calculate the velocity at each point
	vel[0][0] = 0;
	vel[0][1] = 0;
	for (int i = 1; i < NumPts - 1; i++) {
		vel[i][0] = ((timeStance[i + 1] - timeStance[i]) / (timeStance[i + 1] - timeStance[i - 1]))*interVel[i - 1][0] +
			((timeStance[i] - timeStance[i - 1]) / (timeStance[i + 1] - timeStance[i - 1]))*interVel[i][0];
		vel[i][1] = ((timeStance[i + 1] - timeStance[i]) / (timeStance[i + 1] - timeStance[i - 1]))*interVel[i - 1][1] +
			((timeStance[i] - timeStance[i - 1]) / (timeStance[i + 1] - timeStance[i - 1]))*interVel[i][1];
	}
	vel[NumPts - 1][0] = 0;
	vel[NumPts - 1][1] = 0;

	//Draw each subcurve
	for (int i = 0; i < NumPts - 1; i++) {
		//Set the interpolation points of each subcurve
		cntlpts[0][0] = PointArray[i][0];
		cntlpts[1][0] = PointArray[i][0] + (((float)1 / (float)3)*vel[i][0] * (timeStance[i + 1] - timeStance[i]));
		cntlpts[2][0] = PointArray[i + 1][0] - (((float)1 / (float)3)*vel[i + 1][0] * (timeStance[i + 1] - timeStance[i]));
		cntlpts[3][0] = PointArray[i + 1][0];

		cntlpts[0][1] = PointArray[i][1];
		cntlpts[1][1] = PointArray[i][1] + (((float)1 / (float)3)*vel[i][1] * (timeStance[i + 1] - timeStance[i]));
		cntlpts[2][1] = PointArray[i + 1][1] - (((float)1 / (float)3)*vel[i + 1][1] * (timeStance[i + 1] - timeStance[i]));
		cntlpts[3][1] = PointArray[i + 1][1];

		cntlpts[0][2] = 0;
		cntlpts[1][2] = 0;
		cntlpts[2][2] = 0;
		cntlpts[3][2] = 0;

		//Generate the Bezier Curve
		glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, &cntlpts[0][0]);

		//Draw the curve on the screen
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j <= sublines; j++) {
			glEvalCoord1f((float)j / (float)sublines);
		}
		glEnd();
	}
}

/*Overhauser centripetal parametrization draw curve*/
void drawOverCent() {
	float cntlpts[4][3]; //Control points of each subcurve

						 //calculate the u_i's
	timeStance[0] = 0;
	for (int i = 1; i < NumPts; i++) {
		timeStance[i] = timeStance[i - 1] + sqrtf(distPts(i, i - 1));
	}

	//calculate the intermediate velocity between poitns p_i and p_i+1
	for (int i = 0; i < NumPts; i++) {
		interVel[i][0] = (PointArray[i + 1][0] - PointArray[i][0]) / (timeStance[i + 1] - timeStance[i]);
		interVel[i][1] = (PointArray[i + 1][1] - PointArray[i][1]) / (timeStance[i + 1] - timeStance[i]);
	}

	//calculate the velocity at each point
	vel[0][0] = 0;
	vel[0][1] = 0;
	for (int i = 1; i < NumPts - 1; i++) {
		vel[i][0] = ((timeStance[i + 1] - timeStance[i]) / (timeStance[i + 1] - timeStance[i - 1]))*interVel[i - 1][0] +
			((timeStance[i] - timeStance[i - 1]) / (timeStance[i + 1] - timeStance[i - 1]))*interVel[i][0];
		vel[i][1] = ((timeStance[i + 1] - timeStance[i]) / (timeStance[i + 1] - timeStance[i - 1]))*interVel[i - 1][1] +
			((timeStance[i] - timeStance[i - 1]) / (timeStance[i + 1] - timeStance[i - 1]))*interVel[i][1];
	}
	vel[NumPts - 1][0] = 0;
	vel[NumPts - 1][1] = 0;

	//Draw each subcurve
	for (int i = 0; i < NumPts - 1; i++) {
		//Set the interpolation points of each subcurve
		cntlpts[0][0] = PointArray[i][0];
		cntlpts[1][0] = PointArray[i][0] + (((float)1 / (float)3)*vel[i][0] * (timeStance[i + 1] - timeStance[i]));
		cntlpts[2][0] = PointArray[i + 1][0] - (((float)1 / (float)3)*vel[i + 1][0] * (timeStance[i + 1] - timeStance[i]));
		cntlpts[3][0] = PointArray[i + 1][0];

		cntlpts[0][1] = PointArray[i][1];
		cntlpts[1][1] = PointArray[i][1] + (((float)1 / (float)3)*vel[i][1] * (timeStance[i + 1] - timeStance[i]));
		cntlpts[2][1] = PointArray[i + 1][1] - (((float)1 / (float)3)*vel[i + 1][1] * (timeStance[i + 1] - timeStance[i]));
		cntlpts[3][1] = PointArray[i + 1][1];

		cntlpts[0][2] = 0;
		cntlpts[1][2] = 0;
		cntlpts[2][2] = 0;
		cntlpts[3][2] = 0;

		//Generate the Bezier Curve
		glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, &cntlpts[0][0]);

		//Draw the curve on the screen
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j <= sublines; j++) {
			glEvalCoord1f((float)j / (float)sublines);
		}
		glEnd();
	}
}


void displayLines(void)
{

	glClear(GL_COLOR_BUFFER_BIT);

	
	// Draw the line segments
	glColor3f(1.0f, 0.0f, 0.8f);			// Reddish/purple lines

	if (NumPts > 1) {
		if (type == 0) {
			drawCatmullRum();
		}
		else if (type == 1) {
			drawOverChord();
		}
		else if (type == 2) {
			drawOverCent();
		}
	}

	// Draw the interpolated points second.
	glColor3f( 0.0f, 0.0f, 0.0f);			// Draw points in black
	glBegin( GL_POINTS );
	for ( int i=0; i<NumPts; i++ ) {
	   glVertex2f( PointArray[i][0], PointArray[i][1] );
	}
	glEnd();

	glFlush();
}

void initRendering() {
	glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );

	// Make big points and wide lines.  (This may be commented out if desired.)
	glPointSize(8);
	glLineWidth(5);

	// The following commands should induce OpenGL to create round points and 
	//	antialias points and lines.  (This is implementation dependent unfortunately, and
	//  may slow down rendering considerably.)
	//  You may comment these out if you wish.
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);	// Make round points, not square points
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);		// Antialias the lines
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_MAP1_VERTEX_3);
}

void resizeWindow(int w, int h)
{
	WindowHeight = (h>1) ? h : 2;
	WindowWidth = (w>1) ? w : 2;
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0f, 1.0f, 0.0f, 1.0f);  // Always view [0,1]x[0,1].
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB ); 
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);
	initRendering();

	glutDisplayFunc(displayLines);
	glutReshapeFunc(resizeWindow);
	glutKeyboardFunc(myKeyboardFunc);
	glutMouseFunc(myMouseFunc);
	glutMainLoop();

	return 0;					// This line is never reached
}
