/*
* WireFrameScene.c
*
* Author: Samuel R. Buss
*
* Bug reports: Sam Buss, sbuss@ucsd.edu.
*
* Base code for homework project #3, Math 155A, Winter 2012
*
* USAGE:   Please keep these controls in your project!!!
*
*    Press arrow keys to control view position.
*	  left and right keys rotate the viewpoint left and right
*	  Up and down keys rotate viewpoint up and down (up to 80 degrees).
*
*    Press "R" key to make step size bigger (if moving too slowly)
*    Press "r" key to make step size smaller (if moving too fast).
*
*    Press "w" key to toggle wireframe mode on and off
*
*    Press "M" , "m" to increase, decrease Mesh Resolution
*         Should apply to all appropriate objects in the scene,
*			including any glut objects (spheres, etc) and
*			including especially the (sin x)/x surface.
*
*    Press "c" to toggle culling backfaces.
*
*    Press "a" to toggle animation off and on.
*
*	  Press ESCAPE to exit.
*
*/

#include <math.h>			// For math routines (such as sqrt & trig).
#include <stdio.h>
//#include <stdlib.h>		// For the "exit" function
#include <GL/glut.h>		// OpenGL Graphics Utility Library
#include "WireFrameScene.h"
using namespace std;

// The next global variable controls the animation's state and speed.
float RotateAngle = 0.0f;		// Angle in degrees of rotation around y-axis
float Azimuth = 0.0;			// Rotated up or down by this amount
float AngleStepSize = 3.0f;		// Step three degrees at a time
/*
const float AngleStepMax = 10.0f;
const float AngleStepMin = 0.1f;
*/
const float PI = 3.14159265359f;
bool forward = true;

// Some global state variables
int MeshCount = 6;				// The mesh resolution for the mushroom top
int WireFrameOn = 1;			// == 1 for wire frame mode
int CullBackFacesOn = 0;		// == 1 if culling back faces.
const int MeshCountMin = 3;

// Animation
// Your solution for Project #3 may use different animation controls if desired.
const float mint = .75;
const float maxt = 1.25;
int animationOn = 1;
const float constTimeStep = 0.001f;  
float timeStep = constTimeStep;
float t = mint;
const float maxTime = .01f;
const float minTime = .0005f;

// glutKeyboardFunc is called below to set this function to handle
//		all "normal" key presses.
void myKeyboardFunc(unsigned char key, int x, int y)
{
	switch (key) {
	case 'a':
		animationOn = 1 - animationOn;
		break;
	case 'c':
		CullBackFacesOn = 1 - CullBackFacesOn;
		if (CullBackFacesOn) {
			glEnable(GL_CULL_FACE);				// Enable culling of back faces
		}
		else {
			glDisable(GL_CULL_FACE);				// Show all faces (front and back)
		}
		glutPostRedisplay();
		break;
	case 'm':
		MeshCount = (MeshCount>MeshCountMin) ? MeshCount - 1 : MeshCount;
		glutPostRedisplay();
		break;
	case 'M':
		MeshCount++;
		glutPostRedisplay();
		break;
	case 'w':
		WireFrameOn = 1 - WireFrameOn;
		if (WireFrameOn) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);		// Just show wireframes
		}
		else {
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);		// Show solid polygons
		}
		glutPostRedisplay();
		break;
	case 'R':
		timeStep *= 1.5;
		if (timeStep>maxTime) {
			timeStep = maxTime;
		}
		break;
	case 'r':
		timeStep /= 1.5;
		if (timeStep<minTime) {
			timeStep = minTime;
		}
		break;
	case 27:	// Escape key
		exit(1);
	}
}

// glutSpecialFunc is called below to set this function to handle
//		all "special" key presses.  See glut.h for the names of
//		special keys.
void mySpecialKeyFunc(int key, int x, int y)
{
	switch (key) {
	case GLUT_KEY_UP:
		Azimuth += AngleStepSize;
		if (Azimuth>80.0f) {
			Azimuth = 80.0f;
		}
		break;
	case GLUT_KEY_DOWN:
		Azimuth -= AngleStepSize;
		if (Azimuth < -80.0f) {
			Azimuth = -80.0f;
		}
		break;
	case GLUT_KEY_LEFT:
		RotateAngle += AngleStepSize;
		if (RotateAngle > 180.0f) {
			RotateAngle -= 360.0f;
		}
		break;
	case GLUT_KEY_RIGHT:
		RotateAngle -= AngleStepSize;
		if (RotateAngle < -180.0f) {
			RotateAngle += 360.0f;
		}
		break;
	}
	glutPostRedisplay();

}

/*
* drawScene() handles the animation and the redrawing of the
*		graphics window contents.
*/
void drawScene(void)
{
	if (animationOn && forward) {
		t += timeStep;
	}
	else if (animationOn && !forward)
	{
		t -= timeStep;
	}
	// Clear the rendering window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	// Rotate the image
	glMatrixMode(GL_MODELVIEW);			// Current matrix affects objects positions
	glLoadIdentity();						// Initialize to the identity
	glTranslatef(0.0, -2.5, -35.0);				// Translate from origin (in front of viewer)
	glRotatef(Azimuth, 1.0, 0.0, 0.0);			// Set Azimuth angle
	glRotatef(RotateAngle, 0.0, 1.0, 0.0);		// Rotate around y-axis

												// Draw the base plane (white)
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_QUADS);
	glVertex3f(-4.0, 0.0, 4.0);
	glVertex3f(4.0, 0.0, 4.0);
	glVertex3f(4.0, 0.0, -4.0);
	glVertex3f(-4.0, 0.0, -4.0);
	glEnd();

	// Draw the surface of rotation
	myDrawSurfaceOfRotation();			// You will need to write this function!

										// Draw the letter-based structure
	myDrawGeometry();					// You will need to write this function!

										// Flush the pipeline, swap the buffers
	glFlush();
	glutSwapBuffers();

	glutPostRedisplay();


}


void myDrawSurfaceOfRotation()
{
	// For project #3 you need to write this routine.
	// See the homework instructions for information 
	//     on the shape and position of the sombrero.

	// Remove the next four lines of code!
	// Replace them with code that draws
	//     the (cos r)/(5+r*r) surface of revolution.

	float DIAMETER = 3*PI;
	float theta = (2.0 * PI) / MeshCount;
	float Deltar = (DIAMETER / MeshCount);
	float r = Deltar;
	float rP = r + Deltar;
	float x, y, z;
	float theta1,theta2,x1,x2,x3,x4,y1,y2,y3, y4,z1,z2, z3, z4;


	glColor3f(1.0f, 0.3f, 1.0f);

	glPushMatrix();
	glScalef(0.25f, 8.0f, 0.25f);	

	glPushMatrix();
	glTranslatef(3.0*PI/2.0, 0.08f, 3.0 * PI / 2.0);

	glFrontFace(GL_CW);
	glBegin(GL_TRIANGLE_FAN);
	glVertex3f(0, cos(0) / (5 + (0*0)), 0);
	glVertex3f(r, cos(r) / (5 + (r * r)), 0);
	for (int i = 0; i < MeshCount; i++) {
		x = r*cos(theta);
		y = cos(r) / (5 + (r*r));
		z = r*sin(theta);
		glVertex3f(x,y,z);
		theta += (2.0 * PI) / MeshCount;
	}
	glEnd();

	glFrontFace(GL_CCW);
	//reset angle
	theta = (2.0 * PI) / MeshCount;

	for (int i = 1; i <= MeshCount; i++) {
		for (int j = 1; j <= MeshCount-1; j++) {
			theta1 = theta*i;
			theta2 = theta*(i + 1);
			x1 = r*cos(theta1);
			x2 = r*cos(theta2);
			x3 = rP*cos(theta1);
			x4 = rP*cos(theta2);
			z1 = r*sin(theta1);
			z2 = r*sin(theta2);
			z3 = rP*sin(theta1);
			z4 = rP*sin(theta2);

			glBegin(GL_QUAD_STRIP);
			glVertex3f(x1, cos(r) / (5 + (r*r)), z1);
			glVertex3f(x2, cos(r) / (5 + (r*r)), z2);
			glVertex3f(x3, cos(rP) / (5 + (rP*rP)),z3);
			glVertex3f(x4, cos(rP) / (5 + (rP*rP)),z4);
			

			r += Deltar;
			rP += Deltar;
			glEnd();
		}
		r = Deltar;
		rP = r + Deltar;

	}
	


	//
	glPopMatrix();
	glPopMatrix();


	return;

}

void myDrawGeometry()
{
	// You will need to rewrite this routine for Project #3.
	// In the supplied executable, I draw a "S" based on elongated
	//	ellipsoids.  These were formed using glutSolidSphere()'s,
	//  reoriented and scaled to turn the spheres into ellipsoids.
	// You should do something comparable -- about as complicated as
	//  my example.   You are encouraged to include a larger variety
	//  of shapes than just ellipses, and to make the end result
	//  artistically interesting in some aspect.
	// Your scene must include some animation.  

	glColor3f(0.0, 1.0, 0.5);		

	glPushMatrix();
	glRotatef(45.0,0.0f,1.0f,0.0f);

	glPushMatrix();
	glTranslatef(-2.0f, 0.0f, -2.0f);



	glPushMatrix();
	glTranslatef(0.0f, 0.8f, 0.0f);
	glPushMatrix();
	glScalef(0.1f,0.8f,0.1f);
	glutSolidSphere(1.0, MeshCount, MeshCount);
	glPopMatrix();
	glPopMatrix();


	glPushMatrix();
	glTranslatef(0.0f, 2.4f, 0.0f);
	glPushMatrix();
	glScalef(0.1f, 0.8f, 0.1f);
	glutSolidSphere(1.0, MeshCount, MeshCount);
	glPopMatrix();
	glPopMatrix();

	if (t> maxt && forward) {
		forward = !forward;
	}
	if (t < mint && !forward)
	{
		forward = !forward;
	}

	glPushMatrix();
	glTranslatef(0.0f, 2.2f, 1.6*cos(PI / 4)*0.5);
	glPushMatrix();
	glRotatef(45.0f, 1.0f, 0.0f, 0.0f);
	glPushMatrix();
	glScalef(t, t, t);
	glPushMatrix();
	glScalef(0.1f, 0.8f, 0.1f);
	glutSolidSphere(1.0, MeshCount, MeshCount);
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();


	glPushMatrix();
	glTranslatef(0.0f, 1.0f, 1.6*cos(PI / 4)*0.5);
	glPushMatrix();
	glRotatef(-45.0f, 1.0f, 0.0f, 0.0f);
	glPushMatrix();
	glScalef(t, t, t);
	glPushMatrix();
	glScalef(0.1f, 0.8f, 0.1f);
	glutSolidSphere(1.0, MeshCount, MeshCount);
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();


	glPopMatrix();
	glPopMatrix();
	return;

}

// Initialize OpenGL's rendering modes
void initRendering()
{
	glEnable(GL_DEPTH_TEST);	// Depth testing must be turned on

	glCullFace(GL_BACK);		// These two commands will cause backfaces to not be drawn

								// Possibly turn on culling of back faces.
	if (CullBackFacesOn) {
		glEnable(GL_CULL_FACE);
	}

	// Possibly turn on wireframe mode.
	if (WireFrameOn) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);		// Just show wireframes
	}
}

// Called when the window is resized
//		w, h - width and height of the window in pixels.
void resizeWindow(int w, int h)
{
	double aspectRatio;

	// Define the portion of the window used for OpenGL rendering.
	glViewport(0, 0, w, h);	// View port uses whole window

							// Set up the projection view matrix: perspective projection
							// Determine the min and max values for x and y that should appear in the window.
							// The complication is that the aspect ratio of the window may not match the
							//		aspect ratio of the scene we want to view.
	w = (w == 0) ? 1 : w;
	h = (h == 0) ? 1 : h;
	aspectRatio = (double)w / (double)h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(12.0, aspectRatio, 20.0, 50.0);

}


// Main routine
// Set up OpenGL, define the callbacks and start the main loop
int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	// We're going to animate it, so double buffer 
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	// Window position (from top corner), and size (width% and hieght)
	glutInitWindowPosition(10, 60);
	glutInitWindowSize(480, 360);
	glutCreateWindow("Wire Frame Scene Demo");

	// Initialize OpenGL as we like it..
	initRendering();

	// Set up callback functions for key presses
	glutKeyboardFunc(myKeyboardFunc);			// Handles "normal" ascii symbols
	glutSpecialFunc(mySpecialKeyFunc);		// Handles "special" keyboard keys

											// Set up the callback function for resizing windows
	glutReshapeFunc(resizeWindow);

	// Call this for background processing
	// glutIdleFunc( myIdleFunction );

	// call this whenever window needs redrawing
	glutDisplayFunc(drawScene);

	fprintf(stdout, "Arrow keys control viewpoint.n");
	fprintf(stdout, "Press \"w\" to toggle wireframe mode.\n");
	fprintf(stdout, "Press \"c\" to toggle culling of back faces.\n");
	fprintf(stdout, "Press \"M\" or \"m\" to increase or decrease resolution of mushroom cap.\n");
	fprintf(stdout, "Press \"R\" or \"r\" to increase or decrease rate of movement (respectively).\n");
	fprintf(stdout, "Press \"a\" to toggle animation on and off.\n");

	// Start the main loop.  glutMainLoop never returns.
	glutMainLoop();

	return(0);	// This line is never reached.
}