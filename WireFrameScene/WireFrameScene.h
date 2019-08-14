/*
* WireFrameScene.h
*
* Author: Samuel R. Buss
*
* Bug reports: Sam Buss, sbuss@ucsd.edu.
*
* Base code for homework project #3, Math 155A, Winter 2012
*/

// Function prototypes for WireFrameScene.cpp

void myKeyboardFunc(unsigned char key, int x, int y);
void mySpecialKeyFunc(int key, int x, int y);

void drawScene(void);

void initRendering();
void resizeWindow(int w, int h);

void myDrawSurfaceOfRotation();
void myDrawGeometry();