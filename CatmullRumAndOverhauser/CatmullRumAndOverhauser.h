/*
*  CatmullRumAndOverhauser.c
*
*	Adaptation of "ConnectDots.c" from software accompanying the book
*		3D Computer Graphics: A Mathematical Introduction with OpenGL,
*		by S. Buss, Cambridge University Press, 2003.
*
 */

// Function prototypes

void myKeyboardFunc( unsigned char key, int x, int y );
void myMouseFunc( int button, int state, int x, int y );

void displayLines(void);
void removeFirstPoint();
void removeLastPoint();
void addNewPoint( float x, float y );

void initRendering();
void resizeWindow(int w, int h);

void drawCatmullRum();
float distPts(int x, int y);



