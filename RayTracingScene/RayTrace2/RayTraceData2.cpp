/*
 *
 * RayTrace Software Package, release 3.0.  May 3, 2006.
 *
 * Author: Samuel R. Buss
 *
 * Software accompanying the book
 *		3D Computer Graphics: A Mathematical Introduction with OpenGL,
 *		by S. Buss, Cambridge University Press, 2003.
 *
 * Software is "as-is" and carries no warranty.  It may be used without
 *   restriction, but if you modify it, please change the filenames to
 *   prevent confusion between different versions.  Please acknowledge
 *   all use of the software in any publications or products based on it.
 *
 * Bug reports: Sam Buss, sbuss@ucsd.edu.
 * Web page: http://math.ucsd.edu/~sbuss/MathCG
 *
 */

#include "RayTraceData2.h"

#include "../Graphics/Material.h"
#include "../Graphics/ViewableSphere.h"
#include "../Graphics/ViewableEllipsoid.h"
#include "../Graphics/ViewableCone.h"
#include "../Graphics/ViewableTorus.h"
#include "../Graphics/ViewableTriangle.h"
#include "../Graphics/ViewableParallelogram.h"
#include "../Graphics/ViewableCylinder.h"
#include "../Graphics/ViewableParallelepiped.h"
#include "../Graphics/ViewableBezierSet.h"
#include "../Graphics/TextureCheckered.h"
#include "../Graphics/TextureBilinearXform.h"
#include "../Graphics/TextureSequence.h"
#include "../Graphics/Light.h"
#include "../Graphics/CameraView.h"

// Data for views and intersection management

double RedColor[3] = {1.0, 0.0, 0.0};
double BlueColor[3] = {0.0, 0.0, 1.0};
double BlackColor[3] = {0.0, 0.0, 0.0};

// Camera Information   *****************************************

CameraView* MainView;
double Cpos[3] = {0.0,5.0,25.0};	// Position of camera
double Cdir[3] = {0.0,-0.1,-1.0};	// Direction of camera
double Cdist = 25.0;				// Distance to "screen"
double Cdims[2] = {18.0, 10.0};		// Width & height of "screen"
// Use next two to focus on the tori
//double Cdims[2] = {2.0,2.0};		// Width & height of "screen"
//double Cdir[3] = {0.0,-4.9,-20.0};	// Direction of camera

// Data for lights  **********************************************

// Global lighting data
double GlobalAmbientLight[3] = { 0.0, 0.0, 0.0 };
double* BackgroundColor = &BlackColor[0];
VectorR3 GlobalAmbientR3(GlobalAmbientLight[0],GlobalAmbientLight[1],GlobalAmbientLight[2]);
VectorR3 BackgroundColorR3(BackgroundColor[0],BackgroundColor[1],BackgroundColor[2]);

// Array of (pointers to) lights:
const int MAX_LIGHTS = 3;
int NumLights = 0;
Light* LightArray[MAX_LIGHTS];

// Lighting values
float Lt0amb[3] = {0.0f, 0.0f, 0.0f};
float Lt0diff[3] = {1.0f, 1.0f, 1.0f};
float Lt0spec[3] = {1.0f, 1.0f, 1.0f};
float Lt0pos[3] = {7.0, 15.0f, 12.0f};

float Lt1amb[3] = {0.0f, 0.0f, 0.0f};
float Lt1diff[3] = {1.0f, 1.0f, 1.0f};
float Lt1spec[3] = {1.0f, 1.0f, 1.0f};
float Lt1pos[3] = {-7.0f, 25.0, 12.0f};

// Data for materials *******************************************

const int MAX_MATERIALS = 12;
int NumMaterials = 0;
Material* MatArray[MAX_MATERIALS];	// Array of material pointers.

// Material values

// First material for the floor
float Mat0spec[3] = {0.5f, 0.5f, 0.5f}; //greyish
float Mat0nonspec[3] = {1.0f, 0.5f, 0.3f}; //YellowWhitish
float Mat0reflect[3] = {1.0f, 0.5f, 0.3f}; 
float Mat0shiny = 300.0f;

// Second material for the floor
float Mat1spec[3] = {0.5f, 0.5f, 0.5f}; //greyish
float Mat1nonspec[3] = {0.0f, 0.02f, 0.05f};//Navy Blue
float Mat1reflect[3] = { 0.0f, 0.02f, 0.05f};
float Mat1shiny = 300.0f;

// Material for the side walls
float Mat2spec[3] = {0.2f, 0.2f, 0.2f};
float Mat2nonspec[3] = {0.3f, 0.3f, 0.0f};
float Mat2reflect[3] = {0.8f, 0.8f, 0.4f};
float Mat2shiny = 160.0f;

// Red, mixed with some blue, (i.e., magenta) slightly reflective
float Mat3spec[3] = {0.7f, 0.7f, 0.7f};
float Mat3nonspec[3] = {0.9f, 0.0f, 0.6f};
float Mat3reflect[3] = {0.1f, 0.1f, 0.1f};
float Mat3shiny = 512.0f;

float Mat4spec[3] = {1.0f, 1.0f, 1.0f};
float Mat4nonspec[3] = {0.0f, 0.0f, 0.0f};
float Mat4reflect[3] = {0.1f, 0.1f, 0.1f};
float Mat4shiny = 512.0f;

float Mat5spec[3] = {0.6f, 0.6f, 0.6f};
float Mat5nonspec[3] = {0.0f, 0.0f, 0.0f};
float Mat5reflect[3] = {0.3f, 0.3f, 0.3f};
float Mat5trans[3] = {0.8f, 0.8f, 0.8f};
float Mat5shiny = 512.0f;

float Mat6spec[3] = {0.2f, 0.2f, 0.2f};
float Mat6nonspec[3] = {0.0f, 0.2f, 0.8f};
float Mat6reflect[3] = {0.3f, 0.3f, 0.3f};
float Mat6shiny = 160.0f;

// Black!
float Mat7spec[3] = {0.6f, 0.6f, 0.6f};
float Mat7nonspec[3] = {0.0f, 0.0f, 0.0f};
float Mat7reflect[3] = {0.0f, 0.0f, 0.0f};
float Mat7shiny = 160.0f;

// Completely invisible!
float Mat8spec[3] = {0.0f, 0.0f, 0.0f};
float Mat8nonspec[3] = {0.0f, 0.0f, 0.0f};
float Mat8reflect[3] = {0.0f, 0.0f, 0.0f};
float Mat8trans[3] = {1.0f, 1.0f, 1.0f};

// A near perfect mirror
float Mat9spec[3] = {0.95f, 0.95f, 0.95f};
float Mat9nonspec[3] = {0.05f, 0.05f, 0.05f};
float Mat9reflect[3] = {0.95f, 0.95f, 0.95f};
float Mat9shiny = 160.0f;


// Data for Viewable Objects  ************************************

const int MAX_OBJECTS = 50;
int NumObjects = 0;;		// Number of viewable objects
ViewableBase* ViewObj[MAX_OBJECTS];

// Vertices of triangles & rectangles
float par0verts[3][3] = { {-8.0f,0.0f,10.0f}, {8.0f,0.0f,10.0f}, {8.0f,0.0f, -5.0f} };
float par1verts[3][3] = { {-8.0f,0.0f,-5.0f}, {8.0f,0.0f,-5.0f}, {8.0f,10.0f,-5.0f} };
float tri2verts[3][3] = { {-8.0f,0.0f,-5.0f}, {-8.0f,10.0f,-5.0f}, {-8.0f,0.0f,10.0f} };
float tri3verts[3][3] = { {8.0f,0.0f,-5.0f}, {8.0f,0.0f,10.0f}, {8.0f,10.0f,-5.0f} };
 
void SetUpMainView() {
	// Set Up Camera Information
	MainView = new CameraView();
	MainView->SetPosition( Cpos );
	MainView->SetDirection( Cdir );
	MainView->SetScreenDistance( Cdist );
	MainView->SetScreenDimensions( Cdims[0], Cdims[1] );
}

void SetUpMaterials() {
	// Initialize Array of Materials
	MatArray[0] = new Material;
	MatArray[0]->SetColorAmbientDiffuse( Mat0nonspec );
	MatArray[0]->SetColorSpecular( Mat0spec );
	MatArray[0]->SetShininess( Mat0shiny );
	MatArray[0]->SetColorReflective( Mat0reflect );

	MatArray[1] = new Material;
	MatArray[1]->SetColorAmbientDiffuse( Mat1nonspec );
	MatArray[1]->SetColorSpecular( Mat1spec );
	MatArray[1]->SetShininess( Mat1shiny );

	MatArray[2] = new Material;
	MatArray[2]->SetColorAmbientDiffuse( Mat2nonspec );
	MatArray[2]->SetColorSpecular( Mat2spec );
	MatArray[2]->SetColorReflective( Mat2reflect );
	MatArray[2]->SetShininess( Mat2shiny );

	MatArray[3] = new Material;
	MatArray[3]->SetColorAmbientDiffuse( Mat3nonspec );
	MatArray[3]->SetColorSpecular( Mat3spec );
	MatArray[3]->SetColorReflective( Mat3reflect );
	MatArray[3]->SetShininess( Mat3shiny );

	MatArray[4] = new Material;
	MatArray[4]->SetColorAmbientDiffuse( Mat4nonspec );
	MatArray[4]->SetColorSpecular( Mat4spec );
	MatArray[4]->SetColorReflective( Mat4reflect );
	MatArray[4]->SetShininess( Mat4shiny );

	MatArray[5] = new Material;
	MatArray[5]->SetColorAmbientDiffuse( Mat5nonspec );
	MatArray[5]->SetColorSpecular( Mat5spec );
	MatArray[5]->SetColorReflective( Mat5reflect );
	MatArray[5]->SetColorTransmissive( Mat5trans );
	MatArray[5]->SetShininess( Mat5shiny );
	MatArray[5]->SetIndexOfRefraction(1.5);	// Glass!

	MatArray[6] = new Material;
	MatArray[6]->SetColorAmbientDiffuse( Mat6nonspec );
	MatArray[6]->SetColorSpecular( Mat6spec );
	MatArray[6]->SetColorReflective( Mat6reflect );
	MatArray[6]->SetShininess( Mat6shiny );

	MatArray[7] = new Material;
	MatArray[7]->SetColorAmbientDiffuse( Mat7nonspec );
	MatArray[7]->SetColorSpecular( Mat7spec );
	MatArray[7]->SetColorReflective( Mat7reflect );
	MatArray[7]->SetShininess( Mat7shiny );

	// Perfectly invisible with index of
	//		refraction = 0, reflection = 0, transmission = 1.
	//	Use for two facing pieces of glass.
	MatArray[8] = new Material;
	MatArray[8]->SetColorAmbientDiffuse(Mat8nonspec);
	MatArray[8]->SetColorSpecular(Mat8spec);
	MatArray[8]->SetColorReflective( Mat8reflect );
	MatArray[8]->SetColorTransmissive( Mat8trans );

	// A near perfect mirror surface
	MatArray[9] = new Material;
	MatArray[9]->SetColorAmbientDiffuse(Mat9nonspec);
	MatArray[9]->SetColorSpecular(Mat9spec);
	MatArray[9]->SetColorReflective( Mat9reflect );
	MatArray[9]->SetShininess( Mat9shiny );

	NumMaterials = 10;
	assert( NumMaterials<=MAX_MATERIALS );
}

void SetUpLights() {
	// Global ambient light and the background color are set above.
	
	// Initialize Array of Lights
	LightArray[0] = new Light();
	LightArray[0]->SetColorAmbient( Lt0amb );
	LightArray[0]->SetColorDiffuse( Lt0diff );
	LightArray[0]->SetColorSpecular( Lt0spec );
	LightArray[0]->SetPosition( Lt0pos );
	LightArray[1] = new Light();
	LightArray[1]->SetColorAmbient( Lt1amb );
	LightArray[1]->SetColorDiffuse( Lt1diff );
	LightArray[1]->SetColorSpecular( Lt1spec );
	LightArray[1]->SetPosition( Lt1pos );

	NumLights = 2;
	assert ( NumLights<=MAX_LIGHTS );
}

void SetUpViewableObjects() {

	// Initialize array of viewable objects

	// Put in the array of visible objects

	
	// Flat plane (the floor)
	ViewableParallelogram* vp;
	vp = new ViewableParallelogram();
	vp->Init(&par0verts[0][0]);
	vp->SetMaterial(MatArray[1]);
	ViewObj[0] = vp;

	TextureCheckered* txchecked = new TextureCheckered();
	txchecked->SetMaterial1(MatArray[0]);
	txchecked->SetWidths(1.0/15.0,0.0625);
	vp->TextureMap(txchecked);

	
	/*
	ViewableSphere* vs = new ViewableSphere;
	vs->SetCenter(0.0, 0.0, 0.0);
	vs->SetRadius(5.0);
	vs->SetMaterial(MatArray[5]);
	ViewObj[1] = vs;
	*/

	/*
	
	
	// Left wall (triangular)

	ViewableTriangle* vt;
	vt = new ViewableTriangle();
	vt->Init(&tri2verts[0][0]);
	vt->SetMaterial(MatArray[2]);
	ViewObj[2] = vt;

	// Right wall (triangular)

	vt = new ViewableTriangle;
	vt->Init(&tri3verts[0][0]);
	vt->SetMaterial(MatArray[2]);
	ViewObj[3] = vt;

	// Left checkered sphere

	ViewableSphere* vs;
	vs = new ViewableSphere();
	vs->SetCenter(-7.0,0.5,-4.0);
	vs->SetRadius(0.5);
	vs->SetMaterial(MatArray[3]);
	ViewObj[4] = vs;

	TextureCheckered* txc2 = new TextureCheckered();
	txc2->SetMaterial1(MatArray[4]);
	txc2->SetWidths(0.0625,0.125);
	vs->TextureMap(txc2);	
	
	// Right checkered sphere

	vs = new ViewableSphere();
	vs->SetCenter(7.0,0.5,-4.0);
	vs->SetRadius(0.5);
	vs->SetMaterial(MatArray[3]);
	ViewObj[5] = vs;
	vs->TextureMap(txc2);
	vs->SetuvCylindrical();		// SetuvSpherical is the default

	// Two transparent spheres

	double radS = 1.8;
	double zS = 2.0;
	ViewableSphere* vsX = new ViewableSphere();
	vsX->SetRadius( radS );
	vsX->SetMaterial(MatArray[5]);

	// Left transparent sphere

	vs = new ViewableSphere();
	*vs = *vsX;
	vs->SetCenter( -radS, (sqrt(2.0)*radS), zS );
	ViewObj[6] = vs;

	// Right transparent sphere

	vs = new ViewableSphere();
	*vs = *vsX;
	vs->SetCenter( radS, (sqrt(2.0)*radS), zS );
	ViewObj[7] = vs;

	// Make lots of copies of the left small sphere.

	vsX = (ViewableSphere*)ViewObj[4];
	vs = new ViewableSphere;
	*vs = *vsX;
	vs->SetCenter( -5.0, 0.5, -4.0 );
	ViewObj[8] = vs;

	vs = new ViewableSphere;
	*vs = *vsX;
	vs->SetCenter( -3.0, 0.5, -4.0 );
	ViewObj[9] = vs;

	vs = new ViewableSphere;
	*vs = *vsX;
	vs->SetCenter( -1.0, 0.5, -4.0 );
	ViewObj[10] = vs;

	// Make lots of copies of the right small sphere

	vsX = (ViewableSphere*)ViewObj[5];
	vs = new ViewableSphere;
	*vs = *vsX;
	vs->SetCenter( 5.0, 0.5, -4.0 );
	ViewObj[11] = vs;

	vs = new ViewableSphere;
	*vs = *vsX;
	vs->SetCenter( 3.0, 0.5, -4.0 );
	ViewObj[12] = vs;

	vs = new ViewableSphere;
	*vs = *vsX;
	vs->SetCenter( 1.0, 0.5, -4.0 );
	ViewObj[13] = vs;

	// Try out two test tori

	ViewableTorus* vT = new ViewableTorus();
	vT->SetCenter( -0.3, 0.6, 7.0 );
	vT->SetMaterial( MatArray[5] );
	vT->SetRadii(0.6,0.15);
	vT->SetAxis(VectorR3(0.0,1.0,1.0));
	ViewObj[14] = vT;

	vT = new ViewableTorus();
	vT->SetCenter( 0.3, 0.6, 7.0 );
	vT->SetMaterial( MatArray[0] );
	vT->SetRadii(0.6,0.15);
	vT->SetAxis(VectorR3(0.0,1.0,-1.0));
	ViewObj[15] = vT;

	// Upright right cylinder
	VectorR3 firstCylCenter(-4.0,1.5,6.0);
	ViewableCylinder* vcyl = new ViewableCylinder();
	vcyl->SetCenter(firstCylCenter);
	vcyl->SetCenterAxis(0.0,1.0,0.0); // The default
	vcyl->SetHeight(3.0);
	vcyl->SetRadius(0.25);
	vcyl->SetMaterial(MatArray[3]);
	ViewObj[16] = vcyl;

	VectorR3 cylEndPt(-4.0,1.5,8.0);	// Common end pt. of diagonal cylinders
	VectorR3 cylnormalA(1.0,0.0,0.0);
	double coefA = cylnormalA^cylEndPt;
	VectorR3 cylnormalB(0.0,-1.0,0.0);
	vcyl = new ViewableCylinder();
	vcyl->SetCenter(cylEndPt);
	vcyl->SetCenterAxis(1.0,1.0,0.0);
	vcyl->SetRadialAxes(VectorR3(0.0,0.0,1.0),VectorR3(1.0,-1.0,0.0));
	vcyl->SetRadii(0.1,0.5);
	vcyl->SetTopFace( cylnormalA, coefA );
	vcyl->SetBottomFace( cylnormalB, 0.0 );
	vcyl->SetMaterial(MatArray[5]);
	vcyl->SetMaterialTopInner(MatArray[8]);
	vcyl->SetMaterialTopOuter(MatArray[8]);
	ViewObj[17] = vcyl;

	vcyl = new ViewableCylinder();
	vcyl->SetCenter(cylEndPt);
	vcyl->SetCenterAxis(1.0,-1.0,0.0);
	vcyl->SetRadialAxes(VectorR3(0.0,0.0,1.0),VectorR3(1.0,1.0,0.0));
	vcyl->SetRadii(0.1,0.5);
	vcyl->SetBottomFace( -cylnormalA, -coefA );
	vcyl->SetTopFace( cylnormalB, 0.0 );
	vcyl->SetMaterial(MatArray[5]);
	vcyl->SetMaterialTopInner(MatArray[8]);
	vcyl->SetMaterialTopOuter(MatArray[8]);
	ViewObj[18] = vcyl;

	vcyl = new ViewableCylinder();
	vcyl->SetCenter(-4.0,0.2,8.0);
	vcyl->SetCenterAxis(0.0,0.0,-1.0);
	vcyl->SetRadialAxes(VectorR3(0.0,1.0,0.0),VectorR3(1.0,0.0,0.0));
	vcyl->SetRadii(0.2,0.4);
	vcyl->SetHeight(1.0);
	vcyl->SetMaterial(MatArray[3]);
	ViewObj[19] = vcyl;

	TextureCheckered* txc3 = new TextureCheckered();
	txc3->SetMaterial1(MatArray[4]);
	txc3->SetWidths(0.125,0.25);
	ViewObj[19]->TextureMap(txc3);	

	// Horizontal yellow ellipsoid
	ViewableEllipsoid* ve = new ViewableEllipsoid();
	ve->SetCenter(4.0,2.0,8.0);
	ve->SetMaterial(MatArray[0]);
	ve->SetRadii(0.3,0.6,1.0);		// radii along y axis, z axis & x axis
	ViewObj[20] = ve;

	// Diagonal purple and black ellipsoid
	ve = new ViewableEllipsoid();
	ve->SetCenter(4.0,1.0,9.0);
	ve->SetAxes(VectorR3(1.0,1.0,0.0),VectorR3(0.0,0.0,1.0));
	ve->SetRadii(0.5,0.2,0.8);
	ve->SetMaterial(MatArray[3]);
	ViewObj[21] = ve;
	ViewObj[21]->TextureMap(txc2);

	// Vertical, glass ellipsoid
	ve = new ViewableEllipsoid();
	ve->SetCenter(6.5,2.5,6.0);
	ve->SetRadii(2.5,0.1,0.3);
	ve->SetMaterial(MatArray[5]);
	ViewObj[22] = ve;

	// Cone
	ViewableCone* vcone = new ViewableCone();
	vcone->SetApex(3.0,1.0,6.0);
	vcone->SetSlope(2.0);
	vcone->SetMaterial(MatArray[0]);
	ViewObj[23] = vcone;

	vcone = new ViewableCone();
	vcone->SetApex(3.0, 0.0, 9.0);
	vcone->SetCenterAxis(6.0,-1.0,0.0);
	vcone->SetRadialAxis(VectorR3(0.0,0.0,1.0));
	vcone->SetSlopes(1.0, 6.0);
	vcone->SetMaterial(MatArray[3]);
	vcone->TextureMap(txc2);
	ViewObj[24] = vcone;
	
	vcone = new ViewableCone();
	vcone->SetApex(5.0, 1.0, 9.2);
	vcone->SetCenterAxis(1.0,2.0,0.0);
	vcone->SetSlopes(4.0, 4.0);
	vcone->SetBaseFace(VectorR3(0.0,-1.0,0.0), -0.01);
	vcone->SetMaterial(MatArray[5]);
	vcone->SetMaterialBaseInner(MatArray[3]);
	ViewObj[25] = vcone;

	vcone = new ViewableCone();
	vcone->SetApex(5.0, 1.0, 9.2);
	vcone->SetCenterAxis(-1.0,2.0,0.0);
	vcone->SetSlopes(4.0, 4.0);
	vcone->SetBaseFace(VectorR3(0.0,-1.0,0.0), -0.01);
	vcone->SetMaterial(MatArray[5]);
	vcone->SetMaterialBaseInner(MatArray[3]);
	ViewObj[26] = vcone;

	// Yellow cube
	ViewableParallelepiped* vpiped = new ViewableParallelepiped();
	vpiped->SetVertices(VectorR3(-6.0,6.0,4.0), VectorR3(-5.0,6.0,4.0),
						VectorR3(-6.0,7.0,4.0), VectorR3(-6.0,6.0,3.0) );
	vpiped->SetMaterialOuter(MatArray[0]);
	ViewObj[27] = vpiped;

	// Red & Black checked parallelpiped
	vpiped = new ViewableParallelepiped();
	vpiped->SetVertices(VectorR3(-6.6,6.0,4.5), VectorR3(-6.6,1.0,5.0),
						VectorR3(-6.6,5.5,6.0), VectorR3(-6.3,5.5,4.5) );
	vpiped->SetMaterialOuter(MatArray[3]);
	vpiped->TextureMap(txc3);
	ViewObj[28] = vpiped;
	*/

	//I dug my own grave
	float verY = 0.0f;
	float mx = 3.0f;
	float my = 5.0f;
	// Bezier Patches Definitions

	//Curve 1
	double cntlPts1[4][3][4] = {
		{ {0.0f * mx, 0.80f*my+ verY, 0.0f, 1.0f}, {0.0f* mx, 0.0f, 0.0f* mx, 0.0f},{ 0.0f* mx, 0.80f*my + verY, 0.0f, 1.0f } },
		{ {-0.06f*mx, 0.8f*my + verY, 0.0f, 1.0f}, {0.0f* mx, 0.0f, 0.06f* mx, 0.0f },{ 0.06f* mx, 0.8f*my + verY, 0.0f, 1.0f } },
		{ {-0.05f* mx, 0.8f*my + verY, 0.0f, 1.0f}, {0.0f* mx, 0.0f, 0.05f* mx, 0.0f}, {0.05* mx, 0.8f*my + verY, 0.0f, 1.0f} },
		{ {-0.05f* mx, 0.78f*my + verY, 0.0f, 1.0f}, {0.0f* mx, 0.0f, 0.05f* mx, 0.0f},{ -0.05f* mx, 0.78f*my + verY, 0.0f, 1.0f } },
	};



	ViewableBezierSet* vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront( MatArray[5] );	// Glass
	vBezierSet->SetMaterialBack( MatArray[5] );	// Glass
	vBezierSet->AddRationalPatch(3,4,&cntlPts1[0][0][0]);
	ViewObj[1] = vBezierSet;

	//Curve 2
	double cntlPts2[4][3][4] = {
		{ { -0.05f* mx, 0.78f*my + verY, 0.0f, 1.0f },{ 0.0f* mx, 0.0f, 0.05f* mx, 0.0f },{ 0.05f* mx, 0.78f*my + verY, 0.0f, 1.0f } },
		{ { -0.05f* mx, 0.7f*my + verY, 0.0f, 1.0f },{ 0.0f* mx, 0.0f, 0.05f* mx, 0.0f },{ 0.05f* mx, 0.7f*my + verY, 0.0f, 1.0f } },
		{ { -0.2f* mx, 0.8f*my + verY, 0.0f, 1.0f },{ 0.0f* mx, 0.0f, 0.2f* mx, 0.0f },{ 0.2f* mx, 0.8f*my + verY, 0.0f, 1.0f } },
		{ { -0.1f* mx, 0.6f*my + verY, 0.0f, 1.0f },{ 0.0f* mx, 0.0f, 0.1f* mx, 0.0f },{ 0.1f* mx, 0.6f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts2[0][0][0]);
	ViewObj[2] = vBezierSet;

	//Curve 3
	double cntlPts3[4][3][4] = {
		{ { -0.1f* mx, 0.6f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1f* mx, 0.0f },{ 0.1f* mx, 0.6f*my + verY, 0.0f, 1.0f } },
		{ { -0.1f* mx, 0.5999f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1f* mx, 0.0f },{ 0.1f* mx, 0.5999f*my + verY, 0.0f, 1.0f } },
		{ { -0.1077f* mx, 0.5889f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1077f* mx, 0.0f },{ 0.1077f* mx, 0.5889f*my + verY, 0.0f, 1.0f } },
		{ { -0.1f* mx, 0.5868f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1f* mx, 0.0f },{ 0.1f* mx, 0.5868f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts2[0][0][0]);
	ViewObj[3] = vBezierSet;

	//Curve 4
	double cntlPts4[4][3][4] = {
		{ { -0.1f* mx, 0.5868f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1f* mx, 0.0f },{ 0.1f* mx, 0.5868f*my + verY, 0.0f, 1.0f } },
		{ { -0.109f* mx, 0.582f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.109f* mx, 0.0f },{ 0.109f* mx, 0.582f*my + verY, 0.0f, 1.0f } },
		{ { -0.1077f* mx, 0.579f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1077f* mx, 0.0f },{ 0.1077f* mx, 0.579f*my + verY, 0.0f, 1.0f } },
		{ { -0.11025f* mx, 0.57755f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.11025f* mx, 0.0f },{ 0.11025f* mx, 0.57755f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts4[0][0][0]);
	ViewObj[4] = vBezierSet;

	//Curve 5
	double cntlPts5[4][3][4] = {
		{ { -0.11025f* mx, 0.57755f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.11025f* mx, 0.0f },{ 0.11025f* mx, 0.57755f*my + verY, 0.0f, 1.0f } },
		{ { -0.1342f* mx, 0.5688f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1342f* mx, 0.0f },{ 0.1342f* mx, 0.5688f*my + verY, 0.0f, 1.0f } },
		{ { -0.5605f* mx, 0.08f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.5605f* mx, 0.0f },{ 0.5605f* mx, 0.08f*my + verY, 0.0f, 1.0f } },
		{ { -0.08f* mx, 0.555f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.08f* mx, 0.0f },{ 0.08f* mx, 0.555f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts5[0][0][0]);
	ViewObj[5] = vBezierSet;

	//Curve 6
	double cntlPts6[4][3][4] = {
		{ { -0.08f* mx, 0.555f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.08f* mx, 0.0f },{ 0.08f* mx, 0.555f*my + verY, 0.0f, 1.0f } },
		{ { -0.0217f* mx, 0.447f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.0217f* mx, 0.0f },{ 0.0217f* mx, 0.447f*my + verY, 0.0f, 1.0f } },
		{ { -0.1575f* mx, 0.3645f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1575f* mx, 0.0f },{ 0.1575f* mx, 0.3645f*my + verY, 0.0f, 1.0f } },
		{ { -0.1248f* mx, 0.349f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1248f* mx, 0.0f },{ 0.1248f* mx, 0.349f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts6[0][0][0]);
	ViewObj[6] = vBezierSet;

	//Curve 7
	double cntlPts7[4][3][4] = {
		{ { -0.1248f* mx, 0.349f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.1248f* mx, 0.0f },{ 0.1248f* mx, 0.349f*my + verY, 0.0f, 1.0f } },
		{ { -0.133f* mx, 0.228f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.133f* mx, 0.0f },{ 0.133f* mx, 0.228f*my + verY, 0.0f, 1.0f } },
		{ { -0.229f* mx, 0.189f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.229f* mx, 0.0f },{ 0.229f* mx, 0.189f*my + verY, 0.0f, 1.0f } },
		{ { -0.23f* mx, 0.0416f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.23f* mx, 0.0f },{ 0.23f* mx, 0.0416f*my+ verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts7[0][0][0]);
	ViewObj[7] = vBezierSet;

	//Curve 8
	double cntlPts8[4][3][4] = {
		{ { -0.23f* mx, 0.0416f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.23f* mx, 0.0f },{ 0.23f* mx, 0.0416f *my + verY, 0.0f, 1.0f } },
		{ { -0.2117f* mx, 0.034f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.2117* mx, 0.0f },{ 0.2117f* mx, 0.034f*my + verY, 0.0f, 1.0f } },
		{ { -0.318f* mx, 0.0f*my+ verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.318f* mx, 0.0f },{ 0.318f* mx, 0.0f*my + verY, 0.0f, 1.0f } },
		{ { -0.199f* mx, 0.0f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.199f* mx, 0.0f },{ 0.199f* mx, 0.0f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts8[0][0][0]);
	ViewObj[8] = vBezierSet;

	//Curve 9
	double cntlPts9[4][3][4] = {
		{ { 0.0f* mx, 0.80f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, 0.0f* mx, 0.0f },{ 0.0f* mx, 0.80f*my + verY, 0.0f, 1.0f } },
		{ { -0.06f* mx, 0.8f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.06f* mx, 0.0f },{ 0.06f* mx, 0.8f *my + verY, 0.0f, 1.0f } },
		{ { -0.05f* mx, 0.8f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.05f* mx, 0.0f },{ 0.05f* mx, 0.8f *my + verY, 0.0f, 1.0f } },
		{ { -0.05f* mx, 0.78f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.05f* mx, 0.0f },{ 0.05f* mx, 0.78f *my + verY, 0.0f, 1.0f } },
	};



    vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts9[0][0][0]);
	ViewObj[9] = vBezierSet;

	//Curve 10
	double cntlPts10[4][3][4] = {
		{ { -0.05f* mx, 0.78f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.05f* mx, 0.0f },{ 0.05f* mx, 0.78f*my + verY, 0.0f, 1.0f } },
		{ { -0.05f* mx, 0.7f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.05f* mx, 0.0f },{ 0.05f* mx, 0.7f*my + verY, 0.0f, 1.0f } },
		{ { -0.2f* mx, 0.8f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.2f* mx, 0.0f },{ 0.2f* mx, 0.8f*my + verY, 0.0f, 1.0f } },
		{ { -0.1f* mx, 0.6f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1f* mx, 0.0f },{ 0.1f* mx, 0.6f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts10[0][0][0]);
	ViewObj[10] = vBezierSet;

	//Curve 11
	double cntlPts11[4][3][4] = {
		{ { -0.1f* mx, 0.6f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1f* mx, 0.0f },{ 0.1f* mx, 0.6f*my + verY, 0.0f, 1.0f } },
		{ { -0.1f* mx, 0.5999f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1f* mx, 0.0f },{ 0.1f* mx, 0.5999f*my + verY, 0.0f, 1.0f } },
		{ { -0.1077f* mx, 0.5889f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1077f* mx, 0.0f },{ 0.1077f* mx, 0.5889f*my + verY, 0.0f, 1.0f } },
		{ { -0.1f* mx, 0.5868f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1f* mx, 0.0f },{ 0.1f* mx, 0.5868f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts10[0][0][0]);
	ViewObj[11] = vBezierSet;

	//Curve 12
	double cntlPts12[4][3][4] = {
		{ { -0.1f* mx, 0.5868f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1f* mx, 0.0f },{ 0.1f* mx, 0.5868f*my + verY, 0.0f, 1.0f } },
		{ { -0.109f* mx, 0.582f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.109f* mx, 0.0f },{ 0.109f* mx, 0.582f*my + verY, 0.0f, 1.0f } },
		{ { -0.1077f* mx, 0.579f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1077f* mx, 0.0f },{ 0.1077f* mx, 0.579f*my + verY, 0.0f, 1.0f } },
		{ { -0.11025f* mx, 0.57755f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.11025f* mx, 0.0f },{ 0.11025f* mx, 0.57755f *my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts12[0][0][0]);
	ViewObj[12] = vBezierSet;

	//Curve 13
	double cntlPts13[4][3][4] = {
		{ { -0.11025f* mx, 0.57755f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.11025f* mx, 0.0f },{ 0.11025f* mx, 0.57755f*my + verY, 0.0f, 1.0f } },
		{ { -0.1342f* mx, 0.5688f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1342f* mx, 0.0f },{ 0.1342f* mx, 0.5688f*my + verY, 0.0f, 1.0f } },
		{ { -0.5605f* mx, 0.08f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.5605f* mx, 0.0f },{ 0.5605f* mx, 0.08f*my + verY, 0.0f, 1.0f } },
		{ { -0.08f* mx, 0.555f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.08f* mx, 0.0f },{ 0.08f* mx, 0.555f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts13[0][0][0]);
	ViewObj[13] = vBezierSet;

	//Curve 14
	double cntlPts14[4][3][4] = {
		{ { -0.08f* mx, 0.555f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.08f* mx, 0.0f },{ 0.08f* mx, 0.555f*my + verY, 0.0f, 1.0f } },
		{ { -0.0217f* mx, 0.447f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.0217f* mx, 0.0f },{ 0.0217f* mx, 0.447f*my + verY, 0.0f, 1.0f } },
		{ { -0.1575f* mx, 0.3645f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1575f* mx, 0.0f },{ 0.1575f* mx, 0.3645f*my + verY, 0.0f, 1.0f } },
		{ { -0.1248f* mx, 0.349f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1248f* mx, 0.0f },{ 0.1248f*mx ,0.349f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts14[0][0][0]);
	ViewObj[14] = vBezierSet;

	//Curve 15
	double cntlPts15[4][3][4] = {
		{ { -0.1248f* mx, 0.349f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.1248f* mx, 0.0f },{ 0.1248f* mx, 0.349f*my + verY, 0.0f, 1.0f } },
		{ { -0.133f* mx, 0.228f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.133f* mx, 0.0f },{ 0.133f* mx, 0.228f *my + verY, 0.0f, 1.0f } },
		{ { -0.229f* mx, 0.189f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.229f* mx, 0.0f },{ 0.229f* mx, 0.189f*my + verY, 0.0f, 1.0f } },
		{ { -0.23f* mx, 0.0416f*my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.23f* mx, 0.0f },{ 0.23f* mx, 0.0416f*my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts15[0][0][0]);
	ViewObj[15] = vBezierSet;

	//Curve 16
	double cntlPts16[4][3][4] = {
		{ { -0.23f* mx, 0.0416f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.23f* mx, 0.0f },{ 0.23f* mx, 0.0416f*my + verY, 0.0f, 1.0f } },
		{ { -0.2117f* mx, 0.034f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.2117* mx, 0.0f },{ 0.2117f* mx , 0.034f *my + verY, 0.0f, 1.0f } },
		{ { -0.318f* mx, 0.0f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.318f* mx, 0.0f },{ 0.318f* mx, 0.0f *my + verY, 0.0f, 1.0f } },
		{ { -0.199f* mx, 0.0f *my + verY, 0.0f, 1.0f },{ 0.0f, 0.0f, -0.199f* mx, 0.0f },{ 0.199f* mx, 0.0f *my + verY, 0.0f, 1.0f } },
	};

	vBezierSet = new ViewableBezierSet();
	vBezierSet->SetMaterialFront(MatArray[5]);	// Glass
	vBezierSet->SetMaterialBack(MatArray[5]);	// Glass
	vBezierSet->AddRationalPatch(3, 4, &cntlPts16[0][0][0]);
	ViewObj[16] = vBezierSet;

	// Back wall

	vp = new ViewableParallelogram();
	vp->Init(&par1verts[0][0]);
	vp->SetMaterial(MatArray[9]);
	ViewObj[17] = vp;


	// Left wall (triangular)

	ViewableTriangle* vt;
	vt = new ViewableTriangle();
	vt->Init(&tri2verts[0][0]);
	vt->SetMaterial(MatArray[3]);
	ViewObj[18] = vt;

	// Right wall (triangular)

	vt = new ViewableTriangle;
	vt->Init(&tri3verts[0][0]);
	vt->SetMaterial(MatArray[3]);
	ViewObj[19] = vt;


	ViewableTorus* vT = new ViewableTorus();
	vT->SetCenter(0.0, 0.0, 0.0);
	vT->SetMaterial(MatArray[5]);
	vT->SetRadii(1.7, 0.23);
	vT->SetAxis(VectorR3(0.0, 1.0, 0.0));
	ViewObj[20] = vT;


	//Some kind of game

	// Upright right cylinder
	VectorR3 firstCylCenter(-4.0, 1.5, -2.0);
	ViewableCylinder* vcyl = new ViewableCylinder();
	vcyl->SetCenter(firstCylCenter);
	vcyl->SetCenterAxis(0.0, 1.0, 0.0); // The default
	vcyl->SetHeight(5.0);
	vcyl->SetRadius(0.4);
	vcyl->SetMaterial(MatArray[2]);
	ViewObj[21] = vcyl;

	//first ring
	vT = new ViewableTorus();
	vT->SetCenter(-4.0, 0.5, -2.0);
	vT->SetMaterial(MatArray[0]);
	vT->SetRadii(2.0, 0.5);
	vT->SetAxis(VectorR3(0.0, 1.0, 0.0));
	ViewObj[22] = vT;

	//second ring
	vT = new ViewableTorus();
	vT->SetCenter(-4.0, 1.4, -2.0);
	vT->SetMaterial(MatArray[0]);
	vT->SetRadii(1.5, 0.5);
	vT->SetAxis(VectorR3(0.0, 1.0, 0.0));
	ViewObj[23] = vT;

	//third ring
	vT = new ViewableTorus();
	vT->SetCenter(-4.0, 2.3, -2.0);
	vT->SetMaterial(MatArray[0]);
	vT->SetRadii(1.0, 0.5);
	vT->SetAxis(VectorR3(0.0, 1.0, 0.0));
	ViewObj[24] = vT;


	//again


	// Upright right cylinder
	VectorR3 firstCylCenterP(4.0, 1.5, -2.0);
	vcyl = new ViewableCylinder();
	vcyl->SetCenter(firstCylCenterP);
	vcyl->SetCenterAxis(0.0, 1.0, 0.0); // The default
	vcyl->SetHeight(5.0);
	vcyl->SetRadius(0.4);
	vcyl->SetMaterial(MatArray[2]);
	ViewObj[25] = vcyl;

	//first ring
	vT = new ViewableTorus();
	vT->SetCenter(4.0, 0.5, -2.0);
	vT->SetMaterial(MatArray[6]);
	vT->SetRadii(2.0, 0.5);
	vT->SetAxis(VectorR3(0.0, 1.0, 0.0));
	ViewObj[26] = vT;

	//second ring
	vT = new ViewableTorus();
	vT->SetCenter(4.0, 1.4, -2.0);
	vT->SetMaterial(MatArray[6]);
	vT->SetRadii(1.5, 0.5);
	vT->SetAxis(VectorR3(0.0, 1.0, 0.0));
	ViewObj[27] = vT;

	//third ring
	vT = new ViewableTorus();
	vT->SetCenter(4.0, 2.3, -2.0);
	vT->SetMaterial(MatArray[6]);
	vT->SetRadii(1.0, 0.5);
	vT->SetAxis(VectorR3(0.0, 1.0, 0.0));
	ViewObj[28] = vT;

	ViewableSphere* vs = new ViewableSphere;
	vs->SetCenter(4.0, 1.2, 4.0);
	vs->SetRadius(1.2);
	vs->SetMaterial(MatArray[3]);
	ViewObj[29] = vs;

	vs = new ViewableSphere;
	vs->SetCenter(4.0, 1.2, 4.0);
	vs->SetRadius(1.2);
	vs->SetMaterial(MatArray[3]);
	ViewObj[30] = vs;
	


	// Yes we have memory leaks  (sigh!)
	// Leaks occurs for: materials, textures, viewable objects, lights,
	//	and also for PixelArray.
	//
	//	The leaks are not too bad perhaps, but to avoid them, we should
	//	keep track of the textures too, so that they could
	//	be deallocated upon exit.
	
	NumObjects = 31;
	assert( NumObjects <= MAX_OBJECTS );

}

