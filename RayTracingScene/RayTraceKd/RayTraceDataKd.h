/*
 *
 * RayTrace Software Package, release 2.0, February 2004.
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

class VectorR4;
class Light;
class Material;
class ViewableBase;
class CameraView;

// Camera and camera information

extern double Cpos[3];		// Position of camera
extern double Cdir[3];		// Direction of camera
extern double Cdist;		// Distance to "screen"
extern double Cdims[2];		// Width & height of "screen"

// Here are the arrays that hold information about the scene

extern SceneDescription TheScene2;

// Routines that load the data into the scene description:
void SetUpMainView();
void SetUpMaterials();
void SetUpLights();
void SetUpViewableObjects();
