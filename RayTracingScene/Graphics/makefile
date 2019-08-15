
BASEDIR = 
INCDIRS = -I/usr/openwin/include -I/usr/local/include
WARNINGS = -Wall
CFLAGS = -g $(WARNINGS) $(INCDIRS)

libGraphics: BumpMapFunction.cpp BumpMapFunction.h \
		CameraView.cpp CameraView.h DirectLight.cpp \
		DirectLight.h Light.h Material.cpp \
		Material.h MaterialBase.h MaterialCookTorrance.cpp \
		MaterialCookTorrance.h PixelArray.cpp PixelArray.h \
		RgbImage.cpp RgbImage.h TextureAffineXform.cpp \
		TextureAffineXform.h TextureBilinearXform.cpp \
		TextureBilinearXform.h TextureCheckered.cpp \
		TextureCheckered.h TextureMapBase.cpp TextureMapBase.h \
		TextureMultiFaces.cpp TextureMultiFaces.h \
		TextureRgbImage.cpp TextureRgbImage.h \
		TextureSequence.cpp TextureSequence.h \
		ViewableBase.h ViewableBezierSet.cpp ViewableBezierSet.h \
		ViewableCone.cpp ViewableCone.h ViewableCylinder.cpp \
		ViewableCylinder.h ViewableEllipsoid.cpp ViewableEllipsoid.h \
		ViewableParallelepiped.cpp ViewableParallelepiped.h \
		ViewableParallelogram.cpp ViewableParallelogram.h \
		ViewableSphere.cpp ViewableSphere.h ViewableTorus.cpp \
		ViewableTorus.h ViewableTriangle.cpp \
		ViewableTriangle.h VisiblePoint.h 
	gcc -c $(CFLAGS) $(INCDIRS) BumpMapFunction.cpp CameraView.cpp \
		DirectLight.cpp Material.cpp MaterialCookTorrance.cpp \
		PixelArray.cpp RgbImage.cpp TextureAffineXform.cpp \
		TextureBilinearXform.cpp TextureCheckered.cpp \
		TextureMapBase.cpp TextureMultiFaces.cpp TextureRgbImage.cpp \
		TextureSequence.cpp ViewableBezierSet.cpp ViewableCone.cpp \
		ViewableCylinder.cpp ViewableEllipsoid.cpp \
		ViewableParallelepiped.cpp ViewableParallelogram.cpp \
		ViewableSphere.cpp ViewableTorus.cpp ViewableTriangle.cpp
	ar -r libGraphics.a BumpMapFunction.o CameraView.o \
		DirectLight.o Material.o MaterialCookTorrance.o \
		PixelArray.o RgbImage.o TextureAffineXform.o \
		TextureBilinearXform.o TextureCheckered.o \
		TextureMapBase.o TextureMultiFaces.o TextureRgbImage.o \
		TextureSequence.o ViewableBezierSet.o ViewableCone.o \
		ViewableCylinder.o ViewableEllipsoid.o \
		ViewableParallelepiped.o ViewableParallelogram.o \
		ViewableSphere.o ViewableTorus.o ViewableTriangle.o
