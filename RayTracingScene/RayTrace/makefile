
GRAPHICSDIR = ../Graphics
VRMATHDIR = ../VrMath
INCDIRS = -I/usr/openwin/include -I/usr/local/include
WARNINGS = -Wall
CFLAGS = -g $(WARNINGS) $(INCDIRS)
LIBDIRS = -L/usr/openwin/lib -L/usr/local/lib 
LIBS = -lX11 -lXi -lXmu -lglut -lGL -lGLU -lm

RayTrace: RayTrace.cpp RayTraceData.cpp RayTraceData.h \
		$(GRAPHICSDIR)/libGraphics.a $(VRMATHDIR)/libVrMath.a
	g++ $(CFLAGS) $(INCDIRS) RayTrace.cpp RayTraceData.cpp $(GRAPHICSDIR)/libGraphics.a $(VRMATHDIR)/libVrMath.a $(LIBDIRS) $(LIBS) -o RayTrace
