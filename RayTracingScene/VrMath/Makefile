INCDIRS = -I/usr/openwin/include -I/usr/local/include
WARNINGS = -Wall
CFLAGS = -g $(WARNINGS) $(INCDIRS)

libVrMath.a:  Aabb.cpp Aabb.h LinearR2.cpp LinearR2.h LinearR3.cpp LinearR3.h \
		LinearR4.cpp LinearR4.h MathMisc.h  \
		Parallelepiped.cpp Parallelepiped.h \
		PolygonClip.cpp PolygonClip.h \
		PolynomialRC.cpp PolynomialRC.h \
		Quaternion.cpp Quaternion.h \
		MersenneTwister.cpp MersenneTwister.h \
		Numbers.cpp Numbers.h \
		RationalExtrapolation.cpp RationalExtrapolation.h \
		ScalarWithHessian.cpp ScalarWithHessian.h \
		ScalarWithPartials.cpp ScalarWithPartials.h \
		Slerp.cpp Slerp.h \
		Statistics.cpp Statistics.h \
		VectorRn.cpp VectorRn.h \
		Integration.cpp Integration.h \
		MatrixLowerRnn.cpp MatrixLowerRnn.h \
		MatrixRmn.cpp MatrixRmn.h \
		VrData.h 

	g++ -c $(CFLAGS) $(LIBDIRS) Aabb.cpp \
		LinearR2.cpp LinearR3.cpp LinearR4.cpp \
		Parallelepiped.cpp PolygonClip.cpp \
		PolynomialRC.cpp Quaternion.cpp \
		MersenneTwister.cpp Numbers.cpp \
		RationalExtrapolation.cpp ScalarWithHessian.cpp \
		ScalarWithPartials.cpp Slerp.cpp \
		VectorRn.cpp Integration.cpp \
		MatrixLowerRnn.cpp MatrixRmn.cpp \
		Statistics.cpp
	ar -r libVrMath.a Aabb.o LinearR2.o LinearR3.o LinearR4.o \
		Parallelepiped.o PolygonClip.o \
		PolynomialRC.o Quaternion.o \
		MersenneTwister.o Numbers.o \
		RationalExtrapolation.o Slerp.o \
		ScalarWithHessian.o ScalarWithPartials.o \
		Statistics.o VectorRn.o Integration.o \
		MatrixLowerRnn.o MatrixRmn.o


