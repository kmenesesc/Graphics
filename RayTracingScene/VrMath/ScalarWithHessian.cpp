//
// ScalarWithHessian.cpp
//

// Defines a class that permits computation with automatic updates
//   of the partial derivates with respect to multiple variables.
//   This is a non-sparse implementation.

#include "ScalarWithHessian.h"

ScalarWithHessian operator+( const ScalarWithHessian& a, const ScalarWithHessian& b )
{
	long len = a.GetNumberVariables();
	assert ( len==b.GetNumberVariables() );
	ScalarWithHessian ret( len );
	ret = a;
	ret += b;
	return ret;
}

ScalarWithHessian operator-( const ScalarWithHessian& a, const ScalarWithHessian& b )
{
	long len = a.GetNumberVariables();
	assert ( len==b.GetNumberVariables() );
	ScalarWithHessian ret( len );
	ret = a;
	ret -= b;
	return ret;
}

ScalarWithHessian operator*( const ScalarWithHessian& a, const ScalarWithHessian& b )
{
	long len = a.GetNumberVariables();
	assert ( len==b.GetNumberVariables() );
	ScalarWithHessian ret( len );
	ret = a;
	ret *= b;
	return ret;
}

ScalarWithHessian operator/( const ScalarWithHessian& a, const ScalarWithHessian& b )
{
	long len = a.GetNumberVariables();
	assert ( len==b.GetNumberVariables() );
	ScalarWithHessian ret( len );
	ret = a;
	ret /= b;
	return ret;
}

ScalarWithHessian operator+( const ScalarWithHessian& a, double c )
{
	long len = a.GetNumberVariables();
	ScalarWithHessian ret( len );
	ret = a;
	ret += c;
	return ret;
}

ScalarWithHessian operator-( const ScalarWithHessian& a, double c )
{
	long len = a.GetNumberVariables();
	ScalarWithHessian ret( len );
	ret = a;
	ret -= c;
	return ret;
}

ScalarWithHessian operator*( const ScalarWithHessian& a, double c )
{
	long len = a.GetNumberVariables();
	ScalarWithHessian ret( len );
	ret = a;
	ret *= c;
	return ret;
}
ScalarWithHessian operator/( const ScalarWithHessian& a, double c )
{
	assert ( c!=0.0 );
	long len = a.GetNumberVariables();
	ScalarWithHessian ret( len );
	ret = a;
	ret /= c;
	return ret;
}

ScalarWithHessian operator+( double c, const ScalarWithHessian& b )
{
	long len = b.GetNumberVariables();
	ScalarWithHessian ret( len );
	ret = b;
	ret += c;
	return ret;
}

ScalarWithHessian operator-( double c, const ScalarWithHessian& b )
{
	long len = b.GetNumberVariables();
	ScalarWithHessian ret( len );
	ret = b;
	ret.Negate();
	ret += c;
	return ret;
}

ScalarWithHessian operator*( double c, const ScalarWithHessian& b )
{
	long len = b.GetNumberVariables();
	ScalarWithHessian ret( len );
	ret = b;
	ret *= c;
	return ret;
}

ScalarWithHessian operator/( double c, const ScalarWithHessian& b )
{
	assert( b.Value()!=0.0 );
	long len = b.GetNumberVariables();
	ScalarWithHessian ret( len );
	ret = b;
	ret.Invert();
	ret *= c;
	return ret;
}
