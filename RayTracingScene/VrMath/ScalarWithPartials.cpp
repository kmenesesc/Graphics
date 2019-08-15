//
// ScalarWithPartials.cpp
//

// Defines a class that permits computation with automatic updates
//   of the partial derivates with respect to multiple variables.
//   This is a non-sparse implementation.

#include "ScalarWithPartials.h"


ScalarWithPartials operator+( const ScalarWithPartials& a, const ScalarWithPartials& b )
{
	long len = a.GetNumberPartials();
	assert ( len==b.GetNumberPartials() );
	ScalarWithPartials ret( len );
	ret = a;
	ret += b;
	return ret;
}

ScalarWithPartials operator-( const ScalarWithPartials& a, const ScalarWithPartials& b )
{
	long len = a.GetNumberPartials();
	assert ( len==b.GetNumberPartials() );
	ScalarWithPartials ret( len );
	ret = a;
	ret -= b;
	return ret;
}

ScalarWithPartials operator*( const ScalarWithPartials& a, const ScalarWithPartials& b )
{
	long len = a.GetNumberPartials();
	assert ( len==b.GetNumberPartials() );
	ScalarWithPartials ret( len );
	ret = a;
	ret *= b;
	return ret;
}

ScalarWithPartials operator/( const ScalarWithPartials& a, const ScalarWithPartials& b )
{
	long len = a.GetNumberPartials();
	assert ( len==b.GetNumberPartials() );
	ScalarWithPartials ret( len );
	ret = a;
	ret /= b;
	return ret;
}

ScalarWithPartials operator+( const ScalarWithPartials& a, double c )
{
	long len = a.GetNumberPartials();
	ScalarWithPartials ret( len );
	ret = a;
	ret += c;
	return ret;
}

ScalarWithPartials operator-( const ScalarWithPartials& a, double c )
{
	long len = a.GetNumberPartials();
	ScalarWithPartials ret( len );
	ret = a;
	ret -= c;
	return ret;
}

ScalarWithPartials operator*( const ScalarWithPartials& a, double c )
{
	long len = a.GetNumberPartials();
	ScalarWithPartials ret( len );
	ret = a;
	ret *= c;
	return ret;
}
ScalarWithPartials operator/( const ScalarWithPartials& a, double c )
{
	assert ( c!=0.0 );
	long len = a.GetNumberPartials();
	ScalarWithPartials ret( len );
	ret = a;
	ret /= c;
	return ret;
}

ScalarWithPartials operator+( double c, const ScalarWithPartials& b )
{
	long len = b.GetNumberPartials();
	ScalarWithPartials ret( len );
	ret = b;
	ret += c;
	return ret;
}

ScalarWithPartials operator-( double c, const ScalarWithPartials& b )
{
	long len = b.GetNumberPartials();
	ScalarWithPartials ret( len );
	ret = b;
	ret.Negate();
	ret += c;
	return ret;
}

ScalarWithPartials operator*( double c, const ScalarWithPartials& b )
{
	long len = b.GetNumberPartials();
	ScalarWithPartials ret( len );
	ret = b;
	ret *= c;
	return ret;
}

ScalarWithPartials operator/( double c, const ScalarWithPartials& b )
{
	assert( b.Value()!=0.0 );
	long len = b.GetNumberPartials();
	ScalarWithPartials ret( len );
	ret.SetInverse(b);
	ret *= c;
	return ret;
}
