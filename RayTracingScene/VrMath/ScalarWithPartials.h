//
// ScalarWithPartials.h
//

#ifndef SCALAR_WITH_PARTIALS_H
#define SCALAR_WITH_PARTIALS_H

// Defines a class that permits computation with automatic updates
//   of the partial derivates with respect to multiple variables.
//   This is a non-sparse implementation.
class ScalarWithPartials;

#include "VectorRn.h"
#include "MathMisc.h"

class ScalarWithPartials
{
public:
	ScalarWithPartials( long numPartials );

	double Value() const { return TheValue; }
	VectorRn& Partials() { return ThePartials; }
	const VectorRn& Partials() const { return ThePartials; }
	long GetNumberPartials() const { return ThePartials.GetLength(); }

	void SetAsVariable ( double value, long variableNumber );

	void SetSum( const ScalarWithPartials& a, const ScalarWithPartials& b );
	void SetDifference( const ScalarWithPartials& a, const ScalarWithPartials& b );
	void SetProduct( const ScalarWithPartials& a, const ScalarWithPartials& b );
	void SetFraction( const ScalarWithPartials& a, const ScalarWithPartials& b );
	void SetInverse( const ScalarWithPartials& a );
	void SetSquare( const ScalarWithPartials& a );
	void SetCube( const ScalarWithPartials& a );
	void SetSqrt( const ScalarWithPartials& a );

	void Negate();
	void Square();
	void Sqrt();
	void Invert();

	ScalarWithPartials& operator= ( double c );
	ScalarWithPartials& operator += ( double c );
	ScalarWithPartials& operator -= ( double c );
	ScalarWithPartials& operator *= ( double c );
	ScalarWithPartials& operator /= ( double c );

	ScalarWithPartials& operator += ( const ScalarWithPartials& a );
	ScalarWithPartials& operator -= ( const ScalarWithPartials& a );
	ScalarWithPartials& operator *= ( const ScalarWithPartials& a );
	ScalarWithPartials& operator /= ( const ScalarWithPartials& a );
	ScalarWithPartials& SetScaled ( const ScalarWithPartials& a, double factor );
	ScalarWithPartials& AddScaled ( const ScalarWithPartials& a, double factor );

private:
	double TheValue;
	VectorRn ThePartials;
};

ScalarWithPartials operator+( const ScalarWithPartials& a, const ScalarWithPartials& b );
ScalarWithPartials operator-( const ScalarWithPartials& a, const ScalarWithPartials& b );
ScalarWithPartials operator*( const ScalarWithPartials& a, const ScalarWithPartials& b );
ScalarWithPartials operator/( const ScalarWithPartials& a, const ScalarWithPartials& b );
ScalarWithPartials operator+( const ScalarWithPartials& a, double c );
ScalarWithPartials operator-( const ScalarWithPartials& a, double c );
ScalarWithPartials operator*( const ScalarWithPartials& a, double c );
ScalarWithPartials operator/( const ScalarWithPartials& a, double c );
ScalarWithPartials operator+( double c, const ScalarWithPartials& b );
ScalarWithPartials operator-( double c, const ScalarWithPartials& b );
ScalarWithPartials operator*( double c, const ScalarWithPartials& b );
ScalarWithPartials operator/( double c, const ScalarWithPartials& b );

// *******************************************************************
// Inlined member functions
// *******************************************************************

inline ScalarWithPartials::ScalarWithPartials( long numPartials )
: ThePartials(numPartials)
{
	assert ( numPartials>=0 );
	TheValue = 0.0;
	ThePartials.Fill(0.0);
}

inline void ScalarWithPartials::SetAsVariable ( double value, long variableNumber )
{
	assert (variableNumber<GetNumberPartials());
	TheValue = value;
	ThePartials.Fill(0.0);
	ThePartials[variableNumber] = 1.0;
}

inline void ScalarWithPartials::SetSum( const ScalarWithPartials& a, const ScalarWithPartials& b )
{
	assert ( GetNumberPartials()==a.GetNumberPartials() && GetNumberPartials()==b.GetNumberPartials() );
	TheValue = a.TheValue + b.TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	const double* pPtrB = b.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++, pPtrB++ ) {
		*pPtr = (*pPtrA)+(*pPtrB);
	}
}

inline void ScalarWithPartials::SetDifference( const ScalarWithPartials& a, const ScalarWithPartials& b )
{
	assert ( GetNumberPartials()==a.GetNumberPartials() && GetNumberPartials()==b.GetNumberPartials() );
	TheValue = a.TheValue - b.TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	const double* pPtrB = b.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++, pPtrB++ ) {
		*pPtr = (*pPtrA)-(*pPtrB);
	}
}



inline void ScalarWithPartials::SetProduct( const ScalarWithPartials& a, const ScalarWithPartials& b )
{
	assert ( GetNumberPartials()==a.GetNumberPartials() && GetNumberPartials()==b.GetNumberPartials() );
	TheValue = a.TheValue * b.TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	const double* pPtrB = b.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++, pPtrB++ ) {
		*pPtr = (*pPtrA)*b.TheValue + a.TheValue*(*pPtrB);
	}
}

inline void ScalarWithPartials::SetFraction( const ScalarWithPartials& a, const ScalarWithPartials& b )
{
	assert ( GetNumberPartials()==a.GetNumberPartials() && GetNumberPartials()==b.GetNumberPartials() );
	assert ( b.TheValue!=0.0 && ::Square(b.TheValue)!=0.0 );
	TheValue = a.TheValue / b.TheValue;
	double denomSqInv = 1.0/::Square(b.TheValue);
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	const double* pPtrB = b.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++, pPtrB++ ) {
		*pPtr = ((*pPtrA)*b.TheValue - a.TheValue*(*pPtrB))*denomSqInv;
	}
}

inline void ScalarWithPartials::SetInverse( const ScalarWithPartials& a )
{
	assert ( GetNumberPartials()==a.GetNumberPartials() );
	assert ( a.TheValue!=0.0 && ::Square(a.TheValue)!=0.0 );
	TheValue = 1.0 / a.TheValue;
	double denomSqInv = 1.0/::Square(a.TheValue);
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr = -(*pPtrA)*denomSqInv;
	}
}

inline void ScalarWithPartials::SetSquare( const ScalarWithPartials& a )
{
	assert ( GetNumberPartials()==a.GetNumberPartials() );
	TheValue = ::Square(a.TheValue);
	double twoValue = 2.0*a.TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr = twoValue*(*pPtrA);
	}
}

inline void ScalarWithPartials::SetCube( const ScalarWithPartials& a )
{
	assert ( GetNumberPartials()==a.GetNumberPartials() );
	TheValue = Cube(a.TheValue);
	double threeValueSq = 3.0*::Square(a.TheValue);
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr = threeValueSq*(*pPtrA);
	}
}

inline void ScalarWithPartials::SetSqrt( const ScalarWithPartials& a )
{
	assert ( GetNumberPartials()==a.GetNumberPartials() );
	assert ( a.TheValue>0.0 );
	TheValue = sqrt(a.TheValue);
	double halfSqrtInv = 0.5/TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr = halfSqrtInv*(*pPtrA);
	}
}

inline void ScalarWithPartials::Negate()
{
	TheValue = -TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr = -(*pPtr);
	}
}

inline void ScalarWithPartials::Square()
{
	double twoValue = 2.0*TheValue;
	TheValue = ::Square(TheValue);
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= twoValue;
	}
}

inline void ScalarWithPartials::Sqrt()
{
	assert ( TheValue > 0.0 );
	TheValue = sqrt(TheValue);
	double halfSqrtInv = 0.5/TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= halfSqrtInv;
	}
}

inline void ScalarWithPartials::Invert()
{
	assert ( TheValue > 0.0 );
	TheValue = 1.0/TheValue;
	double negDenomInv = -::Square(TheValue);
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= negDenomInv;
	}
}


inline ScalarWithPartials& ScalarWithPartials::operator= ( double c )
{
	TheValue = c;
	ThePartials.SetZero();
	return *this;
}


inline ScalarWithPartials& ScalarWithPartials::operator += ( double c )
{
	TheValue += c;
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::operator -= ( double c )
{
	TheValue -= c;
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::operator *= ( double c )
{
	TheValue *= c;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= c;
	}
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::operator /= ( double c )
{
	assert ( c!=0.0 );
	double cInv = 1.0/c;
	TheValue *= cInv;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= cInv;
	}
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::operator += ( const ScalarWithPartials& a )
{
	TheValue += a.TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr += (*pPtrA);
	}
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::SetScaled( const ScalarWithPartials& a, double factor )
{
	TheValue = factor*a.TheValue;
	ThePartials.SetScaled( a.ThePartials, factor );
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::AddScaled( const ScalarWithPartials& a, double factor )
{
	TheValue += factor*a.TheValue;
	ThePartials.AddScaled( a.ThePartials, factor );
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::operator -= ( const ScalarWithPartials& a )
{
	TheValue -= a.TheValue;
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr -= (*pPtrA);
	}
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::operator *= ( const ScalarWithPartials& a )
{
	long num = ThePartials.GetLength();
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr = (*pPtr)*a.TheValue + TheValue*(*pPtrA);
	}
	TheValue *= a.TheValue;
	return *this;
}

inline ScalarWithPartials& ScalarWithPartials::operator /= ( const ScalarWithPartials& a )
{
	assert ( a.TheValue!=0.0 );
	long num = ThePartials.GetLength();
	double denomInv = 1.0/::Square(a.TheValue);
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr = ((*pPtr)*a.TheValue - TheValue*(*pPtrA)) * denomInv;
	}
	TheValue /= a.TheValue;
	return *this;
}

#endif // SCALAR_WITH_PARTIALS_H
