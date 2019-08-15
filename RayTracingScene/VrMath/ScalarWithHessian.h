//
// ScalarWithHessian.h
//

#ifndef SCALAR_WITH_HESSIAN_H
#define SCALAR_WITH_HESSIAN_H

// Defines a class that permits computation with automatic updates
//   of the partial derivates with respect to multiple variables.
//   This is a non-sparse implementation.
class ScalarWithHessian;

#include "MathMisc.h"
#include "VectorRn.h"
#include "MatrixRmn.h"

class ScalarWithHessian
{
public:
	ScalarWithHessian( long numVariables );

	long GetNumberVariables() const { return ThePartials.GetLength(); }
	double Value() const { return TheValue; }
	VectorRn& Gradient() { return ThePartials; }
	const VectorRn& Gradient() const { return ThePartials; }
	VectorRn& Partials() { return ThePartials; }
	const VectorRn& Partials() const { return ThePartials; }
	MatrixRmn& Hessian() { return TheHessian; }
	const MatrixRmn& Hessian() const { return TheHessian; }

	void SetAsVariable ( double value, long variableNumber );

	void SetSum( const ScalarWithHessian& a, const ScalarWithHessian& b );
	void SetDifference( const ScalarWithHessian& a, const ScalarWithHessian& b );
	void SetProduct( const ScalarWithHessian& a, const ScalarWithHessian& b );
	void SetFraction( const ScalarWithHessian& a, const ScalarWithHessian& b );
	void SetInverse( const ScalarWithHessian& a );
	void SetSquare( const ScalarWithHessian& a );
	void SetCube( const ScalarWithHessian& a );
	void SetSqrt( const ScalarWithHessian& a );

	void Negate();
	void Square();
	void Cube();
	void Sqrt();
	void Invert();

	ScalarWithHessian& operator= ( double c );
	ScalarWithHessian& operator += ( double c );
	ScalarWithHessian& operator -= ( double c );
	ScalarWithHessian& operator *= ( double c );
	ScalarWithHessian& operator /= ( double c );

	ScalarWithHessian& operator += ( const ScalarWithHessian& a );
	ScalarWithHessian& operator -= ( const ScalarWithHessian& a );
	ScalarWithHessian& operator *= ( const ScalarWithHessian& a );
	ScalarWithHessian& operator /= ( const ScalarWithHessian& a );
	ScalarWithHessian& SetScaled ( const ScalarWithHessian& a, double factor );
	ScalarWithHessian& AddScaled ( const ScalarWithHessian& a, double factor );

private:
	double TheValue;
	VectorRn ThePartials;
	MatrixRmn TheHessian;
};

ScalarWithHessian operator+( const ScalarWithHessian& a, const ScalarWithHessian& b );
ScalarWithHessian operator-( const ScalarWithHessian& a, const ScalarWithHessian& b );
ScalarWithHessian operator*( const ScalarWithHessian& a, const ScalarWithHessian& b );
ScalarWithHessian operator/( const ScalarWithHessian& a, const ScalarWithHessian& b );
ScalarWithHessian operator+( const ScalarWithHessian& a, double c );
ScalarWithHessian operator-( const ScalarWithHessian& a, double c );
ScalarWithHessian operator*( const ScalarWithHessian& a, double c );
ScalarWithHessian operator/( const ScalarWithHessian& a, double c );
ScalarWithHessian operator+( double c, const ScalarWithHessian& b );
ScalarWithHessian operator-( double c, const ScalarWithHessian& b );
ScalarWithHessian operator*( double c, const ScalarWithHessian& b );
ScalarWithHessian operator/( double c, const ScalarWithHessian& b );

// *******************************************************************
// Inlined member functions
// *******************************************************************

inline ScalarWithHessian::ScalarWithHessian( long numVariables )
: ThePartials(numVariables), TheHessian( numVariables, numVariables )
{
	assert ( numVariables>=0 );
	TheValue = 0.0;
	ThePartials.SetZero();
	TheHessian.SetZero();
}

inline void ScalarWithHessian::SetAsVariable ( double value, long variableNumber )
{
	assert (variableNumber<GetNumberVariables());
	TheValue = value;
	ThePartials.SetZero();
	ThePartials[variableNumber] = 1.0;
	TheHessian.SetZero();
}

inline void ScalarWithHessian::SetSum( const ScalarWithHessian& a, const ScalarWithHessian& b )
{
	assert ( GetNumberVariables()==a.GetNumberVariables() && GetNumberVariables()==b.GetNumberVariables() );
	*this = a;
	*this += b;
}

inline void ScalarWithHessian::SetDifference( const ScalarWithHessian& a, const ScalarWithHessian& b )
{
	assert ( GetNumberVariables()==a.GetNumberVariables() && GetNumberVariables()==b.GetNumberVariables() );
	*this = a;
	*this -= b;
}

inline void ScalarWithHessian::SetProduct( const ScalarWithHessian& a, const ScalarWithHessian& b )
{
	assert ( GetNumberVariables()==a.GetNumberVariables() && GetNumberVariables()==b.GetNumberVariables() );
	*this = a;
	*this *= b;
}

inline void ScalarWithHessian::SetFraction( const ScalarWithHessian& a, const ScalarWithHessian& b )
{
	assert ( GetNumberVariables()==a.GetNumberVariables() && GetNumberVariables()==b.GetNumberVariables() );
	assert ( b.TheValue!=0.0 && ::Square(b.TheValue)!=0.0 );
	*this = a;
	*this /= b;
}

inline void ScalarWithHessian::SetInverse( const ScalarWithHessian& a )
{
	assert ( GetNumberVariables()==a.GetNumberVariables() );
	assert ( a.TheValue!=0.0 && ::Square(a.TheValue)!=0.0 );
	*this = a;
	this->Invert();
}

inline void ScalarWithHessian::SetSquare( const ScalarWithHessian& a )
{
	assert ( GetNumberVariables()==a.GetNumberVariables() );
	*this = a;
	this->Square();
}

inline void ScalarWithHessian::SetCube( const ScalarWithHessian& a )
{
	assert ( GetNumberVariables()==a.GetNumberVariables() );
	*this = a;
	this->Cube();
}

inline void ScalarWithHessian::SetSqrt( const ScalarWithHessian& a )
{
	assert ( GetNumberVariables()==a.GetNumberVariables() );
	assert ( a.TheValue!=0.0 );
	*this = a;
	this->Sqrt();
}

inline void ScalarWithHessian::Negate()
{
	TheHessian.Negate();
	ThePartials.Negate();
	TheValue = -TheValue;
}

inline void ScalarWithHessian::Square()
{
	long num = ThePartials.GetLength();
	double twoValue = 2.0*TheValue;
	// New Hessian
	double* hPtr = TheHessian.GetPtr();
	double* pPtr1 = ThePartials.GetPtr();
	for ( long j=0; j<num; j++, pPtr1++ ) {
		double firstC = 2.0*(*pPtr1);
		double* pPtr2 = ThePartials.GetPtr();
		for ( long k=0; k<num; k++, pPtr2++ ) {
			*hPtr = firstC*(*pPtr2) + twoValue*(*hPtr);
			hPtr++;
		}
	}
	// New Gradient
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= twoValue;
	}
	// New Value
	TheValue = ::Square(TheValue);
}

inline void ScalarWithHessian::Cube()
{
	long num = ThePartials.GetLength();
	double sixValue = 6.0*TheValue;
	double threeValueSq = 3.0*::Square(TheValue);
	// New Hessian
	double* hPtr = TheHessian.GetPtr();
	double* pPtr1 = ThePartials.GetPtr();
	for ( long j=0; j<num; j++, pPtr1++ ) {
		double firstC = sixValue*(*pPtr1);
		double* pPtr2 = ThePartials.GetPtr();
		for ( long k=0; k<num; k++, pPtr2++ ) {
			*hPtr = firstC*(*pPtr2) + threeValueSq*(*hPtr);
			hPtr++;
		}
	}
	// New Gradient
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= threeValueSq;
	}
	// New Value
	TheValue = ::Square(TheValue);
}

inline void ScalarWithHessian::Sqrt()
{
	assert ( TheValue > 0.0 );
	long num = ThePartials.GetLength();
	double newValue = sqrt(TheValue);
	double halfSqrtInv = 0.5/TheValue;
	double secDenomInv = -halfSqrtInv*0.5/TheValue;
	// New Hessian
	double* hPtr = TheHessian.GetPtr();
	double* pPtr1 = ThePartials.GetPtr();
	for ( long j=0; j<num; j++, pPtr1++ ) {
		double firstC = secDenomInv*(*pPtr1);
		double* pPtr2 = ThePartials.GetPtr();
		for ( long k=0; k<num; k++, pPtr2++ ) {
			*hPtr = firstC*(*pPtr2) + halfSqrtInv*(*hPtr);
			hPtr++;
		}
	}
	// New Gradient
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= halfSqrtInv;
	}
	// New Value
	TheValue = newValue;
}

inline void ScalarWithHessian::Invert()
{
	assert ( TheValue > 0.0 );
	long num = ThePartials.GetLength();
	double newValue = 1.0/TheValue;
	double negDenomInv = -::Square(TheValue);
	double secDenomInv = -2.0*negDenomInv*newValue;
	// New Hessian
	double* hPtr = TheHessian.GetPtr();
	double* pPtr1 = ThePartials.GetPtr();
	for ( long j=0; j<num; j++, pPtr1++ ) {
		double firstC = secDenomInv*(*pPtr1);
		double* pPtr2 = ThePartials.GetPtr();
		for ( long k=0; k<num; k++, pPtr2++ ) {
			*hPtr = firstC*(*pPtr2) + negDenomInv*(*hPtr);
			hPtr++;
		}
	}
	// New Gradient
	double* pPtr = ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++ ) {
		*pPtr *= negDenomInv;
	}
	// New Value
	TheValue = newValue;

}

inline ScalarWithHessian& ScalarWithHessian::operator= ( double c )
{
	TheValue = c;
	ThePartials.SetZero();
	TheHessian.SetZero();
	return *this;
}


inline ScalarWithHessian& ScalarWithHessian::operator += ( double c )
{
	TheValue += c;
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::operator -= ( double c )
{
	TheValue -= c;
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::operator *= ( double c )
{
	TheHessian *= c;
	ThePartials *= c;
	TheValue *= c;
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::operator /= ( double c )
{
	assert ( c!=0.0 );
	double cInv = 1.0/c;
	TheHessian *= cInv;
	ThePartials *= cInv;
	TheValue *= cInv;
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::operator += ( const ScalarWithHessian& a )
{
	TheHessian += a.TheHessian;
	ThePartials += a.ThePartials;
	TheValue += a.TheValue;
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::operator -= ( const ScalarWithHessian& a )
{
	TheHessian -= a.TheHessian;
	ThePartials -= a.ThePartials;
	TheValue -= a.TheValue;
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::operator *= ( const ScalarWithHessian& a )
{
	long num = ThePartials.GetLength();
	// New Hessian
	double* hPtr1 = TheHessian.GetPtr();
	const double* hPtr2 = a.TheHessian.GetPtr();
	double* pPtr1x = ThePartials.GetPtr();
	const double* pPtr2x = a.ThePartials.GetPtr();
	for ( long j=0; j<num; j++, pPtr1x++, pPtr2x++ ) {
		double* pPtr1y = ThePartials.GetPtr();
		const double* pPtr2y = a.ThePartials.GetPtr();
		for ( long k=0; k<num; k++, pPtr1y++, pPtr2y++ ) {
			*hPtr1 = (*hPtr1)*a.TheValue + (*pPtr1x)*(*pPtr2y) + (*pPtr1y)*(*pPtr2x) + TheValue*(*hPtr2);
			hPtr1++;
			hPtr2++;
		}
	}
	// New Gradient
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr = (*pPtr)*a.TheValue + TheValue*(*pPtrA);
	}
	// New Value
	TheValue *= a.TheValue;
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::operator /= ( const ScalarWithHessian& a )
{
	assert ( a.TheValue!=0.0 && ::Cube(a.TheValue)!=0.0 );
	long num = ThePartials.GetLength();
	double aInv = 1.0/a.TheValue;
	double aInvSq = ::Square(aInv);
	double thisOveraSq = TheValue*aInvSq;
	double twoThisOveraCubed = 2.0*thisOveraSq*aInv;
	// New Hessian
	double* hPtr1 = TheHessian.GetPtr();
	const double* hPtr2 = a.TheHessian.GetPtr();
	double* pPtr1x = ThePartials.GetPtr();
	const double* pPtr2x = a.ThePartials.GetPtr();
	for ( long j=0; j<num; j++, pPtr1x++, pPtr2x++ ) {
		double* pPtr1y = ThePartials.GetPtr();
		const double* pPtr2y = a.ThePartials.GetPtr();
		for ( long k=0; k<num; k++, pPtr1y++, pPtr2y++ ) {
			*hPtr1 = (*hPtr1)*aInv - ((*pPtr1x)*(*pPtr2y)+(*pPtr1y)*(*pPtr2x))*aInvSq + twoThisOveraCubed*(*pPtr2x)*(*pPtr2y) - thisOveraSq*(*hPtr2);
			hPtr1++;
			hPtr2++;
		}
	}
	// New Gradient
	double* pPtr = ThePartials.GetPtr();
	const double* pPtrA = a.ThePartials.GetPtr();
	for ( long i=0; i<num; i++, pPtr++, pPtrA++ ) {
		*pPtr = (*pPtr)*aInv - thisOveraSq*(*pPtrA);
	}
	// New Value
	TheValue *= aInv;
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::SetScaled( const ScalarWithHessian& a, double factor )
{
	TheValue = factor*a.TheValue;
	ThePartials.SetScaled( a.ThePartials, factor );
	TheHessian.SetScaled( a.TheHessian, factor );
	return *this;
}

inline ScalarWithHessian& ScalarWithHessian::AddScaled( const ScalarWithHessian& a, double factor )
{
	TheValue += factor*a.TheValue;
	ThePartials.AddScaled( a.ThePartials, factor );
	TheHessian.AddScaled( a.TheHessian, factor );
	return *this;
}

#endif // SCALAR_WITH_HESSIAN_H
