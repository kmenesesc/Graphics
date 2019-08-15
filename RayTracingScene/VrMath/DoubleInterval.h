
// DoubleInterval
//
//  A simple interval arithmetic class for doubles
//
//  WARNING: THIS CODE HAS NEVER BEEN COMPILED OR TESTED AT ALL.
//  BUG ALERT!!!
//
// NEXT TODO: Add more methods for interactions with scalars
//

#ifndef DOUBLE_INTERVAL_H
#define DOUBLE_INTERVAL_H

#include <assert.h>
#include "MathMisc.h"

class DoubleInterval {

public:
	DoubleInterval( double min, double max );
	DoubleInterval( const DoubleInterval& d );

	DoubleInterval& SetScalar( double value );
	
	double GetMin() const { return minValue; }
	double GetMax() const { return maxValue; }

	bool IsInvalid() const { return isInvalidFlag; }

	// Note that <= and > are not opposites!  Nor are >= and <.
	bool operator>= ( const DoubleInterval& d ) const { return minValue>=d.maxValue; }
	bool operator> ( const DoubleInterval& d ) const { return minValue>d.maxValue; }
	bool operator<= ( const DoubleInterval& d ) const { return maxValue<=d.minValue; }
	bool operator< ( const DoubleInterval& d ) const { return maxValue<d.minValue; }

	bool operator>= ( double d ) const { return minValue>=d; }
	bool operator> ( double d ) const { return minValue>d; }
	bool operator<= ( double d ) const { return maxValue<=d; }
	bool operator< ( double d ) const { return maxValue<d; }

	DoubleInterval& operator+= ( const DoubleInterval& d );
	DoubleInterval& operator-= ( const DoubleInterval& d );
	DoubleInterval& operator*= ( const DoubleInterval& d );
	DoubleInterval& operator/= ( const DoubleInterval& d );
	DoubleInterval& operator+= ( double d );
	DoubleInterval& operator-= ( double d );
	DoubleInterval& operator*= ( double d );
	DoubleInterval& operator/= ( double d );
	DoubleInterval& DivideInto ( double d );

private:

	static const bool checkValid = true;
	bool isInvalidFlag;
	double minValue, maxValue;

};

// **************************************************
// Function Prototypes
// **************************************************

DoubleInterval operator+( const DoubleInterval& d, const DoubleInterval& e );
DoubleInterval operator-( const DoubleInterval& d, const DoubleInterval& e );
DoubleInterval operator*( const DoubleInterval& d, const DoubleInterval& e );
DoubleInterval operator/( const DoubleInterval& d, const DoubleInterval& e );

// **************************************************
// Inlined methods for DoubleInterval
// **************************************************

DoubleInterval::DoubleInterval( double min, double max ) {
	assert( min<=max );
	minValue = min;
	maxValue = max;
	isInvalidFlag = (min>max);
	if ( checkValid ) {
		assert(!isInvalidFlag);
	}
}

DoubleInterval::DoubleInterval( const DoubleInterval& d )
{
	minValue = d.minValue;
	maxValue = d.maxValue;
	isInvalidFlag = d.isInvalidFlag;
}

DoubleInterval& DoubleInterval::SetScalar( double value )
{
	minValue = value;
	maxValue = value;
	isInvalidFlag = false;
	return *this;
}

DoubleInterval& DoubleInterval::operator+= ( const DoubleInterval& d )
{
	minValue += d.minValue;
	maxValue += d.maxValue;
	isInvalidFlag = isInvalidFlag||d.isInvalidFlag;
	return *this;
}

DoubleInterval& DoubleInterval::operator+= ( double d )
{
	minValue += d;
	maxValue += d;
	return *this;
}

DoubleInterval& DoubleInterval::operator-= ( const DoubleInterval& d )
{
	minValue -= d.minValue;
	maxValue -= d.maxValue;
	isInvalidFlag = isInvalidFlag||d.isInvalidFlag;
	return *this;
}

DoubleInterval& DoubleInterval::operator-= ( double d )
{
	minValue -= d;
	maxValue -= d;
	return *this;
}

// The logic in operator/= could slightly improved on, but does work as is.
DoubleInterval& DoubleInterval::operator*= ( const DoubleInterval& d )
{
	double minmin = minValue*d.minValue;
	double minmax = minValue*d.maxValue;
	double maxmin = maxValue*d.minValue;
	double maxmax = maxValue*d.maxValue;
	maxValue = Min(Min(minmin,minmax),Min(maxmin,maxmax));
	maxValue = Max(Max(minmin,minmax),Max(maxmin,maxmax));
	isInvalidFlag = isInvalidFlag||d.isInvalidFlag;
	return *this;
}

DoubleInterval& DoubleInterval::operator*= ( double d )
{
	double min = minValue*d;
	double max = maxValue*d;
	if ( d>=0.0 ) {
		minValue = min;
		maxValue = max;
	}
	else {
		minValue = max;
		maxValue = min;
	}
	return *this;
}

// The logic in operator/= could substantially improved on, but does work as is.
DoubleInterval& DoubleInterval::operator/= ( const DoubleInterval& d )
{
	if ( d.minValue<=0.0 && d.maxValue>=0.0 ) {
		isInvalidFlag = true;
		if ( checkValid ) {
			assert(!isInvalidFlag);
		}
		return *this;
	}
	double minmin = minValue/d.minValue;
	double minmax = minValue/d.maxValue;
	double maxmin = maxValue/d.minValue;
	double maxmax = maxValue/d.maxValue;
	maxValue = Min(Min(minmin,minmax),Min(maxmin,maxmax));
	maxValue = Max(Max(minmin,minmax),Max(maxmin,maxmax));
	isInvalidFlag = isInvalidFlag||d.isInvalidFlag;
	return *this;
}

DoubleInterval& DoubleInterval::operator/= ( double d )
{
	if ( d==0 ) {
		isInvalidFlag = true;
		if ( checkValid ) {
			assert(!isInvalidFlag);
		}
		return *this;
	}
	double min = minValue/d;
	double max = maxValue/d;
	if ( d>0.0 ) {
		minValue = min;
		maxValue = max;
	}
	else {
		minValue = max;
		maxValue = min;
	}
	return *this;
}

DoubleInterval& DoubleInterval::DivideInto ( double d )
{
	if ( minValue<=0.0 && maxValue>=0.0 ) {
		isInvalidFlag = true;
		if ( checkValid ) {
			assert(!isInvalidFlag);
		}
		return *this;
	}
	double min = d/minValue;
	double max = d/maxValue;
	if ( min<=max ) {
		minValue = min;
		maxValue = max;
	}
	else {
		minValue = max;
		maxValue = min;
	}
	return *this;
}




#endif // DOUBLE_INTERVAL_H