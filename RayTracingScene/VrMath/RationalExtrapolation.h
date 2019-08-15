//
// Rational Extrapolation
//
// Author: Sam Buss.  All rights and ownership retained. 2005.
//

#ifndef RATIONAL_EXTRAPOLATION_H
#define  RATIONAL_EXTRAPOLATION_H

#include <assert.h>
#include "MathMisc.h"

class RationalExtrapData;		// Holds data for a rational extrapolator

// Rational extrapolation (or, possibly interpolation)
// Intended usage: an extrapolator to help with numerical integration
// The input is:
//    an array of h values, hValues[0],..,hValue[n-1] and
//    an array of function values, funcValue[i] = f(hValue[i]
// They are extrapolated to given as estimated value for f(0).
// This estimated value is returned.
// The returned value "retAccuracy" estimates roughly the
//		the accuracy of the answer.  It equals the difference in
//		accuracy between using n-1 of the hValues versus using all n of them.
// The returned value returnCode is in the range [0,...,n-1] and is the
//		order of the accuracy. In effect, it means that the first (n-returnCode-1)
//		values were ignored.  E.g. order zero is the last value is taken as the best.
//		order one means that the last two values were linearly extrapolated.
// Note that the h values should be decreasing in absolute value in most cases.
//  Furthermore, the implicit assumption is that h's later in the array are the
//	more relevant.
double RationalExtrapolate( int n, double *hValues, double *funcValues, int* returnCode, double* retAccuracy );

// The use of the RationalExtrapData class is for more advanced usage.  
// The extra feature is that is allows you to maintain the old data
//	  for incremental updating of the extrapolated value.
// Usage: The standard usage is as follows:
//	1. Call RationalExtrapData(order.
//		Use order = 1 for linear extrapolation, order = 2 for quadratic, etc.
//		Or, order = (# of values)-1 for maximum order
//  2. Call AddValue repeatedly.  hValues decreasing in absolute value
//  3. Call GetBestEstimate() to get the best extrapolated value
//		and GetEstimateDelta to get an estimate of the error in the extrapolated value.

// For even functions, call AddValue with hValue squared (instead of hValue).

class RationalExtrapData {
public:

	RationalExtrapData( int maxOrder );
	~RationalExtrapData();

	void Reset();

	// Add one more value to the extrapolator.
	// Subsequent values should be decreasing in absolute value.
	double AddValue( double hValue, double funcValue );

	// Get the best estimate for the extrapolation
	double GetBestEstimate() const { return GetExtrapolatedValue(OrderOfBestAccuracy); }
	double GetEstimateDelta() const { return GetDeltaExtrapValue(OrderOfBestAccuracy); }

	// Get the highest order estimate for the extrapolation
	double GetHighestOrderEstimate () const { return GetExtrapolatedValue(OrderOfAccuracy); }
	double GetHighestEstimateDelta () const { return GetDeltaExtrapValue(OrderOfAccuracy); }

	// Specialized function for handling multiple orders
	double GetExtrapolatedValue( int order ) const;
	double GetDeltaExtrapValue( int order ) const;
	int GetOrder() const { return OrderOfAccuracy; };
	int GetBestOrder() const { return OrderOfBestAccuracy; };
	int GetMaxOrder() const { return MaxOrder; }

protected:
	int OrderOfAccuracy;
	int OrderOfBestAccuracy;
	int MaxOrder;

	double* Hvalues;
	double* ExtrapolationValues;
	double* DeltaValues;

};

inline RationalExtrapData::RationalExtrapData( int maxOrder ) {
	assert ( maxOrder>0 );
	MaxOrder = maxOrder;
	Hvalues = new double[maxOrder+1];
	ExtrapolationValues = new double[maxOrder+1];
	DeltaValues = new double[maxOrder+1];

	Reset();
}

inline RationalExtrapData::~RationalExtrapData()
{
	delete[] Hvalues;
	delete[] ExtrapolationValues;
	delete[] DeltaValues;
}

inline void RationalExtrapData::Reset()
{
	OrderOfAccuracy = -1;
	ExtrapolationValues[0] = DBL_NAN;
}

// Get the value at the i-th order of extrapolation
inline double RationalExtrapData::GetExtrapolatedValue( int order ) const
{
	assert( order>=0 && order<=OrderOfAccuracy);
	return ExtrapolationValues[order];
}

// Get the change in extrapolated value at the i-th order of extrapolation
inline double RationalExtrapData::GetDeltaExtrapValue ( int order ) const
{
	assert( order>=0 && order<=OrderOfAccuracy);
	return DeltaValues[order];
}



#endif  // RATIONAL_EXTRAPOLATION_H

