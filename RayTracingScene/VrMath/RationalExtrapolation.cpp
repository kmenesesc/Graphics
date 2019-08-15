//
// Rational Extrapolation
//
// Author: Sam Buss.  All rights and ownership retained. 2005.
//

#include "RationalExtrapolation.h"
#include "assert.h"

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
double RationalExtrapolate( int n, double *hValues, double *funcValues, int* returnCode, double* retAccuracy )
{
	RationalExtrapData ratData(n);

	int i;
	for ( i=0; i<n; i++ ) {
		ratData.AddValue( hValues[i], funcValues[i] );
	}

	*returnCode = ratData.GetOrder();
	*retAccuracy = ratData.GetEstimateDelta();
	return ratData.GetBestEstimate();
}

// Add one more value to the extrapolator.
// Subsequent values should be decreasing in absolute value.
double RationalExtrapData::AddValue( double hValue, double funcValue )
{
	assert ( hValue!=0.0 );

	// Step 1: Update the extrapolated values at all orders
	int k;
	double interResult = funcValue;		// temporary holder of an extrapolation value
	double* oldHvaluePtr = Hvalues+MaxOrder;
	for ( k=0; k<=OrderOfAccuracy; k++ ) {
		double denom = (*(oldHvaluePtr--))-hValue;
		assert ( denom!=0.0 );		// hValues not permitted to be repeated!
		double next = interResult + hValue*(interResult-ExtrapolationValues[k])/denom;
		DeltaValues[k] = fabs(interResult-ExtrapolationValues[k]);
		ExtrapolationValues[k] = interResult;
		interResult = next;
	}
	if ( OrderOfAccuracy<MaxOrder ) {
		// k is equal to OrderOfAccuracy+1
		DeltaValues[k] = fabs(interResult-ExtrapolationValues[k]);
		ExtrapolationValues[k] = interResult;
		OrderOfAccuracy++;
	}

	// Step 2: Find the order that gives the apparently best estimate.
	int bestK = 0;
	double bestDelta = DeltaValues[0];
	for ( k=1; k<=OrderOfAccuracy; k++ ) {
		if ( DeltaValues[k] < bestDelta ) {
			bestDelta = DeltaValues[k];
			bestK = k;
		}
	}
	if ( bestK==OrderOfAccuracy-1 && bestK<MaxOrder ) {
		// Optimistically go for the highest order accuracy if it seems warranted.
		bestK++;
	}
	OrderOfBestAccuracy = bestK;

	if ( OrderOfAccuracy<MaxOrder ) {
		// Copy last h value for next iteration to a higher order
		ExtrapolationValues[OrderOfAccuracy+1] = ExtrapolationValues[OrderOfAccuracy];
	}

	// Slide down h values
	for ( k=MaxOrder-OrderOfAccuracy; k<MaxOrder; k++ ) {
		Hvalues[k] = Hvalues[k+1];
	}
	Hvalues[MaxOrder] = hValue;

	return GetBestEstimate();
}
