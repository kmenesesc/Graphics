//
// Integration.cpp - Author: Sam Buss
//

#include <stdio.h> // For debugging
#include <assert.h>
#include <math.h>
#include "Integration.h"
#include "MathMisc.h"
#include "RationalExtrapolation.h"

// ****************************************************************
// Integrates using the trapezoidal rule.
//    Iterations: if n=number of Iterations, then 2^{n+1}+1 points are sampled
//				at equally spaced intervals.
// ****************************************************************
double IntegrateTrapezoidRule( double a, double b, double f(double),  
							 int minIterations, int maxIterations )
{
	assert ( minIterations<=maxIterations );

	double h = b-a;
	double sumValues = (f(a)+f(b))*0.5;
	int numNewValues = 1;
	double prevResult = sumValues*h;
	double prevResultAbs = fabs(prevResult);
	double newResult;
	int stopIter = maxIterations+1;
	int iterNum;
	for( iterNum=1; iterNum<stopIter; iterNum++ ) {
		h *= 0.5;
		sumValues += f(a+h);
		for ( int i = 1; i<numNewValues; i++ ) {
			sumValues += f(a+(double)(1+(i<<1))*h);
		}
		numNewValues = numNewValues<<1;			// Number of values for next iteration
		newResult = sumValues*h;
		double newResultAbs = fabs(newResult);
		if ( iterNum>=minIterations ) {
			double diff = fabs(newResult-prevResult);
			if ( diff<1.0e-15*prevResultAbs ) {
				stopIter = iterNum;			// Exit now
			}
			else if ( diff<1.0e-12*prevResultAbs
				|| ( newResultAbs<1.0e-13 && prevResultAbs<1.0e-13 ) )
			{
				stopIter = Min(iterNum+3, stopIter);	// Exit after two more iterations
			}
		}
		prevResult = newResult;
		prevResultAbs = newResultAbs;
	}

	fprintf( stdout, "Trapezoidal Rule: %d iterations.\n", iterNum-1 );

	return newResult;
}

// ****************************************************************
// Trapezoidal rule with polynomial extrapolation (Romberg extrapolation)
// Performs integration with the Trapezoidal rule, using Romberg step-size
//		with polynomial extrapolation to speed convergence.
//    Iterations: if n=number of Iterations, then 2^{n+1}+1 points are sampled
//				at equally spaced intervals.
// ****************************************************************
double IntegrateTrapezoidExtrap( double a, double b, double f(double),
							  int minIterations, int maxIterations, 
							  int numRombergLevels )
{
	assert ( minIterations<=maxIterations );

	RationalExtrapData ratData( numRombergLevels );

	double h = b-a;
	double sumValues = (f(a)+f(b))*0.5;
	int numNewValues = 1;
	double prevResult = sumValues*h;
	ratData.AddValue( h*h, prevResult );
	double prevResultAbs = fabs(prevResult);
	double newResult;
	int stopIter = maxIterations+1;
	int iterNum;
	for( iterNum=1; iterNum<stopIter; iterNum++ ) {
		// Step 1: Calculate new trapezoidal estimate as "interResult"
		h *= 0.5;
		sumValues += f(a+h);
		for ( int i = 1; i<numNewValues; i++ ) {
			sumValues += f(a+(double)(1+(i<<1))*h);
		}
		numNewValues = numNewValues<<1;			// Number of values for next iteration

		// Step 2: Update the Romberg values for extrapolation
		ratData.AddValue( h*h, sumValues*h );
		newResult = ratData.GetBestEstimate();
		double newResultAbs = fabs(newResult);
		double diff = ratData.GetEstimateDelta();

		fprintf( stdout, "TrapExtrap: Iteration %2d: %.16lg (Level %d).\n", iterNum, newResult, ratData.GetOrder());	// DEBUG

		// Step 3: Termination conditions
		if ( iterNum>=minIterations ) {
			if ( diff<1.0e-15*newResultAbs ) {
				stopIter = iterNum;			// Exit now
			}
			else if ( diff<1.0e-12*prevResultAbs
				|| ( newResultAbs<1.0e-13 && prevResultAbs<1.0e-13 ) )
			{
				stopIter = Min(iterNum+2, stopIter);	// Exit after one more iteration
			}
		}
		prevResult = newResult;
		prevResultAbs = newResultAbs;
	}

	fprintf( stdout, "TrapezoidExtrap: %d iterations.\n", iterNum-1 );
	
	return newResult;
}

// ****************************************************************
// Trapezoidal rule with rational extrapolation (Romberg extrapolation)
// Performs integration with the Trapezoidal rule, using Romberg step-size
//		with rational extrapolation to speed convergence.
//    Iterations: if n=number of Iterations, then 2^{n+1}+1 points are sampled
//				at equally spaced intervals.
// ****************************************************************
double IntegrateTrapezoidRatExtrap( double a, double b, double f(double),
							  int minIterations, int maxIterations, int numRombergLevels )
{
	assert ( minIterations<=maxIterations );

	const int MaxNumRombergLevels = 12;
	assert ( numRombergLevels<=MaxNumRombergLevels );

	double T[MaxNumRombergLevels];			// Table of values for Romberg interpolation
	double deltaT[MaxNumRombergLevels];		// Changes in Romberg values

	double h = b-a;
	double sumValues = (f(a)+f(b))*0.5;
	int numNewValues = 1;
	double prevResult = sumValues*h;
	T[0] = prevResult;
	int prevNumTs = 0;				// So far, only T[0] is stored.
	double prevResultAbs = fabs(prevResult);
	double newResult;
	int stopIter = maxIterations+1;
	int iterNum;
	for( iterNum=1; iterNum<stopIter; iterNum++ ) {
		// Step 1: Calculate new trapezoidal estimate as "interResult"
		h *= 0.5;
		sumValues += f(a+h);
		for ( int i = 1; i<numNewValues; i++ ) {
			sumValues += f(a+(double)(1+(i<<1))*h);
		}
		numNewValues = numNewValues<<1;			// Number of values for next iteration

		// Step 2: Update the Romberg values for extrapolation
		// interResult - intermediate results, to be stored into T[k].
		double interResult = sumValues*h;		// Result of Trapezoid rule
		// T[0], ... , T[prevNumTs] are presently set.
		// T[0], ... , T[numNewTs] will be used
		int numNewTs = Min(numRombergLevels-1, prevNumTs+1);	
		double fourPowerK = 1.0;
		double prevTk = 0.0;	// T[-1] for rational extrapolation
		for ( int k=0; k<numNewTs; k++ ) {
			fourPowerK *= 4.0;
			// temp equals T[k+1].  prevTk equals old T[k-1] value.
			double denom2 = interResult-prevTk;
			double denom = fourPowerK*(T[k]-prevTk)/denom2 - 1.0;
			if ( denom2==0.0 || denom==0.0 ) {
				// Rational extrapolation failed (degenerate positions).  Stop extrapolation here.
				numNewTs = k;
				prevNumTs = k;
				break;
			}
			double temp = interResult + (interResult-T[k])/denom;
			prevTk = T[k];
			deltaT[k] = interResult - T[k];
			T[k] = interResult;
			interResult = temp;
		}
		deltaT[numNewTs] = interResult - T[numNewTs];
		T[numNewTs] = interResult;
		if ( numNewTs<numRombergLevels-1 ) {
			T[numNewTs+1] = T[numNewTs];	// Duplicate the top entry if there is room
		}
		double diff = fabs(deltaT[numNewTs]);
		int bestM = numNewTs;
		for ( int m = numNewTs-1; m>=0; m-- ) {
			double diffB = fabs(deltaT[m]);
			if ( !(diffB>=diff) ) {
				diff = diffB;
				bestM = m;
			}
		}
		newResult = T[bestM];
		prevNumTs = bestM;	// Use T[0],...,T[prevNumTs] in next iteration.
		double newResultAbs = fabs(newResult);

		fprintf( stdout, "TrapRatExtrap: Iteration %2d: %.16lg (Level %d).\n",iterNum,newResult,bestM);	// DEBUG

		// Step 3: Termination conditions
		if ( iterNum>=minIterations ) {
			if ( diff<1.0e-15*prevResultAbs ) {
				stopIter = iterNum;			// Exit now
			}
			else if ( diff<1.0e-12*prevResultAbs
				|| ( newResultAbs<1.0e-13 && prevResultAbs<1.0e-13 ) )
			{
				stopIter = Min(iterNum+2, stopIter);	// Exit after one more iteration
			}
		}
		prevResult = newResult;
		prevResultAbs = newResultAbs;
	}

	fprintf( stdout, "TrapezoidRatExtrap: %d iterations.\n", iterNum-1 );

	return newResult;
}

// Integrates using Simpson's rule.
//    Iterations: if n=number of Iterations, then 2^{n+1}+1 points are sampled
//				at equally spaced intervals.
double IntegrateSimpsonRule( double a, double b, double f(double),  
							 int minIterations, int maxIterations )
{
	assert ( minIterations<=maxIterations );

	double h = 0.5*(b-a);
	double sumEndPts = f(a) + 2.0*f(a+h) + f(b);
	int numNewValues = 2;
	double prevResult = sumEndPts*h*OneThird;
	double prevResultAbs = fabs(prevResult);
	double newResult;
	int stopIter = maxIterations+1;
	int iterNum;
	for( iterNum=2; iterNum<stopIter; iterNum++ ) {
		h *= 0.5;
		double newValuesSum = f(a+h);
		for ( int i = 1; i<numNewValues; i++ ) {
			newValuesSum += f(a+(double)(1+(i<<1))*h);
		}
		newValuesSum *= 2.0;
		numNewValues = numNewValues<<1;			// Number of values for next iteration
		sumEndPts += newValuesSum;
		newResult = (sumEndPts+newValuesSum)*h*OneThird;
		double newResultAbs = fabs(newResult);

		// fprintf(stdout,"Simpson iteration %2d: %.16lg.\n", iterNum, newResult );

		if ( iterNum>=minIterations ) {
			double diff = fabs(newResult-prevResult);
			if ( diff<1.0e-15*prevResultAbs ) {
				stopIter = iterNum;				// Exit now
			}
			else if ( diff<1.0e-12*prevResultAbs
				|| ( newResultAbs<1.0e-13 && prevResultAbs<1.0e-13 ) )
			{
				stopIter = Min(iterNum+3, stopIter);
			}
		}
		prevResult = newResult;
		prevResultAbs = newResultAbs;
	}

	fprintf(stderr, "Simpson Rule: %d iterations.\n", iterNum-1 );

	return newResult;
}

// Performs integration with Milne's rule.
//    Iterations: if n=number of Iterations, then 2^n + 1 points are sampled
//				at equally spaced intervals.
double IntegrateMilneRule( double a, double b, double f(double),  
						   int minIterations, int maxIterations )
{
	double baseSumValues = f(a)+f(b);

	double prevValuesSum = f( (a+b)*0.5 );
	double h = (b-a)/4.0;
	int numNewValues = 2;
	double newValuesSum = f((3.0*a+b)*0.25)+f((a+3.0*b)*0.25);
	double prevResult = (7.0*baseSumValues+12*prevValuesSum+32.0*newValuesSum)*h*(4.0/90.0);
	double prevResultAbs = fabs(prevResult);
	double newResult;
	// fprintf(stdout, "Milne iteration  2: %.16lg.\n", prevResult );

	int iterNum;
	int stopIter = maxIterations+1;
	for ( iterNum=3; iterNum<stopIter; iterNum++ ) {
		baseSumValues += 2.0*prevValuesSum;
		prevValuesSum = newValuesSum;
		h *= 0.5;
		numNewValues *= 2;
		newValuesSum = 0.0;
		for ( int i = 0; i<numNewValues; i++ ) {
			double c = a + (double)(1+2*i)*h;
			newValuesSum += f(c);
		}
		newResult = (7.0*baseSumValues+12*prevValuesSum+32.0*newValuesSum)*h*(4.0/90.0);
		double newResultAbs = fabs(newResult);

		// fprintf(stdout, "Milne iteration %2d: %.16lg.\n", iterNum, newResult );

		if ( iterNum>=minIterations ) {
			double diff = fabs(newResult-prevResult);
			if ( diff<1.0e-15*prevResultAbs )
			{
				stopIter = iterNum;		// Exit now
			}
			else if ( diff<1.0e-12*prevResultAbs
				|| ( newResultAbs<1.0e-13 && prevResultAbs<1.0e-13 ) )
			{
				stopIter = Min(iterNum+2, stopIter);
			}
		}
		prevResult = newResult;
		prevResultAbs = newResultAbs;
	}

	fprintf(stderr, "Milne Iterations = %d.\n", iterNum-1 );

	return newResult;
}

