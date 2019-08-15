//
// Integration.h - Author: Sam Buss
//

#ifndef INTEGRATION_H
#define INTEGRATION_H

// Trapezoidal rule
double IntegrateTrapezoidRule( double a, double b, double f(double),
							   int minIterations = 2, int maxIterations = 15 );

// Trapezoidal rule with polynomial extrapolation (Romberg extrapolation)
double IntegrateTrapezoidExtrap( double a, double b, double f(double),
							  int minIterations = 2, int maxIterations = 12,
							  int numRombergLevels = 6);

// Trapezoidal rule with rational extrapolation (Romberg extrapolation)
double IntegrateTrapezoidRatExtrap( double a, double b, double f(double),
							  int minIterations = 2, int maxIterations = 12,
							  int numRombergLevels = 6);


// Simpson's rule
double IntegrateSimpsonRule( double a, double b, double f(double),  
							 int minIterations = 2, int maxIterations = 13 );

// double Milne's rule
double IntegrateMilneRule( double a, double b, double f(double),  
							 int minIters = 2, int maxIters = 9 );

#endif // INTEGRATION_H
