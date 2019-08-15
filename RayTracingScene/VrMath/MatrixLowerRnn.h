/*
 *
 * Mathematics Subpackage (VrMath)
 *
 * Author: Samuel R. Buss
 *
 * Software is "as-is" and carries no warranty.  It may be used without
 *   restriction, but if you modify it, please change the filenames to
 *   prevent confusion between different versions.  Please acknowledge
 *   all use of the software in any publications or products based on it.
 *
 * Bug reports: Sam Buss, sbuss@ucsd.edu.
 *
 */

//
// MatrixLowerRmn:  Lower Triangular Square Matrix over reals  
//			(Variable dimensional matrices)
// Lower triangular matrix (with non-zero diagonal elements permitted)
//		Stored in column order
//

#ifndef MATRIX_LOWER_RNN_H
#define MATRIX_LOWER_RNN_H

#include <math.h>
#include <assert.h>
#include "LinearR3.h"
#include "VectorRn.h"
class MatrixRmn;

class MatrixLowerRnn 
{
public:
	MatrixLowerRnn();								// Null constructor
	MatrixLowerRnn( long dimension );				// Constructor with length
	MatrixLowerRnn( const MatrixLowerRnn& init );	// Set equal to another lower triangular matrix
	~MatrixLowerRnn();									// Destructor

	// Matrix size (must be square)
	void SetSize( long dimension );
	long GetDimension() const { return N; }
	
	// Accessors to entries in the matrix
	double* GetPtr() { return x; }
	const double* GetPtr() const { return x; }

	// Return pointer to i-th diagonal entry (i - zero based)
	double* GetDiagPtr( long i ) { return x + ((N*i) - ((i*(i-1))>>1)); }
	const double* GetDiagPtr( long i ) const { return x + ((N*i) - ((i*(i-1))>>1)); }
	double* GetLastPtr() { return x + (((N*(N+1))>>1) - 1); }
	const double* GetLastPtr() const { return x + (((N*(N+1))>>1) - 1); }

	// Set routines
	void SetZero(); 
	void SetIdentity();
	void SetDiagonalEntries( double d );
	void Set( const MatrixRmn& B );				// Set equal to the lower triangular part of B

	// Operations on VectorRn's
	void Multiply( const VectorRn& v, VectorRn& result ) const;				// result = (this)*(v)
	void MultiplyTranspose( const VectorRn& v, VectorRn& result ) const;	// Equivalent to mult by row vector on left
	void MultiplyLLT( const VectorRn& v, VectorRn& result ) const;			// result (this)*(this^transpose)*(v)

	// Solving systems of equations. (Easy since lower triangular)
	// Assumes diagonal entries are non-zero.
	void Solve( const VectorRn& b, VectorRn* x ) const;	  // Solves the equation   (this)*x = b; 
	void SolveTranspose( const VectorRn& b, VectorRn* x ) const;	  // Solves the equation   (this^Transpose)*x = b; 
	void SolveLLT( const VectorRn& b, VectorRn* x) const;		// Solves (this)*(this^Transpose)*x = b

	// Perform a Givens transforms (plane rotation) to column colIdx and column N.
	// "extraCol" is an (N-1) vector that has entries that are to be thought of as
	//		augmenting the matrix with entries in the final column
	//		Entries in the extraCol that must be zero are assumed to be equal to zero (no checking)
	void PostApplyGivens( double c, double s, long colIdx, VectorRn* extraCol );

	// Rank one updates.  Sets this = this + alpha*this*(this^transpose)
	bool UpdateCholeskyRankOne ( double alpha, const VectorRn& v );
	void UpdateCholeskyPositiveRankOne ( double alpha, const VectorRn& v );
	bool UpdateCholeskyNegativeRankOne ( double alpha, const VectorRn& v );

	// Rank Two Cholesky Update: add on u*(v^Transpose)
	// void UpdateCholeskyWithRankTwo( const VectorRn& u, const VectorRn& v );	// TO BE DEBUGGED

private:
	long N;					// Dimension (number of rows equals number of columns) 
	double *x;				// Array of vector entries - stored in column order
	long AllocSize;			// Allocated size of the x array

};

// Simple constructor
inline MatrixLowerRnn::MatrixLowerRnn()
{
	N = 0;
	AllocSize = 0;
	x = 0;
}

inline MatrixLowerRnn::MatrixLowerRnn( long dimension )
{
	assert ( dimension>0 );
	AllocSize = 0;
	x = 0;
	SetSize( dimension );
	N = dimension;
}

// Destructor

inline MatrixLowerRnn::~MatrixLowerRnn()
{
	delete[] x;
}

// Resizing and re-allocating
// No copying of old data is done
inline void MatrixLowerRnn::SetSize( long dimension )
{
	assert ( dimension>=0 );
	N = dimension;
	long neededSize = (N*(N+1)) >> 1;
	if ( neededSize > AllocSize ) {
		delete[] x;
		AllocSize = neededSize;
		x = new double[AllocSize];
	}
}

inline void MatrixLowerRnn::SetZero()
{
	double* toPtr = x;
	for ( long i=(N*(N+1))>>1; i>0; i-- ) {
		*(toPtr++) = 0.0;
	}
}

inline void MatrixLowerRnn::SetIdentity()
{
	SetZero();
	SetDiagonalEntries(1.0);
}

inline void MatrixLowerRnn::SetDiagonalEntries( double d )
{
	double *toPtr = x;
	long toStep = N;
	for ( long i=N; i>0; i-- ) {
		*toPtr = d;
		toPtr += toStep;
		toStep--;
	}
}




#endif // MATRIX_LOWER_RNN_H

