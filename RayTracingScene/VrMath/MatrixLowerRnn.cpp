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
//
//

#include "MatrixLowerRnn.h"
#include "MatrixRmn.h"

// Constructor
MatrixLowerRnn::MatrixLowerRnn( const MatrixLowerRnn& init )
{
	N = init.N;
	AllocSize = (N*(N+1)) >> 1;
	x = new double[AllocSize];
	const double* fromPtr = init.x;
	double* toPtr = x;
	for ( long i = AllocSize; i>0; i-- ) {
		*(toPtr++) = *(fromPtr++);
	}
}

// Set equal to the lower triangular part of B (includes the diagonal)
void MatrixLowerRnn::Set( const MatrixRmn& B )
{
	assert(B.GetNumRows()==B.GetNumColumns());
	N = B.GetNumRows();
	SetSize(N);
	const double* fromPtr = B.GetPtr();
	long extraStep = 0;
	double* toPtr = x;
	for ( long i=N; i>0; i-- ) {
		for ( long j=i; j>0; j-- ) {
			(*(toPtr++)) = (*(fromPtr++));
		}
		extraStep++;
		fromPtr += extraStep;
	}
}


// result = (this)*(v)
void MatrixLowerRnn::Multiply( const VectorRn& v, VectorRn& result ) const
{
	assert(N==v.GetLength() && N==result.GetLength());

	const double* rowPtr = x;				// Pointer to the beginning of the current row
	double* resPtr = result.GetPtr();		// Pointer to current entry of result
	long i = N;
	for ( ; i>0; i-- ) {
		long rowStepSize = N;
		*resPtr = 0.0;
		const double* ePtr = rowPtr;		// Pointer to current entry in this row
		const double* vPtr = v.GetPtr();
		for ( long j=N-i; j>=0; j-- ) {
			*resPtr += (*vPtr)*(*ePtr);
			vPtr++;
			rowStepSize--;
			ePtr += rowStepSize;			// Next entry in this row
		}
		rowPtr++;						// Point to next row
		resPtr++;						// Point to next entry in result
	}
}

// result = (this^Transpose)*(v)
void MatrixLowerRnn::MultiplyTranspose( const VectorRn& v, VectorRn& result ) const
{
	assert(N==v.GetLength() && N==result.GetLength());
	
	const double* ePtr = x;			// Pointer to diagonal entry at top of current column
	double* resPtr = result.GetPtr();	// Pointer to current entry in result
	long diagStepSize = N;
	const double* vBasePtr = v.GetPtr();
	for ( ; diagStepSize>0; diagStepSize-- ) {
		*resPtr = 0.0;
		const double* vPtr = vBasePtr;
		for ( long j=diagStepSize; j>0; j-- ) {
			*resPtr += (*(vPtr++))*(*(ePtr++));
		}
		vBasePtr++;
		resPtr++;
	}
}

// result = (this)*(this^Transpose)*(v)
void MatrixLowerRnn::MultiplyLLT( const VectorRn& v, VectorRn& result ) const
{
	assert(N==v.GetLength() && N==result.GetLength());

	VectorRn& tempVector=VectorRn::GetWorkVector(N);	// Careful!
	MultiplyTranspose( v, tempVector );
	Multiply( tempVector, result );
}

// Solves the equation   (this))&(this^Transpose)*(result) = b; 
// Assumes diagonal entries are non-zero.
void MatrixLowerRnn::SolveLLT( const VectorRn& b, VectorRn* result ) const
{
	VectorRn& tempVector = VectorRn::GetWorkVector(b.GetLength());
	Solve( b, &tempVector );
	SolveTranspose( tempVector, result);
}

// Solves the equation   (this)*(result) = b; 
// Assumes diagonal entries are non-zero.
void MatrixLowerRnn::Solve( const VectorRn& b, VectorRn* result ) const
{
	assert(N==b.GetLength() && N==result->GetLength());

	const double* rowPtr = x;				// Pointer to the beginning of the current row
	const double* bPtr = b.GetPtr();		// Pointer to beginning of b
	for ( long i=N; i>0; i-- ) {
		const double* ePtr = rowPtr;		// Pointer to current entry in this row
		long rowStep = N;
		double* resPtr = result->GetPtr();
		double accum = 0.0;
		for ( long j=N-i; j>0; j-- ) {
			accum += (*resPtr)*(*ePtr);
			rowStep--;
			ePtr += rowStep;			// Next entry in this row
			resPtr++;
		}

		assert ( *ePtr!= 0.0 );				// Diagonal entry must not be zero (Entry number N-rowStepSize+1.)
		*resPtr = ((*bPtr)-accum)/(*ePtr);

		bPtr++;							// Next entry in b
		rowPtr++;						// Point to next row
	}
}

// Solves the equation   (this^Transpose)*(result) = b; 
// Assumes diagonal entries are non-zero.
void MatrixLowerRnn::SolveTranspose( const VectorRn& b, VectorRn* result ) const
{
	assert(N==b.GetLength() && N==result->GetLength());

	const double* colEndPtr = x + (((N*(N+1))>>1) - 1);		// Last entry in last column
	const double* bPtr = b.GetPtr() + (N - 1);				// Last entry in b.
	double* resEndPtr = result->GetPtr() + (N - 1);			// Last entry in result
	long colStep = 1;
	while ( colStep<=N ) {
		const double* ePtr = colEndPtr;
		double* resPtr = resEndPtr;
		double accum = 0.0;
		for ( long j=colStep; j>1; j-- ) {
			accum += (*(ePtr--))*(*(resPtr--));
		}

		assert( (*ePtr)!=0.0 );		// Diagonal entry is zero (entry number N+1-colStep).
		*resPtr = ((*bPtr)-accum)/(*ePtr);

		colEndPtr -= colStep;			// Last entry in next column (leftward)
		colStep++;
		bPtr--;						// Next entry in b
	}
}

// Perform a Givens transforms (plane rotation) to column colIdx and column N.
// "extraCol" is an N-vector that has entries that are to be thought of as
//		augmenting the matrix as an (N+1)-th column (column number N).
//		Entries in the extraCol that must be zero in order to preserve the lower
//		triangular property are assumed to actually be equal to zero (no checking)
// The Givens matrix multiplies on the right ("PostApply") and has the form
//		{ c    -s }
//      { s     c }  
// where zero entries are omitted.
void MatrixLowerRnn::PostApplyGivens( double c, double s, long colIdx, VectorRn* extraCol )
{
	assert ( N == extraCol->GetLength() );
	assert( fabs(c*c+s*s-1.0)<1.0e-12 );	// c*c+s*s should equal 1.0.

	double* colPtr = GetDiagPtr( colIdx );
	double* extraColPtr = extraCol->GetPtr( colIdx );
	
	for ( long i = N-colIdx; i>0; i-- ) {
		double temp = c*(*colPtr) + s*(*extraColPtr);
		*extraColPtr = c*(*extraColPtr) - s*(*colPtr);
		*colPtr = temp;
		colPtr++;
		extraColPtr++;
	}
}

// Rank One Update Lower Triangular: add on u*(v^Transpose)
// Change this matrix so that
//   (newthis)(newthis^Tranpose)  =  (update)*(update^Transpose),
// where
//	(update) = (this + alpa*v*(v^Transpose).
// Inputs are alpha, v.  
//   For this one, alpha can be positive, negative or zero.
// "this" is updated to become "newthis".
bool MatrixLowerRnn::UpdateCholeskyRankOne( double alpha, const VectorRn& v )
{
	assert ( v.GetLength()==N );

	if ( alpha==0 ) {
		return true;
	}
	else if ( alpha>0 ) {
		UpdateCholeskyPositiveRankOne( alpha, v );
		return true;
	}
	else {
		return UpdateCholeskyNegativeRankOne( alpha, v );
	}
}

// Same, but alpha must be positive
void MatrixLowerRnn::UpdateCholeskyPositiveRankOne( double alpha, const VectorRn& v )
{
	assert ( alpha>=0.0 );
	assert ( v.GetLength()==N );

	VectorRn& extraCol = VectorRn::GetWorkVector( N );
	extraCol.Set( v );
	extraCol *= sqrt(alpha);

	for ( long i=0; i<N; i++ ) {
		double a = *GetDiagPtr(i);
		double b = extraCol[i];
		double gamma = sqrt(a*a+b*b);
		double cos = a/gamma;
		double sin = b/gamma;
		PostApplyGivens( cos, sin, i, &extraCol );
	}
}

// Same, but alpha must be negative (or zero, but this is pointless).
bool MatrixLowerRnn::UpdateCholeskyNegativeRankOne( double alpha, const VectorRn& v )
{
	assert( alpha<=0.0 );
	assert ( v.GetLength()==N );

	bool retCode = true;

	VectorRn& pVec = VectorRn::GetWorkVector(N);
	Solve( v, &pVec );
	double normPsq = pVec.NormSq();
	double rhoSq = 1.0 + alpha*normPsq;
	if ( rhoSq<=1.0e-8 ) { 
		retCode = false;
		rhoSq = 1.0e-8;
		alpha = (rhoSq-1.0)/normPsq;	//Adjust alpha to keep positive definite
	}
	double q0 = sqrt(rhoSq);
	VectorRn& qVec = pVec;		// WorkVector of VectorRn (!)
	qVec *= sqrt(-alpha);
	VectorRn& extraGivensCol = VectorRn::GetWorkVector2(N);
	extraGivensCol.SetZero();
	for ( long i=N-1; i>=0; i-- ) {
		// Swap column i and the extra col.
		double qEntry = qVec[i];
		double gamma = sqrt( Square(q0) + Square(qEntry) );
		double sin = qEntry/gamma;
		double cos = q0/gamma;
		PostApplyGivens( cos, sin, i, &extraGivensCol );
		q0 = gamma;
	}

	return retCode;
}

#if 0
// THIS IS OLD AND USES OLD SPIKE CODE THAT NO LONGER EXISTS
// Rank One Update Lower Triangular: add on u*(v^Transpose)
// Change this matrix so that
//   (newthis)(newthis^Tranpose)  =  (update)*(update^Transpose),
// where
//	(update) = (this + u*(v^Transpose).
// Inputs are u and v.
// "this" is updated to become "newthis".
void MatrixLowerRnn::UpdateCholeskyWithRankTwo( const VectorRn& u, const VectorRn& v )
{
	// Needs to be re-written if it is going to be used.
	// This code is buggy -- written with a spike, not an extra column.
	assert(0);		

	assert ( N==u.GetLength() && N==v.GetLength() );

	VectorRn& spike = VectorRn::GetWorkVector(N-1);				// Temporary vector (reusable, but be careful!)
	spike.SetZero();

	// Phase 1: Apply Givens transformations (plane rotations) that transform
	//	vector v to a multiply of the standard unit vector e_n.
	const double* vPtr = v.GetPtr(N-1);
	double vLast = *vPtr;				// Last entry of v.
	for ( long i=N-1; i>=0; i-- ) {		// i is the column number
		vPtr--;
		double thisV = *vPtr;
		if ( thisV==0.0 ) {
			continue;
		}
		double c = vLast;
		double s = thisV;
		double norm = sqrt( c*c + s*s );
		if ( norm==0.0 ) {
			continue;			// Only happens with underflow!
		}
		c /= norm;
		s /= norm;
		vLast = c*vLast + s*thisV;
		PostApplyGivens( c, s, i, &spike );
	}
	assert ( fabs(fabs(vLast)-v.Norm()) < 1.0e-12*v.Norm() );	// DEBUG ONLY

	// Add u(v^Transpose)Q_1 to the spike.
	double* spikePtr = spike.GetPtr();
	const double* uPtr = u.GetPtr();
	for ( long i=0; i<N-1; i++ ) {
		*(spikePtr++) += (*uPtr++)*vLast;
	}

	// Phase 2: Apply Givens transformations to zero out the spike.
	spikePtr = spike.GetPtr();
	double* diagPtr = x;
	long diagStep = N;
	for ( long i=0; i<N-1; i++ ) {
		double c = *diagPtr;
		double s = -(*spikePtr);
		double norm = sqrt(c*c + s*s);
		if ( norm==0.0 ) {
			continue;		// This really should not happen, unless not positive definite
		}
		c /= norm;
		s /= norm;
		PostApplyGivens( c, s, i, &spike );
		
		diagPtr += diagStep;
		diagStep--;
		spikePtr++;
	}
}
#endif