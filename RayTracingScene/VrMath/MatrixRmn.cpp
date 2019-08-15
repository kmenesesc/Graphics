
/*
 *
 * RayTrace Software Package, release 3.1.  December 20, 2006.
 *
 * Mathematics Subpackage (VrMath)
 *
 * Author: Samuel R. Buss
 *
 * Software accompanying the book
 *		3D Computer Graphics: A Mathematical Introduction with OpenGL,
 *		by S. Buss, Cambridge University Press, 2003.
 *
 * Software is "as-is" and carries no warranty.  It may be used without
 *   restriction, but if you modify it, please change the filenames to
 *   prevent confusion between different versions.  Please acknowledge
 *   all use of the software in any publications or products based on it.
 *
 * Bug reports: Sam Buss, sbuss@ucsd.edu.
 * Web page: http://math.ucsd.edu/~sbuss/MathCG
 *
 */

//
// MatrixRmn.cpp:  Matrix over reals  (Variable dimensional matrices)
//
//    Not very sophisticated yet.  Needs more functionality
//		To do: better handling of resizing.
//

#include "MatrixRmn.h"
#include "MatrixLowerRnn.h"
#include <float.h>

MatrixRmn MatrixRmn::WorkMatrix;		// Temporary work matrix
MatrixRmn MatrixRmn::WorkMatrix2;		// Temporary work matrix
MatrixRmn MatrixRmn::WorkMatrix3;		// Temporary work matrix

// Set equal to a lower triangular matrix
void MatrixRmn::Set ( const MatrixLowerRnn& B )
{
	long N = B.GetDimension();
	SetSize(N,N);
	const double* fromPtr = B.GetPtr();
	long extraStep = 0;
	double* toPtr = x;
	for ( long i=N; i>0 ; i-- ) {
		for ( long k=extraStep; k>0; k-- ) {
			*(toPtr++) = 0.0;
		}
		extraStep++;
		for ( long j=i; j>0; j-- ) {
			(*(toPtr++)) = (*(fromPtr++));
		}
	}
}

// Set the entries of the matrix in column order
// Assumes size already set correctly!!
void MatrixRmn::SetByColumns( const double* values )
{
	double* to = GetPtr();			// Point to beginning of matrix values
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		*(to++) = *(values++);
	}
}

// Set d equal to the i-th column of this matrix.
void MatrixRmn::GetColumn(long i, VectorRn* d ) const
{
	d->SetLength(NumRows);
	double* to = d->x;
	const double* from = x+i*NumRows;
	for ( i=NumRows; i>0; i-- ) {
		*(to++) = *(from++);
	}
}


// Set the entries of the matrix in row order
// Assumes size already set correctly!!
void MatrixRmn::SetByRows( const double* values )
{
	double* rowPtr = GetPtr();			// Point to beginning of row
	for ( long i=NumRows; i>0; i-- ) {
		double* to = rowPtr;
		for ( long j=NumCols; j>0; j-- ) {
			*to = *(values++);
			to += NumRows;
		}
        rowPtr++;
	}
}

// Set equal to transpose of M.
void MatrixRmn::SetTranspose( const MatrixRmn& M )
{
	long numRows = M.GetNumColumns();
	long numCols= M.GetNumRows();
	SetSize( numRows, numCols );
	const double* fromRow = M.x;
	double* toPtr = x;
	for ( long i=0; i<numCols; i++ ) {	// Loop over each new column
		const double* fromPtr = fromRow;
		for ( long j=0; j<numRows; j++ ) {
			*(toPtr++) = *fromPtr;
			fromPtr += numCols;
		}
		fromRow++;
	}
}



// Fill the diagonal entries with the value d.  The rest of the matrix is unchanged.
void MatrixRmn::SetDiagonalEntries( double d )
{
	long diagLen = Min( NumRows, NumCols );
	double* dPtr = x;
	for ( ; diagLen>0; diagLen-- ) {
		*dPtr = d;
		dPtr += NumRows+1;
	}
}

// Fill the diagonal entries with values in vector d.  The rest of the matrix is unchanged.
void MatrixRmn::SetDiagonalEntries( const VectorRn& d )
{
	long diagLen = Min( NumRows, NumCols );
	assert ( d.Length == diagLen );
	double* dPtr = x;
	double* from = d.x;
	for ( ; diagLen>0; diagLen-- ) {
		*dPtr = *(from++);
		dPtr += NumRows+1;
	}
}

// Fill the superdiagonal entries with the value d.  The rest of the matrix is unchanged.
void MatrixRmn::SetSuperDiagonalEntries( double d )
{
	long sDiagLen = Min( NumRows, NumCols-1 );
	double* to = x + NumRows;
	for ( ; sDiagLen>0; sDiagLen-- ) {
		*to = d;
		to += NumRows+1;
	}
}

// Fill the superdiagonal entries with values in vector d.  The rest of the matrix is unchanged.
void MatrixRmn::SetSuperDiagonalEntries( const VectorRn& d )
{
	long sDiagLen = Min( NumRows-1, NumCols );
	assert ( sDiagLen == d.Length );
	double* to = x + NumRows;
	double* from = d.x;
	for ( ; sDiagLen>0; sDiagLen-- ) {
		*to = *(from++);
		to += NumRows+1;
	}
}

// Fill the subdiagonal entries with the value d.  The rest of the matrix is unchanged.
void MatrixRmn::SetSubDiagonalEntries( double d )
{
	long sDiagLen = Min( NumRows, NumCols ) - 1;
	double* to = x + 1;
	for ( ; sDiagLen>0; sDiagLen-- ) {
		*to = d;
		to += NumRows+1;
	}
}

// Fill the subdiagonal entries with values in vector d.  The rest of the matrix is unchanged.
void MatrixRmn::SetSubDiagonalEntries( const VectorRn& d )
{
	long sDiagLen = Min( NumRows, NumCols ) - 1;
	assert ( sDiagLen == d.Length );
	double* to = x + 1;
	double* from = d.x;
	for ( ; sDiagLen>0; sDiagLen-- ) {
		*to = *(from++);
		to += NumRows+1;
	}
}

// Set the i-th column equal to d.
void MatrixRmn::SetColumn(long i, const VectorRn& d )
{
	assert ( NumRows==d.GetLength() && 0<=i && i<NumCols );
	double* to = x+i*NumRows;
	const double* from = d.x;
	for ( i=NumRows; i>0; i-- ) {
		*(to++) = *(from++);
	}
}

// Set the i-th column equal to d.
void MatrixRmn::SetRow(long i, const VectorRn& d )
{
	assert ( NumCols==d.GetLength() && 0<=i && i<NumRows );
	double* to = x+i;
	const double* from = d.x;
	for ( i=NumCols; i>0; i-- ) {
		*to = *(from++);
		to += NumRows;
	}
}

// Sets a "linear" portion of the array with the values from a vector d
// The first row and column position are given by startRow, startCol.
// Successive positions are found by using the deltaRow, deltaCol values
//	to increment the row and column indices.  There is no wrapping around.
void MatrixRmn::SetSequence( const VectorRn& d, long startRow, long startCol, long deltaRow, long deltaCol )
{
	long length = d.Length;
	assert( startRow>=0 && startRow<NumRows && startCol>=0 && startCol<NumCols );
	assert( startRow+(length-1)*deltaRow>=0 && startRow+(length-1)*deltaRow<NumRows );
 	assert( startCol+(length-1)*deltaCol>=0 && startCol+(length-1)*deltaCol<NumCols );
	double *to = x + startRow + NumRows*startCol;
	double *from = d.x;
	long stride = deltaRow + NumRows*deltaCol;
	for ( ; length>0; length-- ) {
		*to = *(from++);
		to += stride;
	}
}

// The matrix A is loaded, in into "this" matrix, based at (0,0).
//  The size of "this" matrix must be large enough to accomodate A.
//	The rest of "this" matrix is left unchanged.  It is not filled with zeroes!

void MatrixRmn::LoadAsSubmatrix( const MatrixRmn& A )
{
	assert( A.NumRows<=NumRows && A.NumCols<=NumCols );
	long extraColStep = NumRows - A.NumRows;
	double *to = x;
	double *from = A.x;
	for ( long i=A.NumCols; i>0; i-- ) {			// Copy columns of A, one per time thru loop
		for ( long j=A.NumRows; j>0; j-- ) {		// Copy all elements of this column of A
			*(to++) = *(from++);
		}
		to += extraColStep;
	}	
}

// The matrix A is loaded, in into "this" matrix, based at (i,j).
//  The size of "this" matrix must be large enough to accomodate A.
//	The rest of "this" matrix is left unchanged.  It is not filled with zeroes!

void MatrixRmn::LoadAsSubmatrix( long i, long j, const MatrixRmn& A )
{
	assert( A.NumRows+i<=NumRows && A.NumCols+j<=NumCols );
	long extraColStep = NumRows - A.NumRows;
	double *to = GetPtr(i,j);
	double *from = A.x;
	for ( long i=A.NumCols; i>0; i-- ) {			// Copy columns of A, one per time thru loop
		for ( long j=A.NumRows; j>0; j-- ) {		// Copy all elements of this column of A
			*(to++) = *(from++);
		}
		to += extraColStep;
	}	
}

// The matrix A is loaded, in transposed order into "this" matrix, based at (0,0).
//  The size of "this" matrix must be large enough to accomodate A.
//	The rest of "this" matrix is left unchanged.  It is not filled with zeroes!
void MatrixRmn::LoadAsSubmatrixTranspose( const MatrixRmn& A )
{
	assert( A.NumRows<=NumCols && A.NumCols<=NumRows );
	double* rowPtr = x;
	double* from = A.x;
	for ( long i=A.NumCols; i>0; i-- ) {				// Copy columns of A, once per loop
		double* to = rowPtr;
		for ( long j=A.NumRows; j>0; j-- ) {			// Loop copying values from the column of A
			*to = *(from++);
			to += NumRows;
		}
		rowPtr ++;
	}	
}

// Read out a submatrix into "result"
// The size of "result" determines the size of the submatrix read out.
void MatrixRmn::GetSubmatrix( long i, long j, MatrixRmn* result ) const
{
	assert ( 0<=i && i+result->NumRows <= NumRows );
	assert ( 0<=j && j+result->NumCols <= NumCols );
	const double* fromColPtr = x + j*NumRows + i;
	double* to = result->x;
	for ( long i=result->NumCols; i>0; i-- ) {
		const double* from = fromColPtr;
		for ( long j=result->NumRows; j>0; j-- ) {
			*(to++) = *(from++);
		}
		fromColPtr += NumRows;
	}
}


bool MatrixRmn::CheckDiagPositive() const 
{
	const double* dPtr = GetPtr();
	long nRows = GetNumRows();
	long nCols = GetNumColumns();
	long dJump = nRows+1;
	long i;
	for ( i=Min(nCols,nRows); i>0; i-- ) {
		if ( (*(dPtr)) <= 0.0 ) {
			return false;
		}
		dPtr += dJump;
	}
	return true;
}

double MatrixRmn::MinDiagonalElement() const
{
	const double* dPtr = GetPtr();
	long nRows = GetNumRows();
	long nCols = GetNumColumns();
	long dJump = nRows+1;
	double ret = *dPtr;
	for ( long i=Min(nCols,nRows); i>1; i-- ) {
		dPtr += dJump;
		UpdateMin( *dPtr, ret );
	}
	return ret;
}

double MatrixRmn::MinDiagonalElement( long *indexOfMin ) const
{
	const double* dPtr = GetPtr();
	long nRows = GetNumRows();
	long nCols = GetNumColumns();
	long dJump = nRows+1;
	double ret = *dPtr;
	*indexOfMin = 0;
	long diagLen = Min(nCols, nRows);
	for ( long i=1; i<diagLen; i++ ) {
		dPtr += dJump;
		double thisEntry = *dPtr;
		if ( thisEntry<ret ) {
			ret = thisEntry;
			*indexOfMin = i;
		}
	}
	return ret;
}


// Calculate the Frobenius Norm (square root of sum of squares of entries of the matrix)
double MatrixRmn::FrobeniusNorm() const
{
	return sqrt( FrobeniusNormSq() );
}

// Calculate the maximum sum of absolute values in any single row 
double MatrixRmn::MaxRowSum() const
{
	const double* rowPtr = GetPtr();
	long nRows = GetNumRows();
	long nCols = GetNumColumns();
	double maxRowSum;
	for ( long i=nRows; i>0; i-- ) {
		const double* xPtr = (rowPtr++);
		double sum = 0.0;
		for ( long j=nCols; j>0; j-- ) {
			sum += fabs( *xPtr );
			xPtr += nRows;
		}
		UpdateMax( sum, maxRowSum );
	}
	return maxRowSum;
}

// Returns lower bound on minimum eigenvalue (using rows)
double MatrixRmn::GershgorinEigenvalueLowerBound() const
{
	assert( GetNumRows() == GetNumColumns() );
	long nRows = GetNumRows();
	const double *rowPtr = GetPtr();
	const double *diagPtr = rowPtr;
	double ret = DBL_MAX;
	for ( long i=0; i<nRows; i++ ) {
		const double *rElts = rowPtr;
		double thisLB = *diagPtr;
		for ( long j=0; j<i; j++ ) {
			thisLB -= fabs( *rElts );
			rElts += nRows;
		}
		for ( long j=i+1; j<nRows; j++ ) {
			rElts += nRows;
			thisLB -= fabs( *rElts );
		}
		UpdateMin( thisLB, ret );
		diagPtr += (nRows+1);
		rowPtr ++;
	}
	return ret;
}


// Multiply this matrix by column vector v.  
// Result is column vector "result".  "result" must already have the correct size
void MatrixRmn::Multiply( const VectorRn& v, VectorRn& result ) const
{
	assert ( v.GetLength()==NumCols && result.GetLength()==NumRows );
	double* out = result.GetPtr();				// Points to entry in result vector
	const double* rowPtr = x;					// Points to beginning of next row in matrix
	for ( long j = NumRows; j>0; j-- ) {
		const double* in = v.GetPtr();
		const double* m = rowPtr++;
		*out = 0.0f;
		for ( long i = NumCols; i>0; i-- ) {
			*out += (*(in++)) * (*m);
			m += NumRows;
		}
		out++;
	}
}

// Multiply transpose of this matrix by column vector v.  
//    Result is column vector "result"
// Equivalent to mult by row vector on left
void MatrixRmn::MultiplyTranspose( const VectorRn& v, VectorRn& result ) const	
{
	assert ( v.GetLength()==NumRows && result.GetLength()==NumCols );
	double* out = result.GetPtr();				// Points to entry in result vector
	const double* colPtr = x;					// Points to beginning of next column in matrix
	for ( long i=NumCols; i>0; i-- ) {
		const double* in=v.GetPtr();
		*out = 0.0f;
		for ( long j = NumRows; j>0; j-- ) {
			*out += (*(in++)) * (*(colPtr++));
		}
		out++;
	}
}

// Form the dot product of a vector v with the i-th column of the array
double MatrixRmn::DotProductColumn( const VectorRn& v, long colNum ) const
{
	assert ( v.GetLength()==NumRows );
	double* ptrC = x+colNum*NumRows;
	double* ptrV = v.x;	
	double ret = 0.0;
	for ( long i = NumRows; i>0; i-- ) {
		ret += (*(ptrC++))*(*(ptrV++));
	}
	return ret;
}

// Multiply the strict upper triangular part of this matrix times column vector v.
//		The diagonal or subdiagonal entries are ignored.
void MatrixRmn::MultiplyStrictUpperTriangular( const VectorRn& v, VectorRn& result ) const
{
	assert ( v.GetLength()==NumCols && result.GetLength()==NumRows );
	double* out = result.GetPtr();				// Points to entry in result vector
	long diagStep = NumRows+1;
	const double* rowPtr = x+NumRows;			// Points to first superdiagonal entry in row
	long upperLen = NumCols-1;					// Number of superdiagonal entries in row
	const double* vPtr = v.GetPtr()+1;			// Pointer to v entry for 1st nonzero multiply
	for ( long j = NumRows; j>0; j-- ) {
		const double* in = vPtr;
		const double* m = rowPtr;
		*out = 0.0f;
		for ( long i = upperLen; i>0; i-- ) {
			*out += (*(in++)) * (*m);
			m += NumRows;
		}
		rowPtr += diagStep;
		upperLen--;
		vPtr++;
		out++;
	}
}

// Multiply the lower triangular part of this matrix times column vector v.
//		The diagonal is used, but the superdiagonal entries are ignored.
void MatrixRmn::MultiplyLowerTriangular( const VectorRn& v, VectorRn& result ) const
{
	assert ( v.GetLength()==NumCols && result.GetLength()==NumRows );
	double* out = result.GetPtr();				// Points to entry in result vector
	const double* rowPtr = x;					// Points to beginning of next row in matrix
	long ltLen = 1;								// Number of lower triangular entries in row
	for ( long j = NumRows; j>0; j-- ) {
		const double* in = v.GetPtr();
		const double* m = rowPtr++;
		*out = 0.0f;
		for ( long i = ltLen; i>0; i-- ) {
			*out += (*(in++)) * (*m);
			m += NumRows;
		}
		ltLen++;
		out++;
	}
}


// Add a constant to each entry on the diagonal
MatrixRmn& MatrixRmn::AddToDiagonal( double d )				// Adds d to each diagonal entry
{
	long diagLen = Min( NumRows, NumCols );
	double* dPtr = x;
	for ( ; diagLen>0; diagLen-- ) {
		*dPtr += d;
		dPtr += NumRows+1;
	}
	return *this;
}

// Multiply two MatrixRmn's
MatrixRmn& MatrixRmn::Multiply( const MatrixRmn& A, const MatrixRmn& B, MatrixRmn& dst )
{
	assert( A.NumCols == B.NumRows && A.NumRows == dst.NumRows && B.NumCols == dst.NumCols );
	long length = A.NumCols;

	double *bPtr = B.x;						// Points to beginning of column in B
	double *dPtr = dst.x;
	for ( long i = dst.NumCols; i>0; i-- ) {
		double *aPtr = A.x;					// Points to beginning of row in A
		for ( long j = dst.NumRows; j>0; j-- ) {
			*dPtr = DotArray( length, aPtr, A.NumRows, bPtr, 1 );
			dPtr++;
			aPtr++;
		}
		bPtr += B.NumRows;
	}

	return dst;
}

// Multiply two MatrixRmn's,  Transpose the first matrix before multiplying.  (1st matrix is unchanged)
MatrixRmn& MatrixRmn::TransposeMultiply( const MatrixRmn& A, const MatrixRmn& B, MatrixRmn& dst )
{
	assert( A.NumRows == B.NumRows && A.NumCols == dst.NumRows && B.NumCols == dst.NumCols );
	long length = A.NumRows;

	double *bPtr = B.x;										// bPtr Points to beginning of column in B
	double *dPtr = dst.x;
	for ( long i = dst.NumCols; i>0; i-- ) {				// Loop over all columns of dst
		double *aPtr = A.x;									// aPtr Points to beginning of column in A
		for ( long j = dst.NumRows; j>0; j-- ) {			// Loop over all rows of dst
			*dPtr = DotArray( length, aPtr, 1, bPtr, 1 );
			dPtr ++;
			aPtr += A.NumRows;
		}
		bPtr += B.NumRows;
	}

	return dst;
}

// Multiply two MatrixRmn's.  Transpose the second matrix before multiplying. (2nd matrix is unchanged)
MatrixRmn& MatrixRmn::MultiplyTranspose( const MatrixRmn& A, const MatrixRmn& B, MatrixRmn& dst )
{
	assert( A.NumCols == B.NumCols && A.NumRows == dst.NumRows && B.NumRows == dst.NumCols );
	long length = A.NumCols;

	double *bPtr = B.x;						// Points to beginning of row in B
	double *dPtr = dst.x;
	for ( long i = dst.NumCols; i>0; i-- ) {
		double *aPtr = A.x;					// Points to beginning of row in A
		for ( long j = dst.NumRows; j>0; j-- ) {
			*dPtr = DotArray( length, aPtr, A.NumRows, bPtr, B.NumRows );
			dPtr++;
			aPtr++;
		}
		bPtr++;
	}

	return dst;
}

// Multiply two MatrixRmn's.  The second matrix, B, is tridiagonal.  Rest of B is ignored.
MatrixRmn& MatrixRmn::MultiplyTridiagonal( const MatrixRmn& A, const MatrixRmn& B, MatrixRmn& dst )
{
	assert( A.NumCols == B.NumRows && A.NumRows == dst.NumRows && B.NumCols == dst.NumCols );

	// First colume of product
	double *dPtr = dst.x;
	double b1 = *(B.x);
	double b2 = *(B.x+1);
	const double* aPtr = A.x;
	for ( long j=0; j<dst.NumRows; j++ ) {
		*dPtr = b1*(*aPtr) + b2*(*(aPtr+dst.NumRows));
		dPtr++;
		aPtr++;
	}

	// Middle columns
	// dPtr now points to beginning of second column of dst
	const double *bPtr = B.x+B.NumRows;		// Beginning of second column of B
	long offDiagStep = B.NumRows+1;
	aPtr = A.x;
	for ( long i = dst.NumCols; i>2; i-- ) {
		double b0 = *bPtr;
		double b1 = *(bPtr+1);
		double b2 = *(bPtr+2);
		for ( long j = dst.NumRows; j>0; j-- ) {
			*dPtr = b0*(*aPtr) + b1*(*(aPtr+A.NumCols)) + b2*(*(aPtr+(A.NumCols+A.NumCols)));
			dPtr++;
			aPtr++;
		}
		bPtr += offDiagStep;					// Next superdiagonal entry in B
	}

	// Last column
	double b0 = *bPtr;
	b1 = *(bPtr+1);
	for ( long j = dst.NumRows; j>0; j-- ) {
		*dPtr = b0*(*aPtr) + b1*(*(aPtr+A.NumCols));
		dPtr++;
		aPtr++;
	}

	return dst;
}

// Outer products
MatrixRmn& MatrixRmn::SetOuterProduct( const VectorRn& v )
{
	long n = v.GetLength();
	SetSize(n,n);
	double* toPtr = x;
	const double* fPtr1 = v.x;
	for ( long i=n; i>0; i-- ) {
		const double* fPtr2 = v.x;
		for ( long j=n; j>0; j-- ) {
			*(toPtr++) = (*fPtr1)*(*(fPtr2++));
		}
		fPtr1++;
	}
	return *this;
}

//   Matrix should already be the correct size
MatrixRmn& MatrixRmn::AddOuterProduct( const VectorRn& v )
{
	long n = v.GetLength();
	assert ( n==NumRows && n==NumCols );
	double* toPtr = x;
	const double* fPtr1 = v.x;
	for ( long i=n; i>0; i-- ) {
		const double* fPtr2 = v.x;
		for ( long j=n; j>0; j-- ) {
			*(toPtr++) += (*fPtr1)*(*(fPtr2++));
		}
		fPtr1++;
	}
	return *this;
}



// Solves the equation   (*this)*xVec = b;  
// Uses row operations.  Assumes *this is square and invertible.   
// No error checking for divide by zero or instability (except with asserts)
void MatrixRmn::Solve( const VectorRn& b, VectorRn* xVec ) const
{
	assert ( NumRows==NumCols && NumCols==xVec->GetLength() && NumRows==b.GetLength() );

	// Copy this matrix and b into an Augmented Matrix
	MatrixRmn& AugMat = GetWorkMatrix( NumRows, NumCols+1 );
	AugMat.LoadAsSubmatrix( *this );
	AugMat.SetColumn( NumRows, b );

	// Put into row echelon form with row operations
	AugMat.ConvertToRefNoFree();

	// Solve for x vector values using back substitution
	double* xLast = xVec->x+NumRows-1;				// Last entry in xVec
	double* endRow = AugMat.x+NumRows*NumCols-1;	// Last entry in the current row of the coefficient part of Augmented Matrix
	double* bPtr = endRow+NumRows;				// Last entry in augmented matrix (end of last column, in augmented part)
	for ( long i = NumRows; i>0; i-- ) {
		double accum = *(bPtr--);
		// Next loop computes back substitution terms
		double* rowPtr = endRow;					// Points to entries of the current row for back substitution.
		double* xPtr = xLast;						// Points to entries in the x vector (also for back substitution)
		for ( long j=NumRows-i; j>0; j-- ) {			
			accum -= (*rowPtr)*(*(xPtr--));
			rowPtr -= NumCols;						// Previous entry in the row
		}
		assert( *rowPtr != 0.0 );					// Are not supposed to be any free variables in this matrix
		*xPtr = accum/(*rowPtr);
		endRow--;
	}

}

// SolveLowerDiagonal: Assumes (*this) lower diagonal.  Solves (*this) x  = b;
// All diagonal entries must be nonzero (i.e., matrix must be invertible)
void MatrixRmn::SolveLowerDiagonal( const VectorRn& b, VectorRn* ret ) const
{
	long n = NumRows;
	assert ( NumRows==NumCols );	// Matrix must be square
	assert ( b.GetLength() == n );
	ret->SetLength( n );

	const double* rowPtr = x;
	const double* bPtr = b.x;
	for ( long i=0; i<n; i++ ) {
		double accum = *bPtr;
		double* retPtr = ret->x;
		const double* tPtr = rowPtr;
		for ( long j=0; j<i; j++ ) {
			accum -= (*retPtr)*(*tPtr);
			retPtr++;
			tPtr += n;
		}
		assert ( *tPtr != 0.0 );
		*retPtr = accum / *tPtr;
		
		bPtr++;
		rowPtr++;
	}

}


// SolveLowerDiagonalTranspose: Assumes (*this) lower diagonal.  Solves (*this)^T x  = b;
// All diagonal entries must be nonzero (i.e., matrix must be invertible)
void MatrixRmn::SolveLowerDiagonalTranspose( const VectorRn& b, VectorRn* ret ) const
{
	long n = NumRows;
	assert ( NumRows==NumCols );	// Matrix must be square
	assert ( b.GetLength() == n );
	ret->SetLength( n );
	long diagStep = n+1;

	const double* diagPtr = x + NumRows*NumRows - 1;  // Last entry in diagonal
	const double* colEndPtr = diagPtr;				  // Last entry in last column
	const double* bPtr = b.x + (n - 1);				  // Last entry in b
	double* retEndPtr = ret->x + (n - 1);			  // Last entry in x
	for ( int i = n-1; i>=0; i-- ) {
		double accum = *bPtr;
		const double* colPtr = colEndPtr;
		double* retPtr = retEndPtr;
		for ( int j = n-1; j>i; j-- ) {
			accum -= (*(retPtr--))*(*(colPtr--));
		}
		assert ( (*diagPtr) != 0.0 );				// Diagonal entry not allowed to be zero.
		*retPtr = accum /(*diagPtr);
		bPtr--;										// (i-1)th entry in b
		colEndPtr -= n;								// Last entry of next column
		diagPtr -= diagStep;						// Next entry in diagonal (backwards)
	}
}


// Invert a matrix with row operations (on an augmented matrix)
void MatrixRmn::InverseWithRowOps( MatrixRmn* inverse ) const
{
	assert ( NumRows==NumCols );

	// Copy matrix into an augmented matrix.
	// The other half of the augmented matrix holds the identity
	MatrixRmn& AugMat = GetWorkMatrix( NumRows, 2*NumRows );
	AugMat.LoadAsSubmatrix( *this );
	double* augPart = AugMat.x + NumRows*NumRows;
	double* to = augPart;
	for ( long i = 0; i<NumRows*NumRows; i++ ) {
		*(to++) = 0.0;
	}
	to = augPart;
	long diagStep = NumRows+1;
	for ( long i = 0; i<NumRows; i++ ) {
		*to = 1.0;
		to += diagStep;
	}

	// Do the row operations, with row pivoting
	AugMat.RowConvertToIdentity();

	// Return the inverse from the augmented part of the matrix.
	inverse->SetSize( NumRows, NumCols );
	AugMat.GetSubmatrix( 0, NumRows, inverse );
}


// ConvertToRefNoFree
// Converts the matrix (in place) to row echelon form
// For us, row echelon form allows any non-zero values, not just 1's, in the 
//		position for a lead variable.
// The "NoFree" version operates on the assumption that no free variable will be found.
// Algorithm uses row operations and row pivoting (only).  
// Augmented matrix is correctly accomodated.  Only the first square part participates 
//		in the main work of row operations.
void MatrixRmn::ConvertToRefNoFree()
{
	// Loop over all columns (variables)
	// Find row with most non-zero entry.
	// Swap to the highest active row
	// Subtract appropriately from all the lower rows (row op of type 3)
	long numIters = Min(NumRows,NumCols);		// It would be silly if NumRows were >= NumCols, but OK.
	double* rowPtr1 = x;
	const long diagStep = NumRows+1;
	long lenRowLeft = NumCols;
	for ( ; numIters>1; numIters-- ) {
		// Find row with most non-zero entry.
		double* rowPtr2 = rowPtr1;
		double maxAbs = fabs(*rowPtr1);
		double *rowPivot = rowPtr1;
		for ( long i=numIters-1; i>0; i-- ) {
			const double& newMax = *(++rowPivot);
			if ( newMax > maxAbs ) {
				maxAbs = *rowPivot;
				rowPtr2 = rowPivot;
			}
			else if ( -newMax > maxAbs ) {
				maxAbs = -newMax;
				rowPtr2 = rowPivot;
			}
		}
		// Pivot step: Swap the row with highest entry to the current row
		if ( rowPtr1 != rowPtr2 ) {
			double *to = rowPtr1;
			for ( long i=lenRowLeft; i>0; i-- ) {
				double temp = *to;
				*to = *rowPtr2;
				*rowPtr2 = temp;
				to += NumRows;
				rowPtr2 += NumRows;
			}
		}
		// Subtract this row appropriately from all the lower rows (row operation of type 3)
		rowPtr2 = rowPtr1;
		for ( long i=numIters-1; i>0; i-- ) {
			rowPtr2++;
			double* to = rowPtr2;
			double* from = rowPtr1;
			assert( *from != 0.0 );
			double alpha = (*to)/(*from);
			*to = 0.0;
			for ( long j=lenRowLeft-1; j>0; j-- ) {
				to += NumRows;
				from += NumRows;
				*to -= (*from)*alpha;
			}
		}
		// Update for next iteration of loop
		rowPtr1 += diagStep;
		lenRowLeft--;
	}
}

// ConvertToRrefNoFree
// Converts the matrix (in place) to reduced row echelon form
// For us, reduced row echelon form requires 1's in the 
//		position for a lead variable.
// The "NoFree" version operates on the assumption that no free variable will be found.
// Algorithm uses row operations and row pivoting (only).  
// Augmented matrix is correctly accomodated.  Only the first square part participates 
//		in the main work of row operations.
void MatrixRmn::ConvertToRrefNoFree()
{
	ConvertToRefNoFree();

	// Loop over each row and multiply by the appropriate
	//	factor to make the lead entry equal to one.

	long numIters = Min(NumRows,NumCols);		// It would be silly if NumRows were >= NumCols, but OK.
	double* diagPtr = x;
	long i;
	for (i=0; i<numIters; i++ ) {
		assert ( (*diagPtr) != 0 );			// Should never happen if ConvertToRefNoFree succeeded
		double factor = 1.0/(*diagPtr);
		*diagPtr = 1.0;						// Set diagonal element equal to one
		double* rowPtr = diagPtr+NumRows;
		long k;
		for ( k=i+1; k<NumCols; k++ ) {
			*rowPtr *= factor;
			rowPtr += NumRows;				// Point to next entry in row
		}
		diagPtr += NumRows+1;				// Next diagonal element ptr
	}
}

// Row convert to identity is useful for computing inverses.
// Uses row operations to convert first square part to the identity.
// Algorithm uses row operations and row pivoting (only).  
// Augmented matrix is correctly accomodated.  Only the first square part participates 
//		in the main work of row operations.
// No checking for singularity except thru asserts!
void MatrixRmn::RowConvertToIdentity()
{
	ConvertToRefNoFree();

	// Loop over each row and multiply by the appropriate
	//	factor to make the lead entry equal to one.

	long numIters = Min(NumRows,NumCols);		// It would be silly if NumRows were >= NumCols, but OK.
	double* diagPtr = x + NumRows*(numIters-1) + numIters - 1;
	double* restThisRow = diagPtr+NumRows; // Rest of row after first square block
	long i;
	for (i=0; i<numIters; i++ ) {
		assert ( (*diagPtr) != 0 );			// Should never happen if ConvertToRefNoFree succeeded
		double diagInv = 1.0/(*diagPtr);
		// Subtract multiple of current row from all earlier rows
		for ( long j = 1; j < numIters-i; j++ ) {
			double factor = (*(diagPtr-j))*diagInv;
			*(diagPtr-j) = 0.0;
			double* thisRowPtr = restThisRow;
			double* prevRowPtr = restThisRow - j;
			for ( long k=numIters; k<NumCols; k++ ) {
				*prevRowPtr -= (*thisRowPtr)*factor;
				prevRowPtr += NumRows;
				thisRowPtr += NumRows;
			}

		}
		// Multiply rest of this row by 1/(diagonal value)
		double* thisRowPtr = restThisRow;
		for ( long k=numIters; k<NumCols; k++ ) {
			*thisRowPtr *= diagInv;
			thisRowPtr += NumRows;			// Point to next entry in row
		}
		*diagPtr = 1.0;						// Set diagonal element equal to one
		diagPtr -= (NumRows+1);				// Next diagonal element ptr
		restThisRow --;
	}
}

// ******************************************************************
// Exponentiation
// Compute   exp(A*t)
void MatrixRmn::CalcExponential( MatrixRmn& retExpOfAt, double t ) const
{
	long n = GetNumRows();
	assert ( GetNumRows()==GetNumColumns() && GetNumRows()>0 );		// Makes no sense if matrix is not square
	retExpOfAt.SetSize(n,n);

	// Scale down by power of two to have norm at most one
	// Place result in WorkMatrix
	double tnorm = fabs(t)*MaxRowSum();
	int reps = 0;
	WorkMatrix = *this;
	if ( tnorm>1.0 ) {
		double log2tnorm = log(tnorm)*LnTwoInv;
		reps = 1+(int)log2tnorm;
		double scale = 1.0/pow(2.0, (double)reps);
		tnorm *= scale;			// Use for stopping criteria
		WorkMatrix *= t*scale;
	}
	else { 
		WorkMatrix *= t;
	}

	// Set retExpOfAt equal to exp(WorkMatrix)
	WorkMatrix2 = WorkMatrix;	// WorkMatrix2 holds (WorkMatrix)^i
	WorkMatrix3.SetSize(n,n);	// Preset before call to Multiply below
	retExpOfAt.SetSize(n,n);
	retExpOfAt.SetIdentity();
	double factorialInv=1.0;
	double stopScale = tnorm;
	int i;
	for ( i=1; ; i++ ) {
		factorialInv /= (double)i;
		retExpOfAt.AddScaled( WorkMatrix2, factorialInv );
		if ( 1.0+stopScale*factorialInv == 1.0 ) {
			break;			// A rather conservative test (but simple and fast)
		}
		stopScale *= tnorm;
		Multiply( WorkMatrix2, WorkMatrix, WorkMatrix3 );
		WorkMatrix2 = WorkMatrix3;
	}

	// Use iterated squaring to raise to the right power.
	for ( i=0; i<reps; i++ ) {
		WorkMatrix = retExpOfAt;
		Multiply( WorkMatrix, WorkMatrix, retExpOfAt );
	}
	return;
}

// **************************************************************************
// Cholesky decomposition of positive semidefinite, symmetric matrix
//   Finds G such that  this ==  G * G^T, G lower diagonal.
//   Returns "true" if successful.
bool MatrixRmn::CalcCholesky( MatrixRmn& retCholesky ) const
{
	long successColumn = CalcCholeskyInner( retCholesky );
	return ( NumRows==successColumn );
}

// This version returns a vector v = nonPosExample, such that v^T (this) v = nonPosAmount < 0
// It fails for non-positive matrices and returns false.
// In this case, it also returns a "nonPosExample" vector and a "nonPosAmount".
// The vector u = nonPosExample has the property   u^T (*this) u  =  nonPosAmount/u.NormSq().
//     nonPosExample is usually not a unit vector.
// To make nonnegative definite, must add at least  nonPosAmount/nonPosExample.NormSq() to the diagonal entries.
bool MatrixRmn::CalcCholesky( MatrixRmn& retCholesky, VectorRn* nonPosExample, double* nonPosAmount ) const
{
	long successColumn = CalcCholeskyInner( retCholesky );
	if ( successColumn==NumRows ) {
		return true;		// Succeeded
	}
	nonPosExample->SetLength( NumRows );
	double* nPEbase = nonPosExample->GetPtr(successColumn);
	double* nPEptr = nPEbase;
	*nPEptr = 1.0;
	for ( long j=successColumn+1; j<NumRows; j++ ){
		*(++nPEptr) = 0.0;
	}
	const double* rCbase = retCholesky.GetPtr(successColumn, successColumn);
	*nonPosAmount = *rCbase;
	for ( long i=successColumn-1; i>=0; i-- ) {
		nPEptr = nPEbase;
		rCbase -= NumRows;
		const double* rCptr = rCbase;		// Position (successColumn,i) in retCholesky
		double sum = *rCptr;
		for ( long j = successColumn-1; j>i; j-- ) {
			sum += (*(--rCptr))*(*(--nPEptr));
		}
		(*(--nPEptr)) = -sum/(*(--rCptr));
	}
	return false;
}


long MatrixRmn::CalcCholeskyInner( MatrixRmn& retCholesky ) const
{
	long n = NumRows;
	assert ( n == NumCols && n>0 );		// Should be a square matrix
	retCholesky.SetSize( n, n );
	retCholesky.SetZero();					// Zero so strict upper diagonal part will be zero.
	long diagStep = n+1;

	const double* fromBase = x;
	double* toBase = retCholesky.x;
	for ( long i=0; i<n; i++ ) {
		// Set i-th column of G
		const double* fromPtr = fromBase;
        double* toPtr = toBase;
		for ( long j=i; j<n; j++ ) {
			*(toPtr++) = *(fromPtr++);	// Copy i-th column into G
		}
		// Subtract off amounts from first i-1 columns of G
		fromPtr = retCholesky.x + i;		// Beginning of row i in G
		const double* multPtr = fromPtr;	// Ditto
		for ( long k=0; k<i; k++ ) {
			toPtr = toBase;
			double factor = *multPtr;
			for ( long j=i; j<n; j++ ) {
				*(toPtr++) -= (*(fromPtr++))*factor;
			}
			fromPtr += i;	// Skip first i entries in next column
			multPtr += n;	// Point to next row entry in G
		}
		// Check sign, if positive, proceed, dividing by sqrt
		double diagEntry = *toBase;
		if ( diagEntry<=0 ) {
			return i;
		}
		diagEntry = sqrt(diagEntry);
		*toBase = diagEntry;
		double factor = 1.0/diagEntry;
        toPtr = toBase+1;
		for ( long j=i+1; j<n; j++ ) {
			(*(toPtr++)) *= factor;
		}
		fromBase += diagStep;
		toBase += diagStep;
	}

	return n;		// Return n to indicate success
}


// *****************************************************************************
// Givens transformations
// Calculate the c=cosine and s=sine values for a Givens transformation.
// The matrix M = ( (c, -s), (s, c) ) in row order transforms the
//   column vector (a, b)^T to have y-coordinate zero.
void MatrixRmn::CalcGivensValues( double a, double b, double *c, double *s )
{
	double denomInv = sqrt(a*a + b*b);
	if ( denomInv==0.0 ) {
		*c = 1.0;
		*s = 0.0;
	}
	else {
		denomInv = 1.0/denomInv;
		*c = a*denomInv;
		*s = -b*denomInv;
	}
}

// Applies Givens transform to columns i and i+1.
// Equivalent to postmultiplying by the matrix
//      ( c  -s )
//		( s   c )
// with non-zero entries in rows i and i+1 and columns i and i+1
void MatrixRmn::PostApplyGivens( double c, double s, long idx )
{
	assert ( 0<=idx && idx<NumCols );
	double *colA = x + idx*NumRows;
	double *colB = colA + NumRows;
	for ( long i = NumRows; i>0; i-- ) {
		double temp = *colA;
		*colA = (*colA)*c + (*colB)*s;
		*colB = (*colB)*c - temp*s;
		colA++;
		colB++;
	}
}	

// Applies Givens transform to columns idx1 and idx2.
// Equivalent to postmultiplying by the matrix
//      ( c  -s )
//		( s   c )
// with non-zero entries in rows idx1 and idx2 and columns idx1 and idx2
void MatrixRmn::PostApplyGivens( double c, double s, long idx1, long idx2 )
{
	assert ( idx1!=idx2 && 0<=idx1 && idx1<NumCols && 0<=idx2 && idx2<NumCols );
	double *colA = x + idx1*NumRows;
	double *colB = x + idx2*NumRows;
	for ( long i = NumRows; i>0; i-- ) {
		double temp = *colA;
		*colA = (*colA)*c + (*colB)*s;
		*colB = (*colB)*c - temp*s;
		colA++;
		colB++;
	}
}	


// ********************************************************************************************
// Singular value decomposition.
// Return othogonal matrices U and V and diagonal matrix with diagonal w such that
//     (this) = U * Diag(w) * V^T     (V^T is V-transpose.)
// Diagonal entries have all non-zero entries before all zero entries, but are not
//		necessarily sorted.  (Someday, I will write ComputedSortedSVD that handles
//		sorting the eigenvalues by magnitude.)
// ********************************************************************************************
void MatrixRmn::ComputeSVD( MatrixRmn& U, VectorRn& w, MatrixRmn& V ) const
{
	assert ( U.NumRows==NumRows && V.NumCols==NumCols 
			 && U.NumRows==U.NumCols && V.NumRows==V.NumCols
			 && w.GetLength()==Min(NumRows,NumCols) );

	VectorRn& superDiag = VectorRn::GetWorkVector( w.GetLength()-1 );		// Some extra work space.  Will get passed around.

	// Choose larger of U, V to hold intermediate results
	// If U is larger than V, use U to store intermediate results
	// Otherwise use V.  In the latter case, we form the SVD of A transpose,
	//		(which is essentially identical to the SVD of A).
	MatrixRmn* leftMatrix;
	MatrixRmn* rightMatrix;
	if ( NumRows >= NumCols ) {
		U.LoadAsSubmatrix( *this );				// Copy A into U
		leftMatrix = &U;
		rightMatrix = &V;
	}
	else {
		V.LoadAsSubmatrixTranspose( *this );		// Copy A-transpose into V
		leftMatrix = &V;
		rightMatrix = &U;
	}

	// Do the actual work to calculate the SVD 
	// Now matrix has at least as many rows as columns
	CalcBidiagonal( *leftMatrix, *rightMatrix, w, superDiag );
	ConvertBidiagToDiagonal( *leftMatrix, *rightMatrix, w, superDiag );

}

// ************************************************ CalcBidiagonal **************************
// Helper routine for SVD computation
// U is a matrix to be bidiagonalized.
// On return, U and V are orthonormal and w holds the new diagonal
//	  elements and superDiag holds the super diagonal elements.
	
void MatrixRmn::CalcBidiagonal( MatrixRmn& U, MatrixRmn& V, VectorRn& w, VectorRn& superDiag )
{
	assert ( U.NumRows>=V.NumRows );

	// The diagonal and superdiagonal entries of the bidiagonalized
	//	  version of the U matrix
	//	  are stored in the vectors w and superDiag (temporarily).

	// Apply Householder transformations to U.
	// Householder transformations come in pairs.
	//   First, on the left, we map a portion of a column to zeros
	//   Second, on the right, we map a portion of a row to zeros
	const long rowStep = U.NumCols;
	const long diagStep = U.NumCols+1;
	double *diagPtr = U.x;
	double* wPtr = w.x;
	double* superDiagPtr = superDiag.x;
	long colLengthLeft = U.NumRows;
	long rowLengthLeft = V.NumCols;
	while (true) {
		// Apply a Householder xform on left to zero part of a column
		SvdHouseholder( diagPtr, colLengthLeft, rowLengthLeft, 1, rowStep, wPtr ); 

		if ( rowLengthLeft==2 ) {
			*superDiagPtr = *(diagPtr+rowStep);
			break;
		}
		// Apply a Householder xform on the right to zero part of a row
		SvdHouseholder( diagPtr+rowStep, rowLengthLeft-1, colLengthLeft, rowStep, 1, superDiagPtr ); 

		rowLengthLeft--;
		colLengthLeft--;
		diagPtr += diagStep;
		wPtr++;
		superDiagPtr++;
	}

	long extra = 0;
	diagPtr += diagStep;
	wPtr++;
	if ( colLengthLeft > 2 ) {
		extra = 1;
		// Do one last Householder transformation when the matrix is not square
		colLengthLeft--;
		SvdHouseholder( diagPtr, colLengthLeft, 1, 1, 0, wPtr );
	}
	else {
		*wPtr = *diagPtr;
	}

	// Form U and V from the Householder transformations
	V.ExpandHouseholders( V.NumCols-2, 1, U.x+U.NumRows, U.NumRows, 1 );
	U.ExpandHouseholders( V.NumCols-1+extra, 0, U.x, 1, U.NumRows );

	// Done with bidiagonalization
	return;
}

// Helper routine for CalcBidiagonal
// Performs a series of Householder transformations on a matrix
// Stores results compactly into the matrix:   The Householder vector u (normalized)
//   is stored into the first row/column being transformed.
// The leading term of that row (= plus/minus its magnitude is returned
//	 separately into "retFirstEntry"
void MatrixRmn::SvdHouseholder( double* basePt, 
							    long colLength, long numCols, long colStride, long rowStride, 
								double* retFirstEntry )
{

	// Calc norm of vector u
	double* cPtr = basePt;
	double norm = 0.0;
	long i;
	for ( i=colLength; i>0 ; i-- ) {
		norm += Square( *cPtr );
		cPtr += colStride;
	}
	norm = sqrt(norm);					// Norm of vector to reflect to axis  e_1

	// Handle sign issues		
	double imageVal;					// Choose sign to maximize distance
	if ( (*basePt) < 0.0 ) {
		imageVal = norm;
		norm = 2.0*norm*(norm-(*basePt));	
	}
	else {
		imageVal = -norm;
		norm = 2.0*norm*(norm+(*basePt));
	}
	norm = sqrt(norm);					// Norm is norm of reflection vector

	if ( norm==0.0 ) {			// If the vector being transformed is equal to zero
		// Force to zero in case of roundoff errors
		cPtr = basePt;
		for ( i=colLength; i>0; i-- ) {
			*cPtr = 0.0;
			cPtr += colStride;
		}
		*retFirstEntry = 0.0;
		return;
	}

	*retFirstEntry = imageVal;

	// Set up the normalized Householder vector
	*basePt -= imageVal;					// First component changes. Rest stay the same.
	// Normalize the vector
	norm = 1.0/norm;					// Now it is the inverse norm
	cPtr = basePt;
	for ( i=colLength; i>0 ; i-- ) {
		*cPtr *= norm;
		cPtr += colStride;
	}

	// Transform the rest of the U matrix with the Householder transformation
	double *rPtr = basePt;
	for ( long j=numCols-1; j>0; j-- ) {
		rPtr += rowStride;
		// Calc dot product with Householder transformation vector
		double dotP = DotArray( colLength, basePt, colStride, rPtr, colStride );
		// Transform with I - 2*dotP*(Householder vector)
		AddArrayScale( colLength, basePt, colStride, rPtr, colStride, -2.0*dotP );
	}
}

// ********************************* ExpandHouseholders ********************************************
// The matrix will be square.
//   numXforms = number of Householder transformations to concatenate
//		Each Householder transformation is represented by a unit vector
//		Each successive Householder transformation starts one position later
//			and has one more implied leading zero
//	 basePt = beginning of the first Householder transform
//	 colStride, rowStride: Householder xforms are stored in "columns"
//   numZerosSkipped is the number of implicit zeros on the front each
//			Householder transformation vector (only values supported are 0 and 1).
void MatrixRmn::ExpandHouseholders( long numXforms, long numZerosSkipped, const double* basePt, long colStride, long rowStride )
{
	// Number of applications of the last Householder transform
	//     (That are not trivial!)
	long numToTransform = NumCols-numXforms+1-numZerosSkipped;			
	assert( numToTransform>0 );

	if ( numXforms==0 ) {
		SetIdentity();
		return;
	}

	// Handle the first one separately as a special case,
	// "this" matrix will be treated to simulate being preloaded with the identity
	long hDiagStride = rowStride+colStride;
	const double* hBase = basePt + hDiagStride*(numXforms-1);	// Pointer to the last Householder vector
	const double* hDiagPtr = hBase + colStride*(numToTransform-1);		// Pointer to last entry in that vector
	long i;
	double* diagPtr = x+NumCols*NumRows-1;					// Last entry in matrix (points to diagonal entry)
	double* colPtr = diagPtr-(numToTransform-1);			// Pointer to column in matrix
	for ( i=numToTransform; i>0; i-- ) {
		CopyArrayScale( numToTransform, hBase, colStride, colPtr, 1, -2.0*(*hDiagPtr) );
		*diagPtr += 1.0;						// Add back in 1 to the diagonal entry (since xforming the identity)
		diagPtr -= (NumRows+1);					// Next diagonal entry in this matrix
		colPtr -= NumRows;						// Next column in this matrix
		hDiagPtr -= colStride;
	}

	// Now handle the general case
	// A row of zeros must be in effect added to the top of each old column (in each loop)
	double* colLastPtr = x + NumRows*NumCols - numToTransform - 1; 
	for ( i = numXforms-1; i>0; i-- ) {
		numToTransform++;							// Number of non-trivial applications of this Householder transformation
		hBase -= hDiagStride;						// Pointer to the beginning of the Householder transformation
		colPtr = colLastPtr;
		for ( long j = numToTransform-1; j>0; j-- ) {
			// Get dot product
			double dotProd2N = -2.0*DotArray( numToTransform-1, hBase+colStride, colStride, colPtr+1, 1 );
			*colPtr = dotProd2N*(*hBase);			// Adding onto zero at initial point
			AddArrayScale( numToTransform-1, hBase+colStride, colStride, colPtr+1, 1, dotProd2N ); 
			colPtr -= NumRows;
		}
		// Do last one as a special case (may overwrite the Householder vector)
		CopyArrayScale( numToTransform, hBase, colStride, colPtr, 1, -2.0*(*hBase) );
		*colPtr += 1.0;				// Add back one one as identity
		// Done with this Householder transformation
		colLastPtr --;
	}    

	if ( numZerosSkipped!=0 ) {
		assert( numZerosSkipped==1 );
		// Fill first row and column with identity (More generally: first numZerosSkipped many rows and columns)
		double* d = x;
		*d = 1;
		double* d2 = d;
		for ( i=NumRows-1; i>0; i-- ) {
			*(++d) = 0;
			*(d2+=NumRows) = 0;
		}
	}
}

// **************** ConvertBidiagToDiagonal ***********************************************
// Do the iterative transformation from bidiagonal form to diagonal form using
//		Givens transformation.  (Golub-Reinsch)
// U and V are square.  Size of U less than or equal to that of U.
void MatrixRmn::ConvertBidiagToDiagonal( MatrixRmn& U, MatrixRmn& V, VectorRn& w, VectorRn& superDiag ) const
{
	// These two index into the last bidiagonal block  (last in the matrix, it will be
	//	first one handled.
	long lastBidiagIdx = V.NumRows-1;
	long firstBidiagIdx = 0;
	double eps = 1.0e-15 * Max(w.MaxAbs(), superDiag.MaxAbs());

	while ( true ) {
		bool workLeft = UpdateBidiagIndices( &firstBidiagIdx, &lastBidiagIdx, w, superDiag, eps );
		if ( !workLeft ) {
			break;
		}

		// Get ready for first Givens rotation
		// Push non-zero to M[2,1] with Givens transformation
		double* wPtr = w.x+firstBidiagIdx;
		double* sdPtr = superDiag.x+firstBidiagIdx;
		double extraOffDiag=0.0;
		if ( (*wPtr)==0.0 ) {
			ClearRowWithDiagonalZero( firstBidiagIdx, lastBidiagIdx, U, wPtr, sdPtr, eps );
			if ( firstBidiagIdx>0 ) {
				if ( NearZero( *(--sdPtr), eps ) ) {
					*sdPtr = 0.0;
				}
				else {
					ClearColumnWithDiagonalZero( firstBidiagIdx, V, wPtr, sdPtr, eps );
				}
			}
			continue;
		}

		// Estimate an eigenvalue from bottom four entries of M
		// This gives a lambda value which will shift the Givens rotations
		// Last four entries of M^T * M are  ( ( A, B ), ( B, C ) ).
		double A;
		A = (firstBidiagIdx<lastBidiagIdx-1) ? Square(superDiag[lastBidiagIdx-2]): 0.0;
		double BSq = Square(w[lastBidiagIdx-1]);
		A += BSq;										// The "A" entry of M^T * M
		double C = Square(superDiag[lastBidiagIdx-1]);	
		BSq *= C;										// The squared "B" entry
		C += Square(w[lastBidiagIdx]);					// The "C" entry
		double lambda;									// lambda will hold the estimated eigenvalue
		lambda = sqrt( Square((A-C)*0.5) + BSq );		// Use the lambda value that is closest to C.
		if ( A > C ) {
			lambda = -lambda;
		}
		lambda += (A+C)*0.5;						// Now lambda equals the estimate for the last eigenvalue
		double t11 = Square(w[firstBidiagIdx]);
		double t12 = w[firstBidiagIdx]*superDiag[firstBidiagIdx];

		double c, s;
		CalcGivensValues( t11-lambda, t12, &c, &s );
		ApplyGivensCBTD( c, s, wPtr, sdPtr, &extraOffDiag, wPtr+1 );
		V.PostApplyGivens( c, -s, firstBidiagIdx );
		long i;
		for ( i=firstBidiagIdx; i<lastBidiagIdx-1; i++ ) {
			// Push non-zero from M[i+1,i] to M[i,i+2] 
			CalcGivensValues( *wPtr, extraOffDiag, &c, &s );
			ApplyGivensCBTD( c, s, wPtr, sdPtr, &extraOffDiag, extraOffDiag, wPtr+1, sdPtr+1 );
			U.PostApplyGivens( c, -s, i );
			// Push non-zero from M[i,i+2] to M[1+2,i+1]
			CalcGivensValues( *sdPtr, extraOffDiag, &c, &s );
			ApplyGivensCBTD( c, s, sdPtr, wPtr+1, &extraOffDiag, extraOffDiag, sdPtr+1, wPtr+2 );
			V.PostApplyGivens( c, -s, i+1 );
			wPtr++;
			sdPtr++;
		}
		// Push non-zero value from M[i+1,i] to M[i,i+1] for i==lastBidiagIdx-1
		CalcGivensValues( *wPtr, extraOffDiag, &c, &s );
		ApplyGivensCBTD( c, s, wPtr, &extraOffDiag, sdPtr, wPtr+1 );
		U.PostApplyGivens( c, -s, i );

		// The next line for debugging only
		// DebugCalcBidiagCheck( V, w, superDiag, U );
	}
}

// This is called when there is a zero diagonal entry, with a non-zero superdiagonal entry on the same row.
// We use Givens rotations to "chase" the non-zero entry across the row; when it reaches the last
//	column, it is finally zeroed away.
// wPtr points to the zero entry on the diagonal.  sdPtr points to the non-zero superdiagonal entry on the same row.
void MatrixRmn::ClearRowWithDiagonalZero( long firstBidiagIdx, long lastBidiagIdx, MatrixRmn& U, double *wPtr, double *sdPtr, double eps )
{
	double curSd = *sdPtr;		// Value being chased across the row
	*sdPtr = 0.0;
	long i=firstBidiagIdx+1;
	while (true) { 
		// Rotate row i and row firstBidiagIdx (Givens rotation)
		double c, s;
		CalcGivensValues( *(++wPtr), curSd, &c, &s );
		U.PostApplyGivens( c, -s, i, firstBidiagIdx );
		*wPtr = c*(*wPtr) - s*curSd;
		if ( i==lastBidiagIdx ) {
			break;
		}
		curSd = s*(*(++sdPtr));		// New value pops up one column over to the right
		*sdPtr = c*(*sdPtr);
		i++;
	}
}

// This is called when there is a zero diagonal entry, with a non-zero superdiagonal entry in the same column.
// We use Givens rotations to "chase" the non-zero entry up the column; when it reaches the last
//	column, it is finally zeroed away.
// wPtr points to the zero entry on the diagonal.  sdPtr points to the non-zero superdiagonal entry in the same column.
void MatrixRmn::ClearColumnWithDiagonalZero( long endIdx, MatrixRmn& V, double *wPtr, double *sdPtr, double eps )
{
	double curSd = *sdPtr;		// Value being chased up the column
	*sdPtr = 0.0;
	long i = endIdx-1;
	while ( true ) {
		double c, s;
		CalcGivensValues( *(--wPtr), curSd, &c, &s );
		V.PostApplyGivens( c, -s, i, endIdx );
		*wPtr = c*(*wPtr) - s*curSd;
		if ( i==0 ) {
			break;
		}
		curSd = s*(*(--sdPtr));		// New value pops up one row above
		if ( NearZero( curSd, eps ) ) {
			break;
		}
		*sdPtr = c*(*sdPtr);
		i--;
	}
}



// Matrix A is  ( ( a c ) ( b d ) ), i.e., given in column order.
// Mult's G[c,s]  times  A, replaces A.
void MatrixRmn::ApplyGivensCBTD( double cosine, double sine, double *a, double *b, double *c, double *d )
{
	double temp = *a;
	*a = cosine*(*a) - sine*(*b);
	*b = sine*temp + cosine*(*b);
	temp = *c;
	*c = cosine*(*c) - sine*(*d);
	*d = sine*temp + cosine*(*d);
}

// Now matrix A given in row order, A = ( ( a b c ) ( d e f ) ).
// Return G[c,s] * A, replace A.  d becomes zero, no need to return.
//  Also, it is certain the old *c value is taken to be zero!
void MatrixRmn::ApplyGivensCBTD( double cosine, double sine, double *a, double *b, double *c, 
															 double  d, double *e, double *f )
{
	*a = cosine*(*a) - sine*d;
	double temp = *b;
	*b = cosine*(*b) - sine*(*e);
	*e = sine*temp + cosine*(*e);
	*c =  -sine*(*f);
	*f =  cosine*(*f);
}

// Helper routine for SVD conversion from bidiagonal to diagonal
bool MatrixRmn::UpdateBidiagIndices( long *firstBidiagIdx, long *lastBidiagIdx, VectorRn& w, VectorRn& superDiag, double eps )
{
	long lastIdx = *lastBidiagIdx;
	double* sdPtr = superDiag.GetPtr( lastIdx-1 );		// Entry above the last diagonal entry
	while ( NearZero(*sdPtr, eps) ) {				
		*(sdPtr--) = 0.0;
		lastIdx--;
		if ( lastIdx == 0 ) {
			return false;
		}
	}
	*lastBidiagIdx = lastIdx;
	long firstIdx = lastIdx-1;
	double* wPtr = w.GetPtr( firstIdx );
	while ( firstIdx > 0 ) {
		if ( NearZero( *wPtr, eps ) ) {			// If this diagonal entry (near) zero
			*wPtr = 0.0;
			break;
		}
		if ( NearZero(*(--sdPtr), eps) ) {		// If the entry above the diagonal entry is (near) zero
			*sdPtr = 0.0;
			break;
		}
		wPtr--;
		firstIdx--;
	}
	*firstBidiagIdx = firstIdx;
	return true;
}


// ******************************************DEBUG STUFFF

bool MatrixRmn::DebugCheckSVD( const MatrixRmn& U, const VectorRn& w, const MatrixRmn& V ) const
{
	// Special SVD test code

	MatrixRmn IV( V.GetNumRows(), V.GetNumColumns() );
	IV.SetIdentity();
	MatrixRmn VTV( V.GetNumRows(), V.GetNumColumns() );
	MatrixRmn::TransposeMultiply( V, V, VTV );
	IV -= VTV;
	double error = IV.FrobeniusNorm();

	MatrixRmn IU( U.GetNumRows(), U.GetNumColumns() );
	IU.SetIdentity();
	MatrixRmn UTU( U.GetNumRows(), U.GetNumColumns() );
	MatrixRmn::TransposeMultiply( U, U, UTU );
	IU -= UTU;
	error += IU.FrobeniusNorm();

	MatrixRmn Diag( U.GetNumRows(), V.GetNumRows() );
	Diag.SetZero();
	Diag.SetDiagonalEntries( w );
	MatrixRmn B(U.GetNumRows(), V.GetNumRows() );
	MatrixRmn C(U.GetNumRows(), V.GetNumRows() );
	MatrixRmn::Multiply( U, Diag, B );
	MatrixRmn::MultiplyTranspose( B, V, C );
	C -= *this;
	error += C.FrobeniusNorm();

	bool ret = ( fabs(error)<=1.0e-13*w.MaxAbs() );
	assert ( ret );
	return ret;
}

bool MatrixRmn::DebugCalcBidiagCheck( const MatrixRmn& U, const VectorRn& w, const VectorRn& superDiag, const MatrixRmn& V ) const
{
	// Special SVD test code

	MatrixRmn IV( V.GetNumRows(), V.GetNumColumns() );
	IV.SetIdentity();
	MatrixRmn VTV( V.GetNumRows(), V.GetNumColumns() );
	MatrixRmn::TransposeMultiply( V, V, VTV );
	IV -= VTV;
	double error = IV.FrobeniusNorm();

	MatrixRmn IU( U.GetNumRows(), U.GetNumColumns() );
	IU.SetIdentity();
	MatrixRmn UTU( U.GetNumRows(), U.GetNumColumns() );
	MatrixRmn::TransposeMultiply( U, U, UTU );
	IU -= UTU;
	error += IU.FrobeniusNorm();

	MatrixRmn DiagAndSuper( U.GetNumRows(), V.GetNumRows() );
	DiagAndSuper.SetZero();
	DiagAndSuper.SetDiagonalEntries( w );
	if ( this->GetNumRows()>=this->GetNumColumns() ) {
		DiagAndSuper.SetSequence( superDiag, 0, 1, 1, 1 );
	}
	else {
		DiagAndSuper.SetSequence( superDiag, 1, 0, 1, 1 );
	}
	MatrixRmn B(U.GetNumRows(), V.GetNumRows() );
	MatrixRmn C(U.GetNumRows(), V.GetNumRows() );
	MatrixRmn::Multiply( U, DiagAndSuper, B );
	MatrixRmn::MultiplyTranspose( B, V, C );
	C -= *this;
	error += C.FrobeniusNorm();

	bool ret = ( fabs(error)<1.0e-13*Max(w.MaxAbs(),superDiag.MaxAbs()) );
	assert ( ret );
	return ret;
}

