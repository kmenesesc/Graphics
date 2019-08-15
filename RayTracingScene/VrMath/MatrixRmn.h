/*
 *
 * RayTrace Software Package, release 3.0.  May 3, 2006
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
// MatrixRmn:  Matrix over reals  (Variable dimensional matrices)
//
//    Not very sophisticated yet.  Needs more functionality
//		To do: better handling of resizing.
//

#ifndef MATRIX_RMN_H
#define MATRIX_RMN_H

#include <math.h>
#include <assert.h>
#include "LinearR3.h"
#include "VectorRn.h"
class MatrixLowerRnn;

class MatrixRmn {

public:
	MatrixRmn();								// Null constructor
	MatrixRmn( long numRows, long numCols );	// Constructor with length
	MatrixRmn( const MatrixRmn& init );			// Set equal to another matrix
	~MatrixRmn();								// Destructor

	void SetSize( long numRows, long numCols );
	long GetNumRows() const { return NumRows; } 
	long GetNumColumns() const { return NumCols; } 

	// Return entry in row i and column j.
	double Get( long i, long j ) const;
	void GetTriple( long i, long j, VectorR3 *retValue ) const;
	void GetColumn( long i, VectorRn* retVector ) const;
	VectorRn GetColumn( long i) const { VectorRn ret(NumRows); GetColumn(i,&ret); return ret; };

	// Use GetPtr to get pointer into the array (efficient)
	// Use GetEntry to a reference of an element in the array (same efficiency)
	// They are friendly in that anyone can change the array contents (be careful!)
	// The entries are in column order!!!
	// Use this with care.  You may call GetRowStride and GetColStride to navigate
	//			within the matrix.  I do not expect these values to ever change.
	const double* GetPtr() const;
	double* GetPtr();
	const double* GetPtr( long i, long j ) const;
	double* GetPtr( long i, long j );
	const double& GetEntry( long i, long j ) const;
	double& GetEntry( long i, long j );
	const double* GetColumnPtr( long j ) const;
	double* GetColumnPtr( long j );
	const double* GetRowPtr( long i ) const;
	double* GetRowPtr( long i );
	long GetRowStride() const { return NumRows; }		// Step size (stride) along a row
	long GetColStride() const { return 1; }				// Step size (stide) along a column

	void Set( long i, long j, double val );
	void SetTriple( long i, long c, const VectorR3& u );

	void SetZero(); 
	void SetIdentity();
	void Set( const MatrixRmn& B );				// The two matrices MUST already be the same size!
	void Set( const MatrixLowerRnn& lowerTri );	// Sets equal to a lower triangular matrix
	void SetByColumns( const double* values );		// Assumes size already set correctly -- TO BE WRITTEN
	void SetByRows( const double* values );			// Assumes size already set correctly!!
	void SetByColumns( long nRows, long nCols, const double* values );	// To be written
	void SetByRows( long nRows, long rCols, const double* values );		// Resize and set by rows
	void SetTranspose( const MatrixRmn& M);			// Set equal to transpose of M.
	void SetDiagonalEntries( double d );
	void SetDiagonalEntries( const VectorRn& d );
	void SetSuperDiagonalEntries( double d );
	void SetSuperDiagonalEntries( const VectorRn& d );
	void SetSubDiagonalEntries( double d );
	void SetSubDiagonalEntries( const VectorRn& d );
	void SetColumn(long i, const VectorRn& d );
	void SetRow(long i, const VectorRn& d );
	void SetSequence( const VectorRn& d, long startRow, long startCol, long deltaRow, long deltaCol );

	// Loads matrix in as a sub-matrix.  (i,j) is the base point. Defaults to (0,0).
	// The "Tranpose" versions load the transpose of A.
	void LoadAsSubmatrix( const MatrixRmn& A );
	void LoadAsSubmatrix( long i, long j, const MatrixRmn& A );			// To be written
	void LoadAsSubmatrixTranspose( const MatrixRmn& A );				
	void LoadAsSubmatrixTranspose( long i, long j, const MatrixRmn& A );// To be written

	// Extract (copy out) submatrix
	// The size of "result" should already be set to the correct size!
	void GetSubmatrix( long baseRow, long baseColumn, MatrixRmn* result ) const;

	// Checking properties
	bool CheckDiagPositive() const;
	double MinDiagonalElement() const;
	double MinDiagonalElement( long *indexOfMin ) const;

	// Norms
	double FrobeniusNormSq() const;
	double FrobeniusNorm() const;
	double MaxRowSum() const;		// Max row sum of absolute values

	// Eigenvalues
	double GershgorinEigenvalueLowerBound() const;	// Returns lower bound on minimum eigenvalue (using rows)

	// Operations on VectorRn's
	void Multiply( const VectorRn& v, VectorRn& result ) const;					// result = (this)*(v)
	void MultiplyTranspose( const VectorRn& v, VectorRn& result ) const;		// Equivalent to mult by row vector on left
	double DotProductColumn( const VectorRn& v, long colNum ) const;			// Returns dot product of v with i-th column

	void MultiplyStrictUpperTriangular( const VectorRn& v, VectorRn& result ) const;	// Multiply by super diagonal part
	void MultiplyLowerTriangular( const VectorRn& v, VectorRn& result ) const;	// Multiply by lower diagonal part

	// Operations on MatrixRmn's
	MatrixRmn& operator= ( const MatrixRmn& B );
	MatrixRmn& operator*=( double );
	MatrixRmn& operator*=( const MatrixRmn& B );
	MatrixRmn& operator/=( double d ) { assert(d!=0.0); *this *= (1.0/d); return *this; }
	MatrixRmn& SetScaled( const MatrixRmn& B, double factor );
	MatrixRmn& AddScaled( const MatrixRmn& B, double factor );
	MatrixRmn& operator+=( const MatrixRmn& B );
	MatrixRmn& operator-=( const MatrixRmn& B );
	void Negate();
	static MatrixRmn& Multiply( const MatrixRmn& A, const MatrixRmn& B, MatrixRmn& dst );				// Sets dst = A*B.
	static MatrixRmn& MultiplyTranspose( const MatrixRmn& A, const MatrixRmn& B, MatrixRmn& dst );		// Sets dst = A*(B^tranpose).
	static MatrixRmn& TransposeMultiply( const MatrixRmn& A, const MatrixRmn& B, MatrixRmn& dst );		// Sets dst = (A^transpose)*B.
	// MultiplyTridiagonal: Sets dst = A*B. where B is tridiagonal
	static MatrixRmn& MultiplyTridiagonal( const MatrixRmn& A, const MatrixRmn& B, MatrixRmn& dst );


	// Miscellaneous operation
	MatrixRmn& AddToDiagonal( double d );					// Adds d to each diagonal element

	// Outer products
	MatrixRmn& SetOuterProduct( const VectorRn& v );
	MatrixRmn& AddOuterProduct( const VectorRn& v );

	// Solving systems of linear equations
	void Solve( const VectorRn& b, VectorRn* x ) const;	  // Solves the equation   (*this)*x = b;    Uses row operations.  Assumes *this is invertible.   

	// SolveLowerDiagonal: Assumes (*this) lower diagonal.  Solves (*this) x  = b;
	void SolveLowerDiagonal( const VectorRn& b, VectorRn* ret ) const; 

	// SolveLowerDiagonalTranspose: Assumes (*this) lower diagonal.  Solves (*this)^T x  = b;
	void SolveLowerDiagonalTranspose( const VectorRn& b, VectorRn* ret ) const; 

	// Inversion routines
	void InverseWithRowOps( MatrixRmn* inverse ) const;

	// Row Echelon Form and Reduced Row Echelon Form routines
	// Row echelon form here allows non-negative entries (instead of 1's) in the positions of lead variables.
	// Reduced row echelon form is same, but with only 1's as lead variables.
	void ConvertToRefNoFree();				// Converts the matrix in place to row echelon form -- assumption is no free variables will be found
	void ConvertToRrefNoFree();			// Converts the matrix in place to reduced row echelon form. Assumes no free variables found.
	// Row convert to identity is useful for computing inverses.
	void RowConvertToIdentity();		// Uses row operations to convert first square part to the identity.

	// Next two routines still to be written
	// void ConvertToRef( long numVars);		// Converts the matrix in place to row echelon form -- numVars is number of columns to work with.
	//void ConvertToRef( long numVars, double eps);		// Same, but eps is the measure of closeness to zero

	// Exponentiation - calculate exp( t*A ).
	void CalcExponential( MatrixRmn& retExpOfAt, double t=1.0 ) const;
	MatrixRmn Exponential( double t=1 ) const;

	// Cholesky decomposition of positive semidefinite, symmetric matrix
	//   Finds G such that  this ==  G * G^T, G lower diagonal.
	//   Returns "true" if successful.
	//   The second calling format returns a counter example if it returns false,
	//			namely, a vector v=negativeExample s.t. negativeAmount = n^T*this*n < 0.0;
	bool CalcCholesky( MatrixRmn& retCholesky ) const;
	bool CalcCholesky( MatrixRmn& retCholesky, VectorRn* nonPosExample, double* nonPosAmount ) const;
	// The "inner" version not for general use.  Returns a partially computed Cholesky conversion with column number
	//		at which failure occurs.
	long CalcCholeskyInner( MatrixRmn& retCholesky ) const;

	// Givens transformation
	static void CalcGivensValues( double a, double b, double *c, double *s );
	void PostApplyGivens( double c, double s, long idx );							// Applies Givens transform to columns idx and idx+1.
	void PostApplyGivens( double c, double s, long idx1, long idx2 );				// Applies Givens transform to columns idx1 and idx2.

	// Singular value decomposition
	void ComputeSVD( MatrixRmn& U, VectorRn& w, MatrixRmn& V ) const;
	// Good for debugging SVD computations (I recommend this be used for any new application to check for bugs/instability).
	bool DebugCheckSVD( const MatrixRmn& U, const VectorRn& w, const MatrixRmn& V ) const;

	// Some useful routines for experts who understand the inner workings of these classes.
	inline static double DotArray( long length, const double* ptrA, long strideA, const double* ptrB, long strideB );
	inline static void CopyArrayScale( long length, const double* from, long fromStride, double *to, long toStride, double scale );
	inline static void AddArrayScale( long length, const double* from, long fromStride, double *to, long toStride, double scale );

private:
	long NumRows;				// Number of rows 
	long NumCols;				// Number of columns 
	double *x;					// Array of vector entries - stored in column order
	long AllocSize;				// Allocated size of the x array

	static MatrixRmn WorkMatrix;	// Temporary work matrix
	static MatrixRmn WorkMatrix2;	// Temporary work matrix
	static MatrixRmn WorkMatrix3;	// Temporary work matrix
	static MatrixRmn& GetWorkMatrix() { return WorkMatrix; }
	static MatrixRmn& GetWorkMatrix(long numRows, long numCols) { WorkMatrix.SetSize( numRows, numCols ); return WorkMatrix; }

	// Internal helper routines for SVD calculations
	static void CalcBidiagonal( MatrixRmn& U, MatrixRmn& V, VectorRn& w, VectorRn& superDiag );
	void ConvertBidiagToDiagonal( MatrixRmn& U, MatrixRmn& V, VectorRn& w, VectorRn& superDiag ) const;
	static void SvdHouseholder( double* basePt,
							    long colLength, long numCols, long colStride, long rowStride, 
								double* retFirstEntry );
	void ExpandHouseholders( long numXforms, long numZerosSkipped, const double* basePt, long colStride, long rowStride );
	static bool UpdateBidiagIndices( long *firstDiagIdx, long *lastBidiagIdx, VectorRn& w, VectorRn& superDiag, double eps );
	static void ApplyGivensCBTD( double cosine, double sine, double *a, double *b, double *c, double *d );
	static void ApplyGivensCBTD( double cosine, double sine, double *a, double *b, double *c, 
															 double  d, double *e, double *f );
	static void ClearRowWithDiagonalZero( long firstBidiagIdx, long lastBidiagIdx, 
											MatrixRmn& U, double *wPtr, double *sdPtr, double eps );
	static void ClearColumnWithDiagonalZero( long endIdx, MatrixRmn& V, double *wPtr, double *sdPtr, double eps );
	bool DebugCalcBidiagCheck( const MatrixRmn& U, const VectorRn& w, const VectorRn& superDiag, const MatrixRmn& V ) const;
};

// *********************************************************
//  Some non-member functions
//**********************************************************
inline MatrixRmn operator+( const MatrixRmn& mat1, const MatrixRmn& mat2 );		// Returns the sum of two matrices
inline MatrixRmn operator-( const MatrixRmn& mat1, const MatrixRmn& mat2 );		// Returns the difference of two matrices

inline VectorRn operator*( const MatrixRmn& mat, const VectorRn& vec );		// Returns the product mat*vec.
inline MatrixRmn operator*(const MatrixRmn& A, const MatrixRmn& B );		// Returns the matrix product A*B

inline MatrixRmn operator*( const MatrixRmn& mat, double f );		// Returns the product mat*f.
inline MatrixRmn operator*( double f, const MatrixRmn& mat );		// Returns the product mat*f.
inline MatrixRmn operator/( const MatrixRmn& mat, double d );		// Returns the quotient mat/f.

// *********************************************************
// Inlined member functions
// *********************************************************

inline MatrixRmn::MatrixRmn() 
{
	NumRows = 0;
	NumCols = 0;
	x = 0;
	AllocSize = 0;
}

inline MatrixRmn::MatrixRmn( long numRows, long numCols ) 
{
	NumRows = 0;
	NumCols = 0;
	x = 0;
	AllocSize = 0;
	SetSize( numRows, numCols );
}

inline MatrixRmn::MatrixRmn( const MatrixRmn& init )			// Construct equal to another matrix
{
	NumRows = 0;
	NumCols = 0;
	x = 0;
	AllocSize = 0;
	if ( init.NumRows>0 && init.NumCols>0 ) {
		SetSize( init.NumRows, init.NumCols );
		Set( init );
	}
}

inline MatrixRmn::~MatrixRmn() 
{
	delete[] x;
}

// Resize.  
// If the array resize, the old data will no longer be at the
//		(i,j) entry in the table.
inline void MatrixRmn::SetSize( long numRows, long numCols )
{
	assert ( numRows>0 && numCols>0 );
	long newLength = numRows*numCols;
	if ( newLength>AllocSize ) {
		delete[] x;
		AllocSize = Max(newLength, AllocSize<<1);
		x = new double[AllocSize];
	}
	NumRows = numRows;
	NumCols = numCols;
} 

// Zero out the entire vector
inline void MatrixRmn::SetZero()
{
	double* target = x;
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		*(target++) = 0.0;
	}
}

inline void MatrixRmn::Negate()
{
	double* target = x;
	for ( long i=NumRows*NumCols; i>0; i--, target++ ) {
		*target = -(*target);
	}
}

// The two matrices MUST already be the same size!
inline void MatrixRmn::Set( const MatrixRmn& B )
{
	assert( NumRows==B.NumRows && NumCols==B.NumCols );
	double* toPtr = x;
	const double* fromPtr = B.x;
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		*(toPtr++) = *(fromPtr++);
	}
}


// Resize and load values by columns.
inline void MatrixRmn::SetByColumns( long nRows, long nCols, const double* values )
{
	SetSize( nRows, nCols );
	SetByColumns( values );
}

// Resize and load values by rows.
inline void MatrixRmn::SetByRows( long nRows, long nCols, const double* values )
{
	SetSize( nRows, nCols );
	SetByRows( values );
}

// Return entry in row i and column j.
inline double MatrixRmn::Get( long i, long j ) const { 
	assert ( i<NumRows && j<NumCols );	
	return *(x+j*NumRows+i); 
}

// Return a VectorR3 out of a column.  Starts at row 3*i, in column j.
inline void MatrixRmn::GetTriple( long i, long j, VectorR3 *retValue ) const
{
	long ii = 3*i;
	assert ( 0<=i && ii+2<NumRows && 0<=j && j<NumCols );
	retValue->Load( x+j*NumRows + ii );
}

// Get a pointer to the (0,0) entry.
// The entries are in column order.  
// This version gives read-only pointer
inline const double* MatrixRmn::GetPtr( ) const 
{ 
	return x; 
}

// Get a pointer to the (0,0) entry.
// The entries are in column order.  
inline double* MatrixRmn::GetPtr( ) 
{ 
	return x; 
}

// Get a pointer to the (i,j) entry.
// The entries are in column order.  
// This version gives read-only pointer
inline const double* MatrixRmn::GetPtr( long i, long j ) const 
{ 
	assert ( 0<=i && i<NumRows && 0<=j &&j<NumCols );	
	return (x+j*NumRows+i); 
}

// Get a pointer to the (i,j) entry.
// The entries are in column order.  
// This version gives pointer to writable data
inline double* MatrixRmn::GetPtr( long i, long j )
{ 
	assert ( i<NumRows && j<NumCols );	
	return (x+j*NumRows+i); 
}

// Get a pointer to the (i,j) entry.
// This version gives read-only pointer
// Same as GetPtr but returns a reference rather than a pointer
inline const double& MatrixRmn::GetEntry( long i, long j ) const
{
	assert ( 0<=i && i<NumRows && 0<=j &&j<NumCols );	
	return *((x+j*NumRows+i)); 
}

// Get a reference to the (i,j) entry.
// This version gives pointer to writable data
// Same as GetPtr but returns a reference rather than a pointer
inline double& MatrixRmn::GetEntry( long i, long j )
{
	assert ( 0<=i && i<NumRows && 0<=j &&j<NumCols );	
	return *((x+j*NumRows+i)); 
}


// Get a pointer to the j-th column.
// The entries are in column order.  
// This version gives read-only pointer
inline const double* MatrixRmn::GetColumnPtr( long j ) const
{ 
	assert ( 0<=j && j<NumCols );	
	return (x+j*NumRows); 
}

// Get a pointer to the j-th column.
// This version gives pointer to writable data
inline double* MatrixRmn::GetColumnPtr( long j )
{ 
	assert ( 0<=j && j<NumCols );	
	return (x+j*NumRows); 
}

/// Get a pointer to the i-th row
// The entries are in column order.  
// This version gives read-only pointer
inline const double* MatrixRmn::GetRowPtr( long i ) const
{ 
	assert ( 0<=i && i<NumRows );	
	return (x+i); 
}

// Get a pointer to the i-th row
// This version gives pointer to writable data
inline double* MatrixRmn::GetRowPtr( long i )
{ 
	assert ( 0<=i && i<NumRows );	
	return (x+i); 
}

// Set the (i,j) entry of the matrix
inline void MatrixRmn::Set( long i, long j, double val ) 
{ 
	assert( i<NumRows && j<NumCols );
	*(x+j*NumRows+i) = val;
}

// Set the i-th triple in the j-th column to u's three values
inline void MatrixRmn::SetTriple( long i, long j, const VectorR3& u ) 
{
	long ii = 3*i;
	assert ( 0<=i && ii+2<NumRows && 0<=j && j<NumCols );
	u.Dump( x+j*NumRows + ii );
}

// Set to be equal to the identity matrix
inline void MatrixRmn::SetIdentity()
{
	assert ( NumRows==NumCols );
	SetZero();
	SetDiagonalEntries(1.0);
}

inline MatrixRmn& MatrixRmn::operator= ( const MatrixRmn& B )
{
	SetSize( B.NumRows, B.NumCols );
	Set(B);
	return *this;
}

inline MatrixRmn& MatrixRmn::operator*=( double mult )
{
	double* aPtr = x;
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		(*(aPtr++)) *= mult;
	}
	return (*this);
}

inline MatrixRmn& MatrixRmn::operator*=( const MatrixRmn& B )
{
	WorkMatrix3 = *this;
	Multiply( WorkMatrix3, B, *this );
	return *this;
}

inline MatrixRmn& MatrixRmn::SetScaled( const MatrixRmn& B, double factor )
{
	assert (NumRows==B.NumRows && NumCols==B.NumCols);
	double* aPtr = x;
	double* bPtr = B.x;
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		(*(aPtr++)) = (*(bPtr++))*factor;
	}
	return (*this);
}

inline MatrixRmn& MatrixRmn::AddScaled( const MatrixRmn& B, double factor )
{
	assert (NumRows==B.NumRows && NumCols==B.NumCols);
	double* aPtr = x;
	double* bPtr = B.x;
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		(*(aPtr++)) += (*(bPtr++))*factor;
	}
	return (*this);
}

inline MatrixRmn& MatrixRmn::operator+=( const MatrixRmn& B )
{
	assert (NumRows==B.NumRows && NumCols==B.NumCols);
	double* aPtr = x;
	double* bPtr = B.x;
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		(*(aPtr++)) += *(bPtr++);
	}
	return (*this);
}

inline MatrixRmn& MatrixRmn::operator-=( const MatrixRmn& B )
{
	assert (NumRows==B.NumRows && NumCols==B.NumCols);
	double* aPtr = x;
	double* bPtr = B.x;
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		(*(aPtr++)) -= *(bPtr++);
	}
	return (*this);
}

inline double MatrixRmn::FrobeniusNormSq() const
{
	double* aPtr = x;
	double result = 0.0;
	for ( long i=NumRows*NumCols; i>0; i-- ) {
		result += Square( *(aPtr++) );
	}
	return result;
}

// Helper routine to calculate dot product
inline double MatrixRmn::DotArray( long length, const double* ptrA, long strideA, const double* ptrB, long strideB )
{
	double result = 0.0;
    for ( ; length>0 ; length-- ) {
		result += (*ptrA)*(*ptrB);
		ptrA += strideA;
		ptrB += strideB;
	}
	return result;
}

// Helper routine: copies and scales an array (src and dest may be equal, or overlap)
//    Copying is done in "forward" order.
inline void MatrixRmn::CopyArrayScale( long length, const double* from, long fromStride, double *to, long toStride, double scale )
{
	for ( ; length>0; length-- ) {
		*to = (*from)*scale;
		from += fromStride;
		to += toStride;
	}
}

// Helper routine: adds a scaled array  
//	fromArray = toArray*scale.
inline void MatrixRmn::AddArrayScale( long length, const double* from, long fromStride, double *to, long toStride, double scale )
{
	for ( ; length>0; length-- ) {
		*to += (*from)*scale;
		from += fromStride;
		to += toStride;
	}
}

// *********************************************************
// Inlined non-member functions
// *********************************************************

inline MatrixRmn operator+( const MatrixRmn& mat1, const MatrixRmn& mat2 )
{
	MatrixRmn sum(mat1);
	sum += mat2;
	return sum;
}

inline MatrixRmn operator-( const MatrixRmn& mat1, const MatrixRmn& mat2 )
{
	MatrixRmn diff(mat1);
	diff -= mat2;
	return diff;
}

// Returns the product mat*vec.
inline VectorRn operator*( const MatrixRmn& mat, const VectorRn& vec )
{
	VectorRn product( mat.GetNumRows() );
	mat.Multiply( vec, product );
	return product;
}

// Returns the product mat*f.
inline MatrixRmn operator*( const MatrixRmn& mat, double f )
{
	MatrixRmn ret(mat);
	ret *= f;
	return ret;
}

// Returns the product mat*f.
inline MatrixRmn operator*( double f, const MatrixRmn& mat )
{
	MatrixRmn ret(mat);
	ret *= f;
	return ret;
}

// Returns the quotient mat/d
inline MatrixRmn operator/( const MatrixRmn& mat, double d )
{
	MatrixRmn ret(mat);
	ret /= d;
	return ret;
}

// Returns the matrix product A*B.
inline MatrixRmn operator*(const MatrixRmn& A, const MatrixRmn& B )
{
	MatrixRmn ret(A.GetNumRows(),B.GetNumColumns());
	MatrixRmn::Multiply(A, B, ret);
	return ret;	
}

inline MatrixRmn MatrixRmn::Exponential( double t ) const
{
	MatrixRmn ret;
	CalcExponential( ret, t );
	return ret;
}





#endif // MATRIX_RMN_H