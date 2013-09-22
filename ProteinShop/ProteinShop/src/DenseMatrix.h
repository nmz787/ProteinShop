/***************************************************************************************
*cr								
*cr					Copyright (c) 2004, The Regents of the 
*cr	University of California, through Lawrence Berkeley National Laboratory, 
*cr	Univ. of Calif. at Davis, and Lawrence Livermore National Laboratory 
*cr	(subject to receipt of any required approvals from U.S. Dept. of Energy).  
*cr							All rights reserved.
*cr		
*cr		Please see the accompanying LICENSE file for further information.
*cr
***************************************************************************************/

/***********************************************************************
DenseMatrix - Class for dense matrices and their operations.
***********************************************************************/

#ifndef DENSEMATRIX_INCLUDED
#define DENSEMATRIX_INCLUDED

class DenseMatrix
{
    /* Elements: */
    private:
    int numRows,numColumns; // Size of matrix
    double** rows; // Array of pointers to rows
    double* columns; // Array of column entries
    /* Private methods: */
    void resize(int newNumRows,int newNumColumns); // Resizes a matrix, does not initialize the entries
    /* Constructors and destructors: */
    public:
    DenseMatrix(int sNumRows,int sNumColumns,const double* sEntries =0); // Create (empty) matrix
    DenseMatrix(const DenseMatrix& other); // Copy constructor
    DenseMatrix& operator=(const DenseMatrix& other); // Assignment operator
    ~DenseMatrix(void);
    /* Access methods: */
    int getNumRows(void) const // Returns number of rows
    {
        return numRows;
    };
    int getNumColumns(void) const // Returns number of columns
    {
        return numColumns;
    };
    const double& operator()(int i,int j) const // Returns matrix element at row i and column j
    {
        return rows[i][j];
    };
    double& operator()(int i,int j) // Ditto
    {
        return rows[i][j];
    };
    /* Conversion methods: */
    DenseMatrix getRow(int i) const; // Returns i-th row as (1,numColumns)-matrix
    DenseMatrix getColumn(int j) const; // Returns j-th column as (numRows,1)-matrix
    DenseMatrix& setRow(int i,const DenseMatrix& source); // Copies i-th row from (1,numColumns)-matrix
    DenseMatrix& setColumn(int j,const DenseMatrix& source); // Copies j-th column from (numRows,1)-matrix
    /* Algebra methods and operators: */
    friend DenseMatrix operator-(const DenseMatrix& matrix); // Negation operator (unary minus)
    friend DenseMatrix operator+(const DenseMatrix& matrix1,const DenseMatrix& matrix2); // Addition operator
    friend DenseMatrix operator-(const DenseMatrix& matrix1,const DenseMatrix& matrix2); // Subtraction operator
    friend DenseMatrix operator*(const DenseMatrix& matrix1,const DenseMatrix& matrix2); // Multiplication operator
    friend DenseMatrix& inplaceMultiplication(DenseMatrix& result,const DenseMatrix& matrix1,const DenseMatrix& matrix2); // Multiplies two matrices into existing result matrix
    friend DenseMatrix& inplaceTransposed1Multiplication(DenseMatrix& result,const DenseMatrix& matrix1,const DenseMatrix& matrix2); // Ditto, but the first matrix is transposed
    friend DenseMatrix& inplaceTransposed2Multiplication(DenseMatrix& result,const DenseMatrix& matrix1,const DenseMatrix& matrix2); // Ditto, but the second matrix is transposed
    DenseMatrix& operator*=(const DenseMatrix& matrix2); // Self-multipliciation operator. Performed in-place if matrix2 is numColumns x numColumns
    /* Other matrix methods: */
    DenseMatrix& zero(void); // Sets all entries to zero
    DenseMatrix transpose(void) const; // Returns flipped matrix
    double findColumnPivot(int start,int& pivotI) const; // Finds element of maximal absolute value in lower part of column
    double findFullPivot(int start,int& pivotI,int& pivotJ) const; // Finds element of maximal absolute value in lower right part matrix
    DenseMatrix& swapRows(int i1,int i2); // Swaps two rows
    DenseMatrix& swapColumns(int j1,int j2); // Swaps two columns
    DenseMatrix& scaleRow(int i,double factor); // Multiplies row by given factor
    DenseMatrix& combineRows(int destI,int sourceI,double factor); // Adds multiple of row sourceI to row destI
    double determinant(void) const; // Calculates square matrix' determinant
    int rank(void) const; // Returns rank of a square matrix
    DenseMatrix solveLinearEquations(const DenseMatrix& constants) const; // Solves set of linear equation systems
    void qr(DenseMatrix& q,DenseMatrix& r) const; // Computes the QR factorization of a matrix
};

#endif
