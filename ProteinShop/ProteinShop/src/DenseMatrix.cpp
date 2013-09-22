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

#include <assert.h>
#include <math.h>
#include <string.h>

#include "DenseMatrix.h"

/****************************
Methods of class DenseMatrix:
****************************/

void DenseMatrix::resize(int newNumRows,int newNumColumns)
{
    assert(newNumRows>0&&newNumColumns>0);
    if(newNumRows!=numRows||newNumColumns!=numColumns)
    {
        delete[] rows;
        delete[] columns;
        numRows=newNumRows;
        numColumns=newNumColumns;
        rows=new double*[numRows];
        columns=new double[numRows*numColumns];
        assert(rows!=NULL);
        assert(columns!=NULL);
        for(int i=0;i<numRows;++i)
            rows[i]=columns+numColumns*i;
    }
}

DenseMatrix::DenseMatrix(int sNumRows,int sNumColumns,const double* sEntries)
    :numRows(sNumRows),numColumns(sNumColumns),rows(new double*[numRows]),
     columns(new double[numRows*numColumns])
{
    assert(numRows>0&&numColumns>0);
    assert(rows!=NULL);
    assert(columns!=NULL);
    for(int i=0;i<numRows;++i)
        rows[i]=columns+numColumns*i;
    if(sEntries!=NULL)
        memcpy(columns,sEntries,numRows*numColumns*sizeof(double));
}

DenseMatrix::DenseMatrix(const DenseMatrix& other)
    :numRows(other.numRows),numColumns(other.numColumns),rows(new double*[numRows]),
     columns(new double[numRows*numColumns])
{
    assert(rows!=NULL);
    assert(columns!=NULL);
    for(int i=0;i<numRows;++i)
        rows[i]=columns+numColumns*i;
    memcpy(columns,other.columns,numRows*numColumns*sizeof(double));
}

DenseMatrix& DenseMatrix::operator=(const DenseMatrix& other)
{
    if(&other!=this)
    {
        delete[] rows;
        delete[] columns;
        numRows=other.numRows;
        numColumns=other.numColumns;
        rows=new double*[numRows];
        columns=new double[numRows*numColumns];
        assert(rows!=NULL);
        assert(columns!=NULL);
        for(int i=0;i<numRows;++i)
            rows[i]=columns+numColumns*i;
        memcpy(columns,other.columns,numRows*numColumns*sizeof(double));
    }
    return *this;
}

DenseMatrix::~DenseMatrix(void)
{
    delete[] rows;
    delete[] columns;
}

DenseMatrix DenseMatrix::getRow(int i) const
{
    assert(i>=0&&i<numRows);
    DenseMatrix result(1,numColumns);
    for(int j=0;j<numColumns;++j)
        result.rows[0][j]=rows[i][j];
    return result;
}

DenseMatrix DenseMatrix::getColumn(int j) const
{
    assert(j>=0&&j<numColumns);
    DenseMatrix result(numRows,1);
    for(int i=0;i<numRows;++i)
        result.rows[i][0]=rows[i][j];
    return result;
}

DenseMatrix& DenseMatrix::setRow(int i,const DenseMatrix& source)
{
    assert(i>=0&&i<numRows&&source.numRows==1&&source.numColumns==numColumns);
    for(int j=0;j<numColumns;++j)
        rows[i][j]=source.rows[0][j];
    return *this;
}

DenseMatrix& DenseMatrix::setColumn(int j,const DenseMatrix& source)
{
    assert(j>=0&&j<numColumns&&source.numRows==numRows&&source.numColumns==1);
    for(int i=0;i<numRows;++i)
        rows[i][j]=source.rows[i][0];
    return *this;
}

DenseMatrix operator-(const DenseMatrix& matrix)
{
    DenseMatrix result(matrix.numRows,matrix.numColumns);
    double* rPtr=result.columns;
    const double* ptr=matrix.columns;
    int numElements=result.numRows*result.numColumns;
    for(int i=0;i<numElements;++i,++rPtr,++ptr)
        *rPtr=-*ptr;
    return result;
}

DenseMatrix operator+(const DenseMatrix& matrix1,const DenseMatrix& matrix2)
{
    assert(matrix1.numRows==matrix2.numRows&&matrix1.numColumns==matrix2.numColumns);
    DenseMatrix result(matrix1.numRows,matrix1.numColumns);
    double* rPtr=result.columns;
    const double* ptr1=matrix1.columns;
    const double* ptr2=matrix2.columns;
    int numElements=result.numRows*result.numColumns;
    for(int i=0;i<numElements;++i,++rPtr,++ptr1,++ptr2)
        *rPtr=*ptr1+*ptr2;
    return result;
}

DenseMatrix operator-(const DenseMatrix& matrix1,const DenseMatrix& matrix2)
{
    assert(matrix1.numRows==matrix2.numRows&&matrix1.numColumns==matrix2.numColumns);
    DenseMatrix result(matrix1.numRows,matrix1.numColumns);
    double* rPtr=result.columns;
    const double* ptr1=matrix1.columns;
    const double* ptr2=matrix2.columns;
    int numElements=result.numRows*result.numColumns;
    for(int i=0;i<numElements;++i,++rPtr,++ptr1,++ptr2)
        *rPtr=*ptr1+*ptr2;
    return result;
}

DenseMatrix operator*(const DenseMatrix& matrix1,const DenseMatrix& matrix2)
{
    assert(matrix1.numColumns==matrix2.numRows);
    DenseMatrix result(matrix1.numRows,matrix2.numColumns);
    double* rPtr=result.columns;
    const double* rPtr1=matrix1.columns;
    for(int i=0;i<result.numRows;++i,rPtr1+=matrix1.numColumns)
        for(int j=0;j<result.numColumns;++j,++rPtr)
        {
            double sum=0.0;
            const double* ptr1=rPtr1;
            const double* ptr2=matrix2.columns+j;
            for(int k=0;k<matrix1.numColumns;++k,++ptr1,ptr2+=matrix2.numColumns)
                sum+=*ptr1**ptr2;
            *rPtr=sum;
        }
    return result;
}

DenseMatrix& inplaceMultiplication(DenseMatrix& result,const DenseMatrix& matrix1,const DenseMatrix& matrix2)
{
    assert(matrix1.numColumns==matrix2.numRows);
    result.resize(matrix1.numRows,matrix2.numColumns);
    double* rPtr=result.columns;
    const double* rPtr1=matrix1.columns;
    for(int i=0;i<result.numRows;++i,rPtr1+=matrix1.numColumns)
        for(int j=0;j<result.numColumns;++j,++rPtr)
        {
            double sum=0.0;
            const double* ptr1=rPtr1;
            const double* ptr2=matrix2.columns+j;
            for(int k=0;k<matrix1.numColumns;++k,++ptr1,ptr2+=matrix2.numColumns)
                sum+=*ptr1**ptr2;
            *rPtr=sum;
        }
    return result;
}

DenseMatrix& inplaceTransposed1Multiplication(DenseMatrix& result,const DenseMatrix& matrix1,const DenseMatrix& matrix2)
{
    assert(matrix1.numRows==matrix2.numRows);
    result.resize(matrix1.numColumns,matrix2.numColumns);
    double* rPtr=result.columns;
    for(int i=0;i<result.numRows;++i)
        for(int j=0;j<result.numColumns;++j,++rPtr)
        {
            double sum=0.0;
            const double* ptr1=matrix1.columns+i;
            const double* ptr2=matrix2.columns+j;
            for(int k=0;k<matrix1.numRows;++k,ptr1+=matrix1.numColumns,ptr2+=matrix2.numColumns)
                sum+=*ptr1**ptr2;
            *rPtr=sum;
        }
    return result;
}

DenseMatrix& inplaceTransposed2Multiplication(DenseMatrix& result,const DenseMatrix& matrix1,const DenseMatrix& matrix2)
{
    assert(matrix1.numColumns==matrix2.numColumns);
    result.resize(matrix1.numRows,matrix2.numRows);
    double* rPtr=result.columns;
    const double* rPtr1=matrix1.columns;
    for(int i=0;i<result.numRows;++i,rPtr1+=matrix1.numColumns)
    {
        const double* rPtr2=matrix2.columns;
        for(int j=0;j<result.numColumns;++j,rPtr2+=matrix2.numColumns,++rPtr)
        {
            double sum=0.0;
            const double* ptr1=rPtr1;
            const double* ptr2=rPtr2;
            for(int k=0;k<matrix1.numColumns;++k,++ptr1,++ptr2)
                sum+=*ptr1**ptr2;
            *rPtr=sum;
        }
    }
    return result;
}

DenseMatrix& DenseMatrix::operator*=(const DenseMatrix& matrix2)
{
    assert(numColumns==matrix2.numRows);
    int i,j,k;
    if(numColumns==matrix2.numColumns)
    {
        /* Compute matrix product in-place: */
        double* tempRow=new double[numColumns];
        assert(tempRow!=NULL);
        double* rowPtr1=columns;
        for(i=0;i<numRows;++i,rowPtr1+=numColumns)
        {
            for(j=0;j<numColumns;++j)
            {
                double sum=0.0;
                const double* ptr1=rowPtr1;
                const double* ptr2=matrix2.columns+j;
                for(k=0;k<numColumns;++k,++ptr1,ptr2+=numColumns)
                    sum+=*ptr1**ptr2;
                tempRow[j]=sum;
            }
            memcpy(rowPtr1,tempRow,numColumns*sizeof(double));
        }
        delete[] tempRow;
    }
    else
    {
        /* Compute matrix product into temporary array: */
        double* temp=new double[numRows*matrix2.numColumns];
        assert(temp!=NULL);
        double* rPtr=temp;
        double* rowPtr1=columns;
        for(i=0;i<numRows;++i,rowPtr1+=numColumns)
            for(j=0;j<matrix2.numColumns;++j,++rPtr)
            {
                double sum=0.0;
                const double* ptr1=rowPtr1;
                const double* ptr2=matrix2.columns+j;
                for(k=0;k<numColumns;++k,++ptr1,ptr2+=matrix2.numColumns)
                    sum+=*ptr1**ptr2;
                *rPtr=sum;
            }
        
        /* Change this matrix' size: */
        numColumns=matrix2.numColumns;
        rowPtr1=temp;
        for(i=0;i<numRows;++i,rowPtr1+=numColumns)
            rows[i]=rowPtr1;
        delete[] columns;
        columns=temp;
    }
    return* this;
}

DenseMatrix& DenseMatrix::zero(void)
{
    for(int i=0;i<numRows;++i)
        for(int j=0;j<numColumns;++j)
            rows[i][j]=0.0;
    return *this;
}

DenseMatrix DenseMatrix::transpose(void) const
{
    DenseMatrix result(numColumns,numRows);
    for(int i=0;i<numRows;++i)
        for(int j=0;j<numColumns;++j)
            result.rows[j][i]=rows[i][j];
    return result;
}

double DenseMatrix::findColumnPivot(int start,int& pivotI) const
{
    assert(start>=0&&start<numRows&&start<numColumns);
    double max=0.0;
    for(int i=start;i<numRows;++i)
        if(fabs(rows[i][start])>max)
        {
            max=fabs(rows[i][start]);
            pivotI=i;
        }
    return max;
}

double DenseMatrix::findFullPivot(int start,int& pivotI,int& pivotJ) const
{
    assert(start>=0&&start<numRows-1&&start<numColumns-1);
    double max=0.0;
    for(int i=start;i<numRows;++i)
        for(int j=start;j<numColumns;++j)
            if(fabs(rows[i][j])>max)
            {
                max=fabs(rows[i][j]);
                pivotI=i;
                pivotJ=j;
            }
    return max;
}

DenseMatrix& DenseMatrix::swapRows(int i1,int i2)
{
    assert(i1>=0&&i1<numRows&&i2>=0&&i2<numRows);
    for(int j=0;j<numColumns;++j)
    {
        double t=rows[i1][j];
        rows[i1][j]=rows[i2][j];
        rows[i2][j]=t;
    }
    return *this;
}

DenseMatrix& DenseMatrix::swapColumns(int j1,int j2)
{
    assert(j1>=0&&j1<numColumns&&j2>=0&&j2<numColumns);
    for(int i=0;i<numRows;++i)
    {
        double t=rows[i][j1];
        rows[i][j1]=rows[i][j2];
        rows[i][j2]=t;
    }
    return *this;
}

DenseMatrix& DenseMatrix::scaleRow(int i,double factor)
{
    assert(i>=0&&i<numRows);
    for(int j=0;j<numColumns;++j)
        rows[i][j]*=factor;
    return *this;
}

DenseMatrix& DenseMatrix::combineRows(int destI,int sourceI,double factor)
{
    assert(sourceI>=0&&sourceI<numRows&&destI>=0&&destI<numRows);
    for(int j=0;j<numColumns;++j)
        rows[destI][j]+=rows[sourceI][j]*factor;
    return *this;
}

double DenseMatrix::determinant(void) const
{
    assert(numRows==numColumns);
    DenseMatrix temp(*this);
    double result=1;
    int step;
    for(step=0;step<numRows-1;++step)
    {
        int pivotI,pivotJ;
        if(temp.findFullPivot(step,pivotI,pivotJ)==0.0)
        {
            result=0.0;
            break;
        }
        if(pivotI!=step)
        {
            temp.swapRows(step,pivotI);
            result=-result;
        }
        if(pivotJ!=step)
        {
            temp.swapRows(step,pivotJ);
            result=-result;
        }
        for(int i=step+1;i<numRows;++i)
            temp.combineRows(i,step,-temp.rows[i][step]/temp.rows[step][step]);
        result*=temp.rows[step][step];
    }
    result*=temp.rows[step][step];
    return result;
}

int DenseMatrix::rank(void) const
{
    int i,j;
    assert(numRows==numColumns);
    DenseMatrix temp=*this;
    int step;
    for(step=0;step<numRows;++step)
    {
        int pivotI;
        double max=temp.findColumnPivot(step,pivotI);
        if(max==0.0)
            break;
        if(pivotI!=step)
            temp.swapRows(step,pivotI);
        for(i=step+1;i<numRows;++i)
            temp.combineRows(i,step,-temp.rows[i][step]/temp.rows[step][step]);
        temp.scaleRow(step,1.0/temp.rows[step][step]);
    }
    return step;
}

DenseMatrix DenseMatrix::solveLinearEquations(const DenseMatrix& constants) const
{
    int i,j;
    assert(numRows==numColumns&&constants.numRows==numRows);
    DenseMatrix temp(numRows,numColumns+constants.numColumns);
    DenseMatrix result(numRows,constants.numColumns);
    for(i=0;i<numRows;++i)
    {
        for(j=0;j<numColumns;++j)
            temp.rows[i][j]=rows[i][j];
        for(j=0;j<constants.numColumns;++j)
            temp.rows[i][numColumns+j]=constants.rows[i][j];
    }
    int step;
    for(step=0;step<numRows-1;++step)
    {
        int pivotI;
        double max=temp.findColumnPivot(step,pivotI);
        assert(max>0.0);
        if(pivotI!=step)
            temp.swapRows(step,pivotI);
        for(i=step+1;i<numRows;++i)
            temp.combineRows(i,step,-temp.rows[i][step]/temp.rows[step][step]);
        temp.scaleRow(step,1.0/temp.rows[step][step]);
    }
    assert(temp.rows[numRows-1][numColumns-1]!=0.0);
    temp.scaleRow(numRows-1,1.0/temp.rows[numRows-1][numColumns-1]);
    for(step=numRows-1;step>=1;--step)
        for(i=step-1;i>=0;--i)
            temp.combineRows(i,step,-temp.rows[i][step]/temp.rows[step][step]);
    for(i=0;i<numRows;++i)
        for(j=0;j<result.numColumns;++j)
            result.rows[i][j]=temp.rows[i][numColumns+j];
    return result;
}

void DenseMatrix::qr(DenseMatrix& q,DenseMatrix& r) const
{
}
