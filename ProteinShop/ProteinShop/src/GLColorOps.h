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
GLColorOps - Operations on RGBA color values.
***********************************************************************/

#ifndef GLCOLOROPS_INCLUDED
#define GLCOLOROPS_INCLUDED

#include <GLTypes.h>

template <class ScalarType,int numComponents>
GLColor<ScalarType,numComponents>& operator+=(GLColor<ScalarType,numComponents>& col1,const GLColor<ScalarType,numComponents>& col2)
{
    for(int i=0;i<numComponents;++i)
        col1[i]+=col2[i];
    
    return col1;
}

template <class ScalarType,int numComponents>
GLColor<ScalarType,numComponents> operator+(const GLColor<ScalarType,numComponents>& col1,const GLColor<ScalarType,numComponents>& col2)
{
    GLColor<ScalarType,numComponents> result;
    for(int i=0;i<numComponents;++i)
        result[i]=col1[i]+col2[i];
    
    return result;
}

template <class ScalarType,int numComponents>
GLColor<ScalarType,numComponents>& operator*=(GLColor<ScalarType,numComponents>& col,ScalarType factor)
{
    for(int i=0;i<numComponents;++i)
        col[i]*=factor;
    
    return col;
}

template <class ScalarType,int numComponents>
GLColor<ScalarType,numComponents> operator*(const GLColor<ScalarType,numComponents>& col,ScalarType factor)
{
    GLColor<ScalarType,numComponents> result;
    for(int i=0;i<numComponents;++i)
        result[i]=col[i]*factor;
    
    return result;
}

template <class ScalarType,int numComponents>
GLColor<ScalarType,numComponents> operator*(ScalarType factor,const GLColor<ScalarType,numComponents>& col)
{
    GLColor<ScalarType,numComponents> result;
    for(int i=0;i<numComponents;++i)
        result[i]=factor*col[i];
    
    return result;
}

template <class ScalarType,int numComponents>
GLColor<ScalarType,numComponents>& operator*=(GLColor<ScalarType,numComponents>& col1,const GLColor<ScalarType,numComponents>& col2)
{
    for(int i=0;i<numComponents;++i)
        col1[i]*=col2[i];
    
    return col1;
}

template <class ScalarType,int numComponents>
GLColor<ScalarType,numComponents> operator*(const GLColor<ScalarType,numComponents>& col1,const GLColor<ScalarType,numComponents>& col2)
{
    GLColor<ScalarType,numComponents> result;
    for(int i=0;i<numComponents;++i)
        result[i]=col1[i]*col2[i];
    
    return result;
}

template <int numComponents>
GLColor<GLfloat,numComponents>& clamp(GLColor<GLfloat,numComponents>& col)
{
    for(int i=0;i<numComponents;++i)
    {
        if(col[i]<0.0f)
            col[i]=0.0f;
        else if(col[i]>1.0f)
            col[i]=1.0f;
    }
    
    return col;
}

template <int numComponents>
GLColor<GLdouble,numComponents>& clamp(GLColor<GLdouble,numComponents>& col)
{
    for(int i=0;i<numComponents;++i)
    {
        if(col[i]<0.0)
            col[i]=0.0;
        else if(col[i]>1.0)
            col[i]=1.0;
    }
    
    return col;
}

#endif
