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
GLGeometry - Wrappers to use Templatized Geometry Library classes in
combination with OpenGL functions.
***********************************************************************/

#ifndef GLGEOMETRY_INCLUDED
#define GLGEOMETRY_INCLUDED

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <GLTypes.h>

#include <Geometry/ComponentArray.h>
#include <Geometry/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/HVector.h>

/*****************************************************************
Conversion routines between geometry points/vectors and GLVectors:
*****************************************************************/

template <class ScalarType>
inline GLVector<ScalarType,3> toGLVector(const Geometry::ComponentArray<ScalarType,3>& ca)
{
    return GLVector<ScalarType,3>(ca.getComponents());
}

template <class ScalarType>
inline GLVector<ScalarType,4> toGLVector(const Geometry::ComponentArray<ScalarType,4>& ca)
{
    return GLVector<ScalarType,4>(ca.getComponents());
}

template <class ScalarType>
inline GLVector<ScalarType,4> toGLVector(const Geometry::Vector<ScalarType,3>& v)
{
    return GLVector<ScalarType,4>(v[0],v[1],v[2],ScalarType(0));
}

template <class ScalarType>
inline GLVector<ScalarType,4> toGLVector (const Geometry::Point<ScalarType,3>& p)
{
    return GLVector<ScalarType,4>(p[0],p[1],p[2],ScalarType(1));
}

/*********************************
Overloaded versions of glTexCoord:
*********************************/

template <class ScalarType>
inline void glTexCoord(const Geometry::ComponentArray<ScalarType,1>& tc)
{
    glTexCoord1(tc.getComponents());
}

template <class ScalarType>
inline void glTexCoord(const Geometry::ComponentArray<ScalarType,2>& tc)
{
    glTexCoord2(tc.getComponents());
}

template <class ScalarType>
inline void glTexCoord(const Geometry::ComponentArray<ScalarType,3>& tc)
{
    glTexCoord3(tc.getComponents());
}

template <class ScalarType>
inline void glTexCoord(const Geometry::ComponentArray<ScalarType,4>& tc)
{
    glTexCoord4(tc.getComponents());
}

/****************************************
Overloaded versions of glTexCoordPointer:
****************************************/

template <class ScalarType,int numComponents>
inline void glTexCoordPointer(GLsizei stride,const Geometry::ComponentArray<ScalarType,numComponents>* pointer)
{
    glTexCoordPointer(numComponents,stride,pointer[0].getComponents());
}

/*******************************
Overloaded versions of glNormal:
*******************************/

template <class ScalarType>
inline void glNormal(const Geometry::Vector<ScalarType,3>& n)
{
    glNormal3(n.getComponents());
}

/**************************************
Overloaded versions of glNormalPointer:
**************************************/

template <class ScalarType>
inline void glNormalPointer(GLsizei stride,const Geometry::Vector<ScalarType,3>* pointer)
{
    glNormalPointer(stride,pointer[0].getComponents());
}

/*******************************
Overloaded versions of glVertex:
*******************************/

template <class ScalarType>
inline void glVertex(const Geometry::Point<ScalarType,2>& v)
{
    glVertex2(v.getComponents());
}

template <class ScalarType>
inline void glVertex(const Geometry::Point<ScalarType,3>& v)
{
    glVertex3(v.getComponents());
}

template <class ScalarType>
inline void glVertex(const Geometry::Point<ScalarType,4>& v)
{
    glVertex4(v.getComponents());
}

/**************************************
Overloaded versions of glVertexPointer:
**************************************/

template <class ScalarType,int numComponents>
inline void glVertexPointer(GLsizei stride,const Geometry::Point<ScalarType,numComponents>* pointer)
{
    glVertexPointer(numComponents,stride,pointer[0].getComponents());
}

/************************************
Overloaded versions of glLight calls:
************************************/

inline void glLight(GLenum light,GLenum param,const Geometry::ComponentArray<GLint,4>& ca)
{
    glLightiv(light,param,ca.getComponents());
}

inline void glLight(GLenum light,GLenum param,const Geometry::ComponentArray<GLfloat,4>& ca)
{
    glLightfv(light,param,ca.getComponents());
}

/***********************************
Overloaded versions of matrix calls:
***********************************/

template <class ScalarParam>
inline void glTranslate(const Geometry::Vector<ScalarParam,3>& t)
{
    glTranslate(t.getComponents());
}

template <class ScalarParam>
inline void glRotate(ScalarParam angle,const Geometry::Vector<ScalarParam,3>& axis)
{
    glRotate(angle,axis.getComponents());
}

template <class ScalarParam>
inline void glScale(const Geometry::ComponentArray<ScalarParam,3>& s)
{
    glScale(s.getComponents());
}

#endif
