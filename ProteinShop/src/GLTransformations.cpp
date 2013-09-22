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
GLTransformations - Wrappers to use Templatized Geometry Library
transformation classes in combination with OpenGL transformation
functions.
***********************************************************************/

#define GLTRANSFORMATIONS_IMPLEMENTATION

#include <GLGeometry.h>

#include "GLTransformations.h"

/***********************************
Overloaded versions of matrix calls:
***********************************/

template <class ScalarParam>
void glRotate(const Geometry::Rotation<ScalarParam,3>& r)
{
    glRotate(Math::deg(r.getAngle()),r.getAxis());
}

template <class ScalarParam>
void glLoadMatrix(const Geometry::Rotation<ScalarParam,3>& r)
{
    glLoadIdentity();
    glRotate(Math::deg(r.getAngle()),r.getAxis());
}

template <class ScalarParam>
void glMultMatrix(const Geometry::Rotation<ScalarParam,3>& r)
{
    glRotate(Math::deg(r.getAngle()),r.getAxis());
}

template <class ScalarParam>
void glLoadMatrix(const Geometry::OrthonormalTransformation<ScalarParam,3>& t)
{
    glLoadIdentity();
    glTranslate(t.getTranslation());
    glRotate(t.getRotation());
}

template <class ScalarParam>
void glMultMatrix(const Geometry::OrthonormalTransformation<ScalarParam,3>& t)
{
    glTranslate(t.getTranslation());
    glRotate(t.getRotation());
}

template <class ScalarParam>
void glLoadMatrix(const Geometry::OrthogonalTransformation<ScalarParam,3>& t)
{
    glLoadIdentity();
    glTranslate(t.getTranslation());
    glRotate(t.getRotation());
    glScale(t.getScaling(),t.getScaling(),t.getScaling());
}

template <class ScalarParam>
void glMultMatrix(const Geometry::OrthogonalTransformation<ScalarParam,3>& t)
{
    glTranslate(t.getTranslation());
    glRotate(t.getRotation());
    glScale(t.getScaling(),t.getScaling(),t.getScaling());
}

template <class ScalarParam>
void glLoadMatrix(const Geometry::AffineTransformation<ScalarParam,3>& t)
{
    /* Copy the transformation coefficients into a temporary array: */
    ScalarParam temp[16];
    const typename Geometry::AffineTransformation<ScalarParam,3>::Matrix& m=t.getMatrix();
    ScalarParam* tPtr=temp;
    for(int j=0;j<4;++j)
    {
        for(int i=0;i<3;++i,++tPtr)
            *tPtr=m(i,j);
        *tPtr=ScalarParam(0);
        ++tPtr;
    }
    temp[15]=ScalarParam(1);
    
    /* Upload the temporary array: */
    glLoadMatrix(temp);
}

template <class ScalarParam>
void glMultMatrix(const Geometry::AffineTransformation<ScalarParam,3>& t)
{
    /* Copy the transformation coefficients into a temporary array: */
    ScalarParam temp[16];
    const typename Geometry::AffineTransformation<ScalarParam,3>::Matrix& m=t.getMatrix();
    ScalarParam* tPtr=temp;
    for(int j=0;j<4;++j)
    {
        for(int i=0;i<3;++i,++tPtr)
            *tPtr=m(i,j);
        *tPtr=ScalarParam(0);
        ++tPtr;
    }
    temp[15]=ScalarParam(1);
    
    /* Upload the temporary array: */
    glMultMatrix(temp);
}

template <class ScalarParam>
void glLoadMatrix(const Geometry::ProjectiveTransformation<ScalarParam,3>& t)
{
    /* Copy the transformation coefficients into a temporary array: */
    ScalarParam temp[16];
    const typename Geometry::ProjectiveTransformation<ScalarParam,3>::Matrix& m=t.getMatrix();
    ScalarParam* tPtr=temp;
    for(int j=0;j<4;++j)
        for(int i=0;i<4;++i,++tPtr)
            *tPtr=m(i,j);
    
    /* Upload the temporary array: */
    glLoadMatrix(temp);
}

template <class ScalarParam>
void glMultMatrix(const Geometry::ProjectiveTransformation<ScalarParam,3>& t)
{
    /* Copy the transformation coefficients into a temporary array: */
    ScalarParam temp[16];
    const typename Geometry::ProjectiveTransformation<ScalarParam,3>::Matrix& m=t.getMatrix();
    ScalarParam* tPtr=temp;
    for(int j=0;j<4;++j)
        for(int i=0;i<4;++i,++tPtr)
            *tPtr=m(i,j);
    
    /* Upload the temporary array: */
    glMultMatrix(temp);
}

template <class ScalarParam>
Geometry::ProjectiveTransformation<ScalarParam,3> glGetMatrix(GLenum whichMatrix)
{
    ScalarParam temp[16];
    glGet(whichMatrix,temp);
    return Geometry::ProjectiveTransformation<ScalarParam,3>::fromColumnMajor(temp);
}

#if !defined(NONSTANDARD_TEMPLATES)

/********************************************************
Force instantiation of all standard GLGeometry functions:
********************************************************/

template void glRotate(const Geometry::Rotation<float,3>&);
template void glRotate(const Geometry::Rotation<double,3>&);

template void glLoadMatrix(const Geometry::Rotation<float,3>&);
template void glLoadMatrix(const Geometry::Rotation<double,3>&);
template void glMultMatrix(const Geometry::Rotation<float,3>&);
template void glMultMatrix(const Geometry::Rotation<double,3>&);

template void glLoadMatrix(const Geometry::OrthonormalTransformation<float,3>&);
template void glLoadMatrix(const Geometry::OrthonormalTransformation<double,3>&);
template void glMultMatrix(const Geometry::OrthonormalTransformation<float,3>&);
template void glMultMatrix(const Geometry::OrthonormalTransformation<double,3>&);

template void glLoadMatrix(const Geometry::OrthogonalTransformation<float,3>&);
template void glLoadMatrix(const Geometry::OrthogonalTransformation<double,3>&);
template void glMultMatrix(const Geometry::OrthogonalTransformation<float,3>&);
template void glMultMatrix(const Geometry::OrthogonalTransformation<double,3>&);

template void glLoadMatrix(const Geometry::AffineTransformation<float,3>&);
template void glLoadMatrix(const Geometry::AffineTransformation<double,3>&);
template void glMultMatrix(const Geometry::AffineTransformation<float,3>&);
template void glMultMatrix(const Geometry::AffineTransformation<double,3>&);

template void glLoadMatrix(const Geometry::ProjectiveTransformation<float,3>&);
template void glLoadMatrix(const Geometry::ProjectiveTransformation<double,3>&);
template void glMultMatrix(const Geometry::ProjectiveTransformation<float,3>&);
template void glMultMatrix(const Geometry::ProjectiveTransformation<double,3>&);
template Geometry::ProjectiveTransformation<float,3> glGetMatrix<float>(GLenum whichMatrix);
template Geometry::ProjectiveTransformation<double,3> glGetMatrix<double>(GLenum whichMatrix);

#endif
