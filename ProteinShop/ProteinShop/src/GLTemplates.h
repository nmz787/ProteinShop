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
GLTemplates - Helper functions to simplify OpenGL programming.
***********************************************************************/

#ifndef GLTEMPLATES_INCLUDED
#define GLTEMPLATES_INCLUDED

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif


/******************************
Overloaded versions of glColor:
******************************/

inline void glColor(GLbyte r,GLbyte g,GLbyte b)
{
    glColor3b(r,g,b);
}

inline void glColor(GLubyte r,GLubyte g,GLubyte b)
{
    glColor3ub(r,g,b);
}

inline void glColor(GLshort r,GLshort g,GLshort b)
{
    glColor3s(r,g,b);
}

inline void glColor(GLushort r,GLushort g,GLushort b)
{
    glColor3us(r,g,b);
}

inline void glColor(GLint r,GLint g,GLint b)
{
    glColor3i(r,g,b);
}

inline void glColor(GLuint r,GLuint g,GLuint b)
{
    glColor3ui(r,g,b);
}

inline void glColor(GLfloat r,GLfloat g,GLfloat b)
{
    glColor3f(r,g,b);
}

inline void glColor(GLdouble r,GLdouble g,GLdouble b)
{
    glColor3d(r,g,b);
}

inline void glColor3(const GLbyte* components)
{
    glColor3bv(components);
}

inline void glColor3(const GLubyte* components)
{
    glColor3ubv(components);
}

inline void glColor3(const GLshort* components)
{
    glColor3sv(components);
}

inline void glColor3(const GLushort* components)
{
    glColor3usv(components);
}

inline void glColor3(const GLint* components)
{
    glColor3iv(components);
}

inline void glColor3(const GLuint* components)
{
    glColor3uiv(components);
}

inline void glColor3(const GLfloat* components)
{
    glColor3fv(components);
}

inline void glColor3(const GLdouble* components)
{
    glColor3dv(components);
}

inline void glColor(GLbyte r,GLbyte g,GLbyte b,GLbyte a)
{
    glColor4b(r,g,b,a);
}

inline void glColor(GLubyte r,GLubyte g,GLubyte b,GLubyte a)
{
    glColor4ub(r,g,b,a);
}

inline void glColor(GLshort r,GLshort g,GLshort b,GLshort a)
{
    glColor4s(r,g,b,a);
}

inline void glColor(GLushort r,GLushort g,GLushort b,GLushort a)
{
    glColor4us(r,g,b,a);
}

inline void glColor(GLint r,GLint g,GLint b,GLint a)
{
    glColor4i(r,g,b,a);
}

inline void glColor(GLuint r,GLuint g,GLuint b,GLuint a)
{
    glColor4ui(r,g,b,a);
}

inline void glColor(GLfloat r,GLfloat g,GLfloat b,GLfloat a)
{
    glColor4f(r,g,b,a);
}

inline void glColor(GLdouble r,GLdouble g,GLdouble b,GLdouble a)
{
    glColor4d(r,g,b,a);
}

inline void glColor4(const GLbyte* components)
{
    glColor4bv(components);
}

inline void glColor4(const GLubyte* components)
{
    glColor4ubv(components);
}

inline void glColor4(const GLshort* components)
{
    glColor4sv(components);
}

inline void glColor4(const GLushort* components)
{
    glColor4usv(components);
}

inline void glColor4(const GLint* components)
{
    glColor4iv(components);
}

inline void glColor4(const GLuint* components)
{
    glColor4uiv(components);
}

inline void glColor4(const GLfloat* components)
{
    glColor4fv(components);
}

inline void glColor4(const GLdouble* components)
{
    glColor4dv(components);
}

/*********************************
Overloaded versions of glTexCoord:
*********************************/

inline void glTexCoord(GLshort s)
{
    glTexCoord1s(s);
};

inline void glTexCoord(GLint s)
{
    glTexCoord1i(s);
};

inline void glTexCoord(GLfloat s)
{
    glTexCoord1f(s);
};

inline void glTexCoord(GLdouble s)
{
    glTexCoord1d(s);
};

inline void glTexCoord1(const GLshort* components)
{
    glTexCoord1sv(components);
}

inline void glTexCoord1(const GLint* components)
{
    glTexCoord1iv(components);
}

inline void glTexCoord1(const GLfloat* components)
{
    glTexCoord1fv(components);
}

inline void glTexCoord1(const GLdouble* components)
{
    glTexCoord1dv(components);
}

inline void glTexCoord(GLshort s,GLshort t)
{
    glTexCoord2s(s,t);
};

inline void glTexCoord(GLint s,GLint t)
{
    glTexCoord2i(s,t);
};

inline void glTexCoord(GLfloat s,GLfloat t)
{
    glTexCoord2f(s,t);
};

inline void glTexCoord(GLdouble s,GLdouble t)
{
    glTexCoord2d(s,t);
};

inline void glTexCoord2(const GLshort* components)
{
    glTexCoord2sv(components);
}

inline void glTexCoord2(const GLint* components)
{
    glTexCoord2iv(components);
}

inline void glTexCoord2(const GLfloat* components)
{
    glTexCoord2fv(components);
}

inline void glTexCoord2(const GLdouble* components)
{
    glTexCoord2dv(components);
}

inline void glTexCoord(GLshort s,GLshort t,GLshort r)
{
    glTexCoord3s(s,t,r);
};

inline void glTexCoord(GLint s,GLint t,GLint r)
{
    glTexCoord3i(s,t,r);
};

inline void glTexCoord(GLfloat s,GLfloat t,GLfloat r)
{
    glTexCoord3f(s,t,r);
};

inline void glTexCoord(GLdouble s,GLdouble t,GLdouble r)
{
    glTexCoord3d(s,t,r);
};

inline void glTexCoord3(const GLshort* components)
{
    glTexCoord3sv(components);
}

inline void glTexCoord3(const GLint* components)
{
    glTexCoord3iv(components);
}

inline void glTexCoord3(const GLfloat* components)
{
    glTexCoord3fv(components);
}

inline void glTexCoord3(const GLdouble* components)
{
    glTexCoord3dv(components);
}

inline void glTexCoord(GLshort s,GLshort t,GLshort r,GLshort q)
{
    glTexCoord4s(s,t,r,q);
};

inline void glTexCoord(GLint s,GLint t,GLint r,GLint q)
{
    glTexCoord4i(s,t,r,q);
};

inline void glTexCoord(GLfloat s,GLfloat t,GLfloat r,GLfloat q)
{
    glTexCoord4f(s,t,r,q);
};

inline void glTexCoord(GLdouble s,GLdouble t,GLdouble r,GLdouble q)
{
    glTexCoord4d(s,t,r,q);
};

inline void glTexCoord4(const GLshort* components)
{
    glTexCoord4sv(components);
}

inline void glTexCoord4(const GLint* components)
{
    glTexCoord4iv(components);
}

inline void glTexCoord4(const GLfloat* components)
{
    glTexCoord4fv(components);
}

inline void glTexCoord4(const GLdouble* components)
{
    glTexCoord4dv(components);
}

/*******************************
Overloaded versions of glNormal:
*******************************/

inline void glNormal(GLbyte x,GLbyte y,GLbyte z)
{
    glNormal3b(x,y,z);
};

inline void glNormal(GLshort x,GLshort y,GLshort z)
{
    glNormal3s(x,y,z);
};

inline void glNormal(GLint x,GLint y,GLint z)
{
    glNormal3i(x,y,z);
};

inline void glNormal(GLfloat x,GLfloat y,GLfloat z)
{
    glNormal3f(x,y,z);
};

inline void glNormal(GLdouble x,GLdouble y,GLdouble z)
{
    glNormal3d(x,y,z);
};

inline void glNormal3(const GLbyte* components)
{
    glNormal3bv(components);
}

inline void glNormal3(const GLshort* components)
{
    glNormal3sv(components);
}

inline void glNormal3(const GLint* components)
{
    glNormal3iv(components);
}

inline void glNormal3(const GLfloat* components)
{
    glNormal3fv(components);
}

inline void glNormal3(const GLdouble* components)
{
    glNormal3dv(components);
}

/*******************************
Overloaded versions of glVertex:
*******************************/

inline void glVertex(GLshort x,GLshort y)
{
    glVertex2s(x,y);
};

inline void glVertex(GLint x,GLint y)
{
    glVertex2i(x,y);
};

inline void glVertex(GLfloat x,GLfloat y)
{
    glVertex2f(x,y);
};

inline void glVertex(GLdouble x,GLdouble y)
{
    glVertex2d(x,y);
};

inline void glVertex2(const GLshort* components)
{
    glVertex2sv(components);
}

inline void glVertex2(const GLint* components)
{
    glVertex2iv(components);
}

inline void glVertex2(const GLfloat* components)
{
    glVertex2fv(components);
}

inline void glVertex2(const GLdouble* components)
{
    glVertex2dv(components);
}

inline void glVertex(GLshort x,GLshort y,GLshort z)
{
    glVertex3s(x,y,z);
};

inline void glVertex(GLint x,GLint y,GLint z)
{
    glVertex3i(x,y,z);
};

inline void glVertex(GLfloat x,GLfloat y,GLfloat z)
{
    glVertex3f(x,y,z);
};

inline void glVertex(GLdouble x,GLdouble y,GLdouble z)
{
    glVertex3d(x,y,z);
};

inline void glVertex3(const GLshort* components)
{
    glVertex3sv(components);
}

inline void glVertex3(const GLint* components)
{
    glVertex3iv(components);
}

inline void glVertex3(const GLfloat* components)
{
    glVertex3fv(components);
}

inline void glVertex3(const GLdouble* components)
{
    glVertex3dv(components);
}

inline void glVertex(GLshort x,GLshort y,GLshort z,GLshort w)
{
    glVertex4s(x,y,z,w);
};

inline void glVertex(GLint x,GLint y,GLint z,GLint w)
{
    glVertex4i(x,y,z,w);
};

inline void glVertex(GLfloat x,GLfloat y,GLfloat z,GLfloat w)
{
    glVertex4f(x,y,z,w);
};

inline void glVertex(GLdouble x,GLdouble y,GLdouble z,GLdouble w)
{
    glVertex4d(x,y,z,w);
};

inline void glVertex4(const GLshort* components)
{
    glVertex4sv(components);
}

inline void glVertex4(const GLint* components)
{
    glVertex4iv(components);
}

inline void glVertex4(const GLfloat* components)
{
    glVertex4fv(components);
}

inline void glVertex4(const GLdouble* components)
{
    glVertex4dv(components);
}

/**************************************
Overloaded versions of glColorPointer:
**************************************/

inline void glColorPointer(GLint numComponents,GLsizei stride,const GLbyte* pointer)
{
    glColorPointer(numComponents,GL_BYTE,stride,pointer);
}

inline void glColorPointer(GLint numComponents,GLsizei stride,const GLubyte* pointer)
{
    glColorPointer(numComponents,GL_UNSIGNED_BYTE,stride,pointer);
}

inline void glColorPointer(GLint numComponents,GLsizei stride,const GLshort* pointer)
{
    glColorPointer(numComponents,GL_SHORT,stride,pointer);
}

inline void glColorPointer(GLint numComponents,GLsizei stride,const GLushort* pointer)
{
    glColorPointer(numComponents,GL_UNSIGNED_SHORT,stride,pointer);
}

inline void glColorPointer(GLint numComponents,GLsizei stride,const GLint* pointer)
{
    glColorPointer(numComponents,GL_INT,stride,pointer);
}

inline void glColorPointer(GLint numComponents,GLsizei stride,const GLuint* pointer)
{
    glColorPointer(numComponents,GL_UNSIGNED_INT,stride,pointer);
}

inline void glColorPointer(GLint numComponents,GLsizei stride,const GLfloat* pointer)
{
    glColorPointer(numComponents,GL_FLOAT,stride,pointer);
}

inline void glColorPointer(GLint numComponents,GLsizei stride,const GLdouble* pointer)
{
    glColorPointer(numComponents,GL_DOUBLE,stride,pointer);
}

/****************************************
Overloaded versions of glTexCoordPointer:
****************************************/

inline void glTexCoordPointer(GLint numComponents,GLsizei stride,const GLshort* pointer)
{
    glTexCoordPointer(numComponents,GL_SHORT,stride,pointer);
}

inline void glTexCoordPointer(GLint numComponents,GLsizei stride,const GLint* pointer)
{
    glTexCoordPointer(numComponents,GL_INT,stride,pointer);
}

inline void glTexCoordPointer(GLint numComponents,GLsizei stride,const GLfloat* pointer)
{
    glTexCoordPointer(numComponents,GL_FLOAT,stride,pointer);
}

inline void glTexCoordPointer(GLint numComponents,GLsizei stride,const GLdouble* pointer)
{
    glTexCoordPointer(numComponents,GL_DOUBLE,stride,pointer);
}

/**************************************
Overloaded versions of glNormalPointer:
**************************************/

inline void glNormalPointer(GLsizei stride,const GLbyte* pointer)
{
    glNormalPointer(GL_BYTE,stride,pointer);
}

inline void glNormalPointer(GLsizei stride,const GLshort* pointer)
{
    glNormalPointer(GL_SHORT,stride,pointer);
}

inline void glNormalPointer(GLsizei stride,const GLint* pointer)
{
    glNormalPointer(GL_INT,stride,pointer);
}

inline void glNormalPointer(GLsizei stride,const GLfloat* pointer)
{
    glNormalPointer(GL_FLOAT,stride,pointer);
}

inline void glNormalPointer(GLsizei stride,const GLdouble* pointer)
{
    glNormalPointer(GL_DOUBLE,stride,pointer);
}

/**************************************
Overloaded versions of glVertexPointer:
**************************************/

inline void glVertexPointer(GLint numComponents,GLsizei stride,const GLshort* pointer)
{
    glVertexPointer(numComponents,GL_SHORT,stride,pointer);
}

inline void glVertexPointer(GLint numComponents,GLsizei stride,const GLint* pointer)
{
    glVertexPointer(numComponents,GL_INT,stride,pointer);
}

inline void glVertexPointer(GLint numComponents,GLsizei stride,const GLfloat* pointer)
{
    glVertexPointer(numComponents,GL_FLOAT,stride,pointer);
}

inline void glVertexPointer(GLint numComponents,GLsizei stride,const GLdouble* pointer)
{
    glVertexPointer(numComponents,GL_DOUBLE,stride,pointer);
}

/*********************************
Overloaded versions of glMaterial:
*********************************/

inline void glMaterial(GLenum face,GLenum pname,GLint param)
{
    glMateriali(face,pname,param);
}

inline void glMaterial(GLenum face,GLenum pname,GLfloat param)
{
    glMaterialf(face,pname,param);
}

inline void glMaterial(GLenum face,GLenum pname,const GLint* params)
{
    glMaterialiv(face,pname,params);
}

inline void glMaterial(GLenum face,GLenum pname,const GLfloat* params)
{
    glMaterialfv(face,pname,params);
}

/******************************
Overloaded versions of glLight:
******************************/

inline void glLight(GLenum light,GLenum pname,GLint param)
{
    glLighti(light,pname,param);
}

inline void glLight(GLenum light,GLenum pname,GLfloat param)
{
    glLightf(light,pname,param);
}

inline void glLight(GLenum light,GLenum pname,const GLint* params)
{
    glLightiv(light,pname,params);
}

inline void glLight(GLenum light,GLenum pname,const GLfloat* params)
{
    glLightfv(light,pname,params);
}

/*******************************
Overloaded versions of glTexEnv:
*******************************/

inline void glTexEnv(GLenum target,GLenum pname,GLint param)
{
    glTexEnvi(target,pname,param);
}

inline void glTexEnv(GLenum target,GLenum pname,GLfloat param)
{
    glTexEnvf(target,pname,param);
}

inline void glTexEnv(GLenum target,GLenum pname,const GLint* params)
{
    glTexEnviv(target,pname,params);
}

inline void glTexEnv(GLenum target,GLenum pname,const GLfloat* params)
{
    glTexEnvfv(target,pname,params);
}

/****************************
Overloaded versions of glFog:
****************************/

inline void glFog(GLenum pname,GLint param)
{
    glFogi(pname,param);
}

inline void glFog(GLenum pname,GLfloat param)
{
    glFogf(pname,param);
}

inline void glFog(GLenum pname,const GLint* params)
{
    glFogiv(pname,params);
}

inline void glFog(GLenum pname,const GLfloat* params)
{
    glFogfv(pname,params);
}

/***********************************
Overloaded versions of matrix calls:
***********************************/

inline void glTranslate(float x,float y,float z)
{
    glTranslatef(x,y,z);
}

inline void glTranslate(double x,double y,double z)
{
    glTranslated(x,y,z);
}

inline void glTranslate(const float t[3])
{
    glTranslatef(t[0],t[1],t[2]);
}

inline void glTranslate(const double t[3])
{
    glTranslated(t[0],t[1],t[2]);
}

inline void glRotate(float angle,float x,float y,float z)
{
    glRotatef(angle,x,y,z);
}

inline void glRotate(double angle,double x,double y,double z)
{
    glRotated(angle,x,y,z);
}

inline void glRotate(float angle,const float a[3])
{
    glRotatef(angle,a[0],a[1],a[2]);
}

inline void glRotate(double angle,const double a[3])
{
    glRotated(angle,a[0],a[1],a[2]);
}

inline void glScale(float scaleX,float scaleY,float scaleZ)
{
    glScalef(scaleX,scaleY,scaleZ);
}

inline void glScale(double scaleX,double scaleY,double scaleZ)
{
    glScaled(scaleX,scaleY,scaleZ);
}

inline void glScale(const float s[3])
{
    glScalef(s[0],s[1],s[2]);
}

inline void glScale(const double s[3])
{
    glScaled(s[0],s[1],s[2]);
}

inline void glLoadMatrix(const GLfloat* array)
{
    glLoadMatrixf(array);
}

inline void glLoadMatrix(const GLdouble* array)
{
    glLoadMatrixd(array);
}

inline void glMultMatrix(const GLfloat* array)
{
    glMultMatrixf(array);
}

inline void glMultMatrix(const GLdouble* array)
{
    glMultMatrixd(array);
}

/**********************************
Overloaded versions of glGet calls:
**********************************/

inline void glGet(GLenum pname,GLboolean* params)
{
    glGetBooleanv(pname,params);
}

inline void glGet(GLenum pname,GLint* params)
{
    glGetIntegerv(pname,params);
}

inline void glGet(GLenum pname,GLfloat* params)
{
    glGetFloatv(pname,params);
}

inline void glGet(GLenum pname,GLdouble* params)
{
    glGetDoublev(pname,params);
}

template <class ScalarParam>
inline ScalarParam glGet(GLenum pname) // glGet for a single value
{
    ScalarParam result;
    glGet(pname,&result);
    return result;
}

inline void glGetLight(GLenum light,GLenum pname,GLint* params)
{
    glGetLightiv(light,pname,params);
}

inline void glGetLight(GLenum light,GLenum pname,GLfloat* params)
{
    glGetLightfv(light,pname,params);
}

inline void glGetLight(GLenum light,GLenum pname,GLdouble* params)
{
    GLfloat fParams[4];
    glGetLightfv(light,pname,fParams);
    
    /* Determine the amount of data to copy based on the queried parameter: */
    int numComponents;
    switch(pname)
    {
        case GL_SPOT_EXPONENT:
        case GL_SPOT_CUTOFF:
        case GL_CONSTANT_ATTENUATION:
        case GL_LINEAR_ATTENUATION:
        case GL_QUADRATIC_ATTENUATION:
            numComponents=1;
            break;
        
        case GL_SPOT_DIRECTION:
            numComponents=3;
            break;
        
        case GL_AMBIENT:
        case GL_DIFFUSE:
        case GL_SPECULAR:
        case GL_POSITION:
            numComponents=4;
            break;
        
        default:
            numComponents=0;
    }
    for(int i=0;i<numComponents;++i)
        params[i]=GLdouble(fParams[i]);
}

template <class ScalarParam>
inline ScalarParam glGetLight(GLenum light,GLenum pname) // glGetLight for a single value
{
    ScalarParam result;
    glGetLight(light,pname,&result);
    return result;
}

#endif
