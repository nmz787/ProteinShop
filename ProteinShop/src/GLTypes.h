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
GLTypes - Helper classes to simplify OpenGL programming.
***********************************************************************/

#ifndef GLTYPES_INCLUDED
#define GLTYPES_INCLUDED

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLTemplates.h>

/***********************************************
GLColor: Generic class for colors in RGBA format
***********************************************/

template <class ScalarType =GLubyte,int numComponents =4>
class GLColor
{
};

/*********************************************
Specialized GLColor class without alpha value:
*********************************************/

template <class ScalarType>
class GLColor<ScalarType,3>
{
    /* Elements: */
    public:
    ScalarType rgba[3]; // RGB components
    
    /* Constructors and destructors: */
    GLColor(void)
    {
    };
    GLColor(ScalarType r,ScalarType g,ScalarType b)
    {
        rgba[0]=r;
        rgba[1]=g;
        rgba[2]=b;
    };
    template <class InputScalarType>
    GLColor(const InputScalarType sRgba[3])
    {
        for(int i=0;i<3;++i)
            rgba[i]=ScalarType(sRgba[i]);
    };
    template <class InputScalarType>
    GLColor(const GLColor<InputScalarType,3>& source) // Copy constructor with type conversion
    {
        for(int i=0;i<3;++i)
            rgba[i]=ScalarType(source.rgba[i]);
    };
    
    /* Methods: */
    const ScalarType* getComponents(void) const // Returns components as array
    {
        return rgba;
    };
    ScalarType* getComponents(void) // Ditto
    {
        return rgba;
    };
    ScalarType operator[](int index) const // Returns a single component
    {
        return rgba[index];
    };
    ScalarType& operator[](int index) // Ditto
    {
        return rgba[index];
    };
    friend void glColor(const GLColor& color) // Passes a color to OpenGL
    {
        glColor3(color.rgba);
    };
    friend void glColorPointer(GLsizei stride,const GLColor* colorArray) // Passes an array of colors to OpenGL
    {
        glColorPointer(3,stride,&colorArray[0].rgba[0]);
    };
};

/* Specialized versions of copy constructor: */

template <>
template <>
inline GLColor<GLfloat,3>::GLColor(const GLColor<GLbyte,3>& source)
{
    for(int i=0;i<3;++i)
        rgba[i]=GLfloat(source.rgba[i])/127.0;
}

template <>
template <>
inline GLColor<GLfloat,3>::GLColor(const GLColor<GLubyte,3>& source)
{
    for(int i=0;i<3;++i)
        rgba[i]=GLfloat(source.rgba[i])/255.0;
}

template <class ScalarType>
class GLColor<ScalarType,4>
{
    /* Elements: */
    public:
    ScalarType rgba[4]; // RGBA components
    
    /* Constructors and destructors: */
    GLColor(void)
    {
    };
    GLColor(ScalarType r,ScalarType g,ScalarType b)
    {
        rgba[0]=r;
        rgba[1]=g;
        rgba[2]=b;
        rgba[3]=ScalarType(1);
    };
    GLColor(ScalarType r,ScalarType g,ScalarType b,ScalarType a)
    {
        rgba[0]=r;
        rgba[1]=g;
        rgba[2]=b;
        rgba[3]=a;
    };
    template <class InputScalarType>
    GLColor(const InputScalarType sRgba[4])
    {
        for(int i=0;i<4;++i)
            rgba[i]=sRgba[i];
    };
    
    /* Methods: */
    const ScalarType* getComponents(void) const // Returns components as array
    {
        return rgba;
    };
    ScalarType* getComponents(void) // Ditto
    {
        return rgba;
    };
    ScalarType operator[](int index) const // Returns a single component
    {
        return rgba[index];
    };
    ScalarType& operator[](int index) // Ditto
    {
        return rgba[index];
    };
    friend void glColor(const GLColor& color) // Passes a color to OpenGL
    {
        glColor4(color.rgba);
    };
    friend void glColorPointer(GLsizei stride,const GLColor* colorArray) // Passes an array of colors to OpenGL
    {
        glColorPointer(4,stride,&colorArray[0].rgba[0]);
    };
};

/* Specialized versions of the GLColor constructor: */

template <>
inline GLColor<GLbyte,4>::GLColor(GLbyte r,GLbyte g,GLbyte b)
{
    rgba[0]=r;
    rgba[1]=g;
    rgba[2]=b;
    rgba[3]=GLbyte(127);
}

template <>
inline GLColor<GLubyte,4>::GLColor(GLubyte r,GLubyte g,GLubyte b)
{
    rgba[0]=r;
    rgba[1]=g;
    rgba[2]=b;
    rgba[3]=GLubyte(255U);
}

template <>
inline GLColor<GLshort,4>::GLColor(GLshort r,GLshort g,GLshort b)
{
    rgba[0]=r;
    rgba[1]=g;
    rgba[2]=b;
    rgba[3]=GLshort(32767);
}

template <>
inline GLColor<GLushort,4>::GLColor(GLushort r,GLushort g,GLushort b)
{
    rgba[0]=r;
    rgba[1]=g;
    rgba[2]=b;
    rgba[3]=GLushort(65535U);
}

template <>
inline GLColor<GLint,4>::GLColor(GLint r,GLint g,GLint b)
{
    rgba[0]=r;
    rgba[1]=g;
    rgba[2]=b;
    rgba[3]=GLint(2147483647);
}

template <>
inline GLColor<GLuint,4>::GLColor(GLuint r,GLuint g,GLuint b)
{
    rgba[0]=r;
    rgba[1]=g;
    rgba[2]=b;
    rgba[3]=GLuint(4294967295U);
}

/* Comparison operators: */

template <class ScalarType,int numComponents>
bool operator==(const GLColor<ScalarType,numComponents>& c1,const GLColor<ScalarType,numComponents>& c2)
{
    bool result=true;
    for(int i=0;i<numComponents&&result;++i)
        result=c1.rgba[i]==c2.rgba[i];
    return result;
}

template <class ScalarType,int numComponents>
bool operator!=(const GLColor<ScalarType,numComponents>& c1,const GLColor<ScalarType,numComponents>& c2)
{
    bool result=false;
    for(int i=0;i<numComponents&&!result;++i)
        result=c1.rgba[i]!=c2.rgba[i];
    return result;
}

/* Overloaded versions of glMaterial/glLight/glTexEnv/glFog calls: */

template <class ScalarParam>
inline GLColor<ScalarParam,4> glGetColor(GLenum pname)
{
    GLColor<ScalarParam,4> result;
    glGet(pname,result.rgba);
    return result;
}

template <class ScalarParam>
inline void glMaterial(GLenum face,GLenum pname,const GLColor<ScalarParam,4>& color)
{
    glMaterial(face,pname,color.rgba);
}

template <class ScalarParam>
inline void glLight(GLenum light,GLenum pname,const GLColor<ScalarParam,4>& color)
{
    glLight(light,pname,color.rgba);
}

template <class ScalarParam>
inline void glTexEnv(GLenum target,GLenum pname,const GLColor<ScalarParam,4>& color)
{
    glTexEnv(target,pname,color.rgba);
}

template <class ScalarParam>
inline void glFog(GLenum pname,const GLColor<ScalarParam,4>& color)
{
    glFog(pname,color.rgba);
}

/* Assorted functions using colors: */

template <class ScalarParam>
inline void glClearColor(const GLColor<ScalarParam,4>& color)
{
    glClearColor(color[0],color[1],color[2],color[3]);
}

/******************************************************************************************
GLVector: Class for generic vectors (texture coordinates, normal vectors, vertex positions)
******************************************************************************************/

template <class ScalarType =GLfloat,int numComponents =4>
class GLVector
{
};

template <class ScalarType>
class GLVector<ScalarType,1>
{
    /* Elements: */
    public:
    ScalarType xyzw[1]; // xyzw components
    
    /* Constructors and destructors: */
    GLVector(void)
    {
    };
    GLVector(ScalarType x)
    {
        xyzw[0]=x;
    };
    template <class InputScalarType>
    GLVector(const InputScalarType sXyzw[1])
    {
        for(int i=0;i<1;++i)
            xyzw[i]=sXyzw[i];
    };
    
    /* Methods: */
    const ScalarType* getComponents(void) const // Returns components as array
    {
        return xyzw;
    };
    ScalarType* getComponents(void) // Ditto
    {
        return xyzw;
    };
    ScalarType operator[](int index) const // Returns a single component
    {
        return xyzw[index];
    };
    ScalarType& operator[](int index) // Ditto
    {
        return xyzw[index];
    };
    friend void glTexCoord(const GLVector& texCoord) // Passes a texture coordinate to OpenGL
    {
        glTexCoord1(texCoord.xyzw);
    };
    friend void glTexCoordPointer(GLsizei stride,const GLVector* texCoordArray) // Passes an array of texture coordinates to OpenGL
    {
        glTexCoordPointer(1,stride,&texCoordArray[0].xyzw[0]);
    };
};

template <class ScalarType>
class GLVector<ScalarType,2>
{
    /* Elements: */
    public:
    ScalarType xyzw[2]; // xyzw components
    
    /* Constructors and destructors: */
    GLVector(void)
    {
    };
    GLVector(ScalarType x)
    {
        xyzw[0]=x;
        xyzw[1]=ScalarType(0);
    };
    GLVector(ScalarType x,ScalarType y)
    {
        xyzw[0]=x;
        xyzw[1]=y;
    };
    template <class InputScalarType>
    GLVector(const InputScalarType sXyzw[2])
    {
        for(int i=0;i<2;++i)
            xyzw[i]=sXyzw[i];
    };
    
    /* Methods: */
    const ScalarType* getComponents(void) const // Returns components as array
    {
        return xyzw;
    };
    ScalarType* getComponents(void) // Ditto
    {
        return xyzw;
    };
    ScalarType operator[](int index) const // Returns a single component
    {
        return xyzw[index];
    };
    ScalarType& operator[](int index) // Ditto
    {
        return xyzw[index];
    };
    friend void glTexCoord(const GLVector& texCoord) // Passes a texture coordinate to OpenGL
    {
        glTexCoord2(texCoord.xyzw);
    };
    friend void glTexCoordPointer(GLsizei stride,const GLVector* texCoordArray) // Passes an array of texture coordinates to OpenGL
    {
        glTexCoordPointer(2,stride,&texCoordArray[0].xyzw[0]);
    };
    friend void glVertex(const GLVector& vertex) // Passes a vertex position to OpenGL
    {
        glVertex2(vertex.xyzw);
    };
    friend void glVertexPointer(GLsizei stride,const GLVector* vertexArray) // Passes an array of vertex positions to OpenGL
    {
        glVertexPointer(2,stride,&vertexArray[0].xyzw[0]);
    };
};

template <class ScalarType>
class GLVector<ScalarType,3>
{
    /* Elements: */
    public:
    ScalarType xyzw[3]; // xyzw components
    
    /* Constructors and destructors: */
    GLVector(void)
    {
    };
    GLVector(ScalarType x)
    {
        xyzw[0]=x;
        xyzw[1]=ScalarType(0);
        xyzw[2]=ScalarType(0);
    };
    GLVector(ScalarType x,ScalarType y)
    {
        xyzw[0]=x;
        xyzw[1]=y;
        xyzw[2]=ScalarType(0);
    };
    GLVector(ScalarType x,ScalarType y,ScalarType z)
    {
        xyzw[0]=x;
        xyzw[1]=y;
        xyzw[2]=z;
    };
    template <class InputScalarType>
    GLVector(const InputScalarType sXyzw[3])
    {
        for(int i=0;i<3;++i)
            xyzw[i]=sXyzw[i];
    };
    
    /* Methods: */
    const ScalarType* getComponents(void) const // Returns components as array
    {
        return xyzw;
    };
    ScalarType* getComponents(void) // Ditto
    {
        return xyzw;
    };
    ScalarType operator[](int index) const // Returns a single component
    {
        return xyzw[index];
    };
    ScalarType& operator[](int index) // Ditto
    {
        return xyzw[index];
    };
    friend void glTexCoord(const GLVector& texCoord) // Passes a texture coordinate to OpenGL
    {
        glTexCoord3(texCoord.xyzw);
    };
    friend void glTexCoordPointer(GLsizei stride,const GLVector* texCoordArray) // Passes an array of texture coordinates to OpenGL
    {
        glTexCoordPointer(3,stride,&texCoordArray[0].xyzw[0]);
    };
    friend void glNormal(const GLVector& normal) // Passes a normal vector to OpenGL
    {
        glNormal3(normal.xyzw);
    };
    friend void glNormalPointer(GLsizei stride,const GLVector* normalArray) // Passes an array of normal vectors to OpenGL
    {
        glNormalPointer(stride,&normalArray[0].xyzw[0]);
    };
    friend void glVertex(const GLVector& vertex) // Passes a vertex position to OpenGL
    {
        glVertex3(vertex.xyzw);
    };
    friend void glVertexPointer(GLsizei stride,const GLVector* vertexArray) // Passes an array of vertex positions to OpenGL
    {
        glVertexPointer(3,stride,&vertexArray[0].xyzw[0]);
    };
};

template <class ScalarType>
class GLVector<ScalarType,4>
{
    /* Elements: */
    public:
    ScalarType xyzw[4]; // xyzw components
    
    /* Constructors and destructors: */
    GLVector(void)
    {
    };
    GLVector(ScalarType x)
    {
        xyzw[0]=x;
        xyzw[1]=ScalarType(0);
        xyzw[2]=ScalarType(0);
        xyzw[3]=ScalarType(1);
    };
    GLVector(ScalarType x,ScalarType y)
    {
        xyzw[0]=x;
        xyzw[1]=y;
        xyzw[2]=ScalarType(0);
        xyzw[3]=ScalarType(1);
    };
    GLVector(ScalarType x,ScalarType y,ScalarType z)
    {
        xyzw[0]=x;
        xyzw[1]=y;
        xyzw[2]=z;
        xyzw[3]=ScalarType(1);
    };
    GLVector(ScalarType x,ScalarType y,ScalarType z,ScalarType w)
    {
        xyzw[0]=x;
        xyzw[1]=y;
        xyzw[2]=z;
        xyzw[3]=w;
    };
    template <class InputScalarType>
    GLVector(const InputScalarType sXyzw[4])
    {
        for(int i=0;i<4;++i)
            xyzw[i]=sXyzw[i];
    };
    
    /* Methods: */
    const ScalarType* getComponents(void) const // Returns components as array
    {
        return xyzw;
    };
    ScalarType* getComponents(void) // Ditto
    {
        return xyzw;
    };
    ScalarType operator[](int index) const // Returns a single component
    {
        return xyzw[index];
    };
    ScalarType& operator[](int index) // Ditto
    {
        return xyzw[index];
    };
    friend void glTexCoord(const GLVector& texCoord) // Passes a texture coordinate to OpenGL
    {
        glTexCoord4(texCoord.xyzw);
    };
    friend void glTexCoordPointer(GLsizei stride,const GLVector* texCoordArray) // Passes an array of texture coordinates to OpenGL
    {
        glTexCoordPointer(4,stride,&texCoordArray[0].xyzw[0]);
    };
    friend void glVertex(const GLVector& vertex) // Passes a vertex position to OpenGL
    {
        glVertex4(vertex.xyzw);
    };
    friend void glVertexPointer(GLsizei stride,const GLVector* vertexArray) // Passes an array of vertex positions to OpenGL
    {
        glVertexPointer(4,stride,&vertexArray[0].xyzw[0]);
    };
};

/* Comparison operators: */

template <class ScalarType,int numComponents>
bool operator==(const GLVector<ScalarType,numComponents>& c1,const GLVector<ScalarType,numComponents>& c2)
{
    bool result=true;
    for(int i=0;i<numComponents&&result;++i)
        result=c1.xyzw[i]==c2.xyzw[i];
    return result;
}

template <class ScalarType,int numComponents>
bool operator!=(const GLVector<ScalarType,numComponents>& c1,const GLVector<ScalarType,numComponents>& c2)
{
    bool result=false;
    for(int i=0;i<numComponents&&!result;++i)
        result=c1.xyzw[i]!=c2.xyzw[i];
    return result;
}

/* Overloaded versions of glLight calls: */

template <class ScalarParam,int numComponents>
inline void glLight(GLenum light,GLenum pname,const GLVector<ScalarParam,numComponents>& vector)
{
    glLight(light,pname,vector.xyzw);
}

#endif
