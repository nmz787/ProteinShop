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
GLVertexArrays - Helper functions to simplify vertex array rendering.
***********************************************************************/

#ifndef GLVERTEXARRAYS_INCLUDED
#define GLVERTEXARRAYS_INCLUDED

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLTypes.h>

/*************************
Helper class declarations:
*************************/

class GLVertexParts // Class to select vertex information passed to OpenGL
{
    /* Embedded classes: */
    public:
    enum PartsMasks // Masks to select vertex information passed to OpenGL
    {
        Noone=0x0,Color=0x1,Normal=0x2,TexCoord=0x4
    };
    
    /* Methods: */
    static void enableClientState(GLenum vertexPartsMask)
    {
        /* Enable the selected arrays: */
        if(vertexPartsMask&Color)
            glEnableClientState(GL_COLOR_ARRAY);
        if(vertexPartsMask&Normal)
            glEnableClientState(GL_NORMAL_ARRAY);
        if(vertexPartsMask&TexCoord)
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glEnableClientState(GL_VERTEX_ARRAY);
    };
    static void disableClientState(GLenum vertexPartsMask)
    {
        /* Disable the selected arrays: */
        if(vertexPartsMask&Color)
            glDisableClientState(GL_COLOR_ARRAY);
        if(vertexPartsMask&Normal)
            glDisableClientState(GL_NORMAL_ARRAY);
        if(vertexPartsMask&TexCoord)
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    };
};


/**************************************
Specialized versions of class GLVertex:
**************************************/

template <class TexCoordScalarType =GLfloat,int texCoordComponents =2,
          class ColorScalarType =GLubyte,int colorComponents =4,
          class NormalScalarType =GLfloat,
          class VertexScalarType =GLfloat,int vertexComponents =3>
class GLVertex // Class for fully-specified GL vertices
{
    /* Embedded classes: */
    public:
    typedef GLVector<TexCoordScalarType,texCoordComponents> TexCoord;
    typedef GLColor<ColorScalarType,colorComponents> Color;
    typedef GLVector<NormalScalarType,3> Normal;
    typedef GLVector<VertexScalarType,vertexComponents> Vertex;
    
    /* Elements: */
    TexCoord texCoord; // Vertex texture coordinates
    Color color; // Vertex color
    Normal normal; // Vertex normal
    Vertex vertex; // Vertex position
    
    /* Constructors and destructors: */
    GLVertex(void)
    {
    };
    GLVertex(const TexCoord& sTexCoord,const Color& sColor,const Normal& sNormal,const Vertex& sVertex)
        :texCoord(sTexCoord),color(sColor),normal(sNormal),vertex(sVertex)
    {
    };
};

template <class TexCoordScalarType,int texCoordComponents,
          class NormalScalarType,
          class VertexScalarType,int vertexComponents>
class GLVertex<TexCoordScalarType,texCoordComponents,void,0,NormalScalarType,VertexScalarType,vertexComponents>
{
    /* Embedded classes: */
    public:
    typedef GLVector<TexCoordScalarType,texCoordComponents> TexCoord;
    typedef GLVector<NormalScalarType,3> Normal;
    typedef GLVector<VertexScalarType,vertexComponents> Vertex;
    
    /* Elements: */
    TexCoord texCoord; // Vertex texture coordinates
    Normal normal; // Vertex normal
    Vertex vertex; // Vertex position
    
    /* Constructors and destructors: */
    GLVertex(void)
    {
    };
    GLVertex(const TexCoord& sTexCoord,const Normal& sNormal,const Vertex& sVertex)
        :texCoord(sTexCoord),normal(sNormal),vertex(sVertex)
    {
    };
};

template <class TexCoordScalarType,int texCoordComponents,
          class ColorScalarType,int colorComponents,
          class VertexScalarType,int vertexComponents>
class GLVertex<TexCoordScalarType,texCoordComponents,ColorScalarType,colorComponents,void,VertexScalarType,vertexComponents>
{
    /* Embedded classes: */
    public:
    typedef GLVector<TexCoordScalarType,texCoordComponents> TexCoord;
    typedef GLColor<ColorScalarType,colorComponents> Color;
    typedef GLVector<VertexScalarType,vertexComponents> Vertex;
    
    /* Elements: */
    TexCoord texCoord; // Vertex texture coordinates
    Color color; // Vertex color
    Vertex vertex; // Vertex position
    
    /* Constructors and destructors: */
    GLVertex(void)
    {
    };
    GLVertex(const TexCoord& sTexCoord,const Color& sColor,const Vertex& sVertex)
        :texCoord(sTexCoord),color(sColor),vertex(sVertex)
    {
    };
};

template <class ColorScalarType,int colorComponents,
          class NormalScalarType,
          class VertexScalarType,int vertexComponents>
class GLVertex<void,0,ColorScalarType,colorComponents,NormalScalarType,VertexScalarType,vertexComponents>
{
    /* Embedded classes: */
    public:
    typedef GLColor<ColorScalarType,colorComponents> Color;
    typedef GLVector<NormalScalarType,3> Normal;
    typedef GLVector<VertexScalarType,vertexComponents> Vertex;
    
    /* Elements: */
    Color color; // Vertex color
    Normal normal; // Vertex normal
    Vertex vertex; // Vertex position
    
    /* Constructors and destructors: */
    GLVertex(void)
    {
    };
    GLVertex(const Color& sColor,const Normal& sNormal,const Vertex& sVertex)
        :color(sColor),normal(sNormal),vertex(sVertex)
    {
    };
};

template <class TexCoordScalarType,int texCoordComponents,
          class VertexScalarType,int vertexComponents>
class GLVertex<TexCoordScalarType,texCoordComponents,void,0,void,VertexScalarType,vertexComponents>
{
    /* Embedded classes: */
    public:
    typedef GLVector<TexCoordScalarType,texCoordComponents> TexCoord;
    typedef GLVector<VertexScalarType,vertexComponents> Vertex;
    
    /* Elements: */
    TexCoord texCoord; // Vertex texture coordinates
    Vertex vertex; // Vertex position
    
    /* Constructors and destructors: */
    GLVertex(void)
    {
    };
    GLVertex(const TexCoord& sTexCoord,const Vertex& sVertex)
        :texCoord(sTexCoord),vertex(sVertex)
    {
    };
};

template <class NormalScalarType,
          class VertexScalarType,int vertexComponents>
class GLVertex<void,0,void,0,NormalScalarType,VertexScalarType,vertexComponents>
{
    /* Embedded classes: */
    public:
    typedef GLVector<NormalScalarType,3> Normal;
    typedef GLVector<VertexScalarType,vertexComponents> Vertex;
    
    /* Elements: */
    Normal normal; // Vertex normal
    Vertex vertex; // Vertex position
    
    /* Constructors and destructors: */
    GLVertex(void)
    {
    };
    GLVertex(const Normal& sNormal,const Vertex& sVertex)
        :normal(sNormal),vertex(sVertex)
    {
    };
};

template <class ColorScalarType,int colorComponents,
          class VertexScalarType,int vertexComponents>
class GLVertex<void,0,ColorScalarType,colorComponents,void,VertexScalarType,vertexComponents>
{
    /* Embedded classes: */
    public:
    typedef GLColor<ColorScalarType,colorComponents> Color;
    typedef GLVector<VertexScalarType,vertexComponents> Vertex;
    
    /* Elements: */
    Color color; // Vertex color
    Vertex vertex; // Vertex position
    
    /* Constructors and destructors: */
    GLVertex(void)
    {
    };
    GLVertex(const Color& sColor,const Vertex& sVertex)
        :color(sColor),vertex(sVertex)
    {
    };
};

template <class VertexScalarType,int vertexComponents>
class GLVertex<void,0,void,0,void,VertexScalarType,vertexComponents>
{
    /* Embedded classes: */
    public:
    typedef GLVector<VertexScalarType,vertexComponents> Vertex;
    
    /* Elements: */
    Vertex vertex; // Vertex position
    
    /* Constructors and destructors: */
    GLVertex(void)
    {
    };
    GLVertex(const Vertex& sVertex)
        :vertex(sVertex)
    {
    };
};

/**********************************
Specialized versions of glVertex():
**********************************/

template <class TexCoordScalarType,int texCoordComponents,
          class ColorScalarType,int colorComponents,
          class NormalScalarType,
          class VertexScalarType,int vertexComponents>
inline void glVertex(const GLVertex<TexCoordScalarType,texCoordComponents,ColorScalarType,colorComponents,NormalScalarType,VertexScalarType,vertexComponents>& v)
{
    glTexCoord(v.texCoord);
    glColor(v.color);
    glNormal(v.normal);
    glVertex(v.vertex);
}

template <class TexCoordScalarType,int texCoordComponents,
          class NormalScalarType,
          class VertexScalarType,int vertexComponents>
inline void glVertex(const GLVertex<TexCoordScalarType,texCoordComponents,void,0,NormalScalarType,VertexScalarType,vertexComponents>& v)
{
    glTexCoord(v.texCoord);
    glNormal(v.normal);
    glVertex(v.vertex);
}

template <class TexCoordScalarType,int texCoordComponents,
          class ColorScalarType,int colorComponents,
          class VertexScalarType,int vertexComponents>
inline void glVertex(const GLVertex<TexCoordScalarType,texCoordComponents,ColorScalarType,colorComponents,void,VertexScalarType,vertexComponents>& v)
{
    glTexCoord(v.texCoord);
    glColor(v.color);
    glVertex(v.vertex);
}

template <class ColorScalarType,int colorComponents,
          class NormalScalarType,
          class VertexScalarType,int vertexComponents>
inline void glVertex(const GLVertex<void,0,ColorScalarType,colorComponents,NormalScalarType,VertexScalarType,vertexComponents>& v)
{
    glColor(v.color);
    glNormal(v.normal);
    glVertex(v.vertex);
}

template <class TexCoordScalarType,int texCoordComponents,
          class VertexScalarType,int vertexComponents>
inline void glVertex(const GLVertex<TexCoordScalarType,texCoordComponents,void,0,void,VertexScalarType,vertexComponents>& v)
{
    glTexCoord(v.texCoord);
    glVertex(v.vertex);
}

template <class NormalScalarType,
          class VertexScalarType,int vertexComponents>
inline void glVertex(const GLVertex<void,0,void,0,NormalScalarType,VertexScalarType,vertexComponents>& v)
{
    glNormal(v.normal);
    glVertex(v.vertex);
}

template <class ColorScalarType,int colorComponents,
          class VertexScalarType,int vertexComponents>
inline void glVertex(const GLVertex<void,0,ColorScalarType,colorComponents,void,VertexScalarType,vertexComponents>& v)
{
    glColor(v.color);
    glVertex(v.vertex);
}

template <class VertexScalarType,int vertexComponents>
inline void glVertex(const GLVertex<void,0,void,0,void,VertexScalarType,vertexComponents>& v)
{
    glVertex(v.vertex);
}

/*****************************************
Specialized versions of glVertexPointer():
*****************************************/

template <class TexCoordScalarType,int texCoordComponents,
          class ColorScalarType,int colorComponents,
          class NormalScalarType,
          class VertexScalarType,int vertexComponents>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<TexCoordScalarType,texCoordComponents,ColorScalarType,colorComponents,NormalScalarType,VertexScalarType,vertexComponents>* vertexArray)
{
    /* Set the selected arrays: */
    if(vertexPartsMask&GLVertexParts::TexCoord)
        glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
    if(vertexPartsMask&GLVertexParts::Color)
        glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
    if(vertexPartsMask&GLVertexParts::Normal)
        glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
    glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
}

template <class TexCoordScalarType,int texCoordComponents,
          class NormalScalarType,
          class VertexScalarType,int vertexComponents>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<TexCoordScalarType,texCoordComponents,void,0,NormalScalarType,VertexScalarType,vertexComponents>* vertexArray)
{
    /* Set the selected arrays: */
    if(vertexPartsMask&GLVertexParts::TexCoord)
        glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
    if(vertexPartsMask&GLVertexParts::Normal)
        glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
    glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
}

template <class TexCoordScalarType,int texCoordComponents,
          class ColorScalarType,int colorComponents,
          class VertexScalarType,int vertexComponents>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<TexCoordScalarType,texCoordComponents,ColorScalarType,colorComponents,void,VertexScalarType,vertexComponents>* vertexArray)
{
    /* Set the selected arrays: */
    if(vertexPartsMask&GLVertexParts::TexCoord)
        glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
    if(vertexPartsMask&GLVertexParts::Color)
        glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
    glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
}

template <class ColorScalarType,int colorComponents,
          class NormalScalarType,
          class VertexScalarType,int vertexComponents>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,ColorScalarType,colorComponents,NormalScalarType,VertexScalarType,vertexComponents>* vertexArray)
{
    /* Set the selected arrays: */
    if(vertexPartsMask&GLVertexParts::Color)
        glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
    if(vertexPartsMask&GLVertexParts::Normal)
        glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
    glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
}

template <class TexCoordScalarType,int texCoordComponents,
          class VertexScalarType,int vertexComponents>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<TexCoordScalarType,texCoordComponents,void,0,void,VertexScalarType,vertexComponents>* vertexArray)
{
    /* Set the selected arrays: */
    if(vertexPartsMask&GLVertexParts::TexCoord)
        glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
    glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
}

template <class NormalScalarType,
          class VertexScalarType,int vertexComponents>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,void,0,NormalScalarType,VertexScalarType,vertexComponents>* vertexArray)
{
    /* Set the selected arrays: */
    if(vertexPartsMask&GLVertexParts::Normal)
        glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
    glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
}

template <class ColorScalarType,int colorComponents,
          class VertexScalarType,int vertexComponents>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,ColorScalarType,colorComponents,void,VertexScalarType,vertexComponents>* vertexArray)
{
    /* Set the selected arrays: */
    if(vertexPartsMask&GLVertexParts::Color)
        glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
    glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
}

template <class VertexScalarType,int vertexComponents>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,void,0,void,VertexScalarType,vertexComponents>* vertexArray)
{
    /* Set the selected arrays: */
    glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<GLfloat,4,GLfloat,4,GLfloat,GLfloat,4>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::TexCoord|GLVertexParts::Color|GLVertexParts::Normal)
        glInterleavedArrays(GL_T4F_C4F_N3F_V4F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::TexCoord)
            glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
        if(vertexPartsMask&GLVertexParts::Color)
            glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
        if(vertexPartsMask&GLVertexParts::Normal)
            glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<GLfloat,2,GLfloat,4,GLfloat,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::TexCoord|GLVertexParts::Color|GLVertexParts::Normal)
        glInterleavedArrays(GL_T2F_C4F_N3F_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::TexCoord)
            glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
        if(vertexPartsMask&GLVertexParts::Color)
            glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
        if(vertexPartsMask&GLVertexParts::Normal)
            glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<GLfloat,2,void,0,GLfloat,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::TexCoord|GLVertexParts::Normal)
        glInterleavedArrays(GL_T2F_N3F_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::TexCoord)
            glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
        if(vertexPartsMask&GLVertexParts::Normal)
            glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<GLfloat,2,GLfloat,3,void,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::TexCoord|GLVertexParts::Color)
        glInterleavedArrays(GL_T2F_C3F_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::TexCoord)
            glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
        if(vertexPartsMask&GLVertexParts::Color)
            glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<GLfloat,2,GLubyte,4,void,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::TexCoord|GLVertexParts::Color)
        glInterleavedArrays(GL_T2F_C4UB_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::TexCoord)
            glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
        if(vertexPartsMask&GLVertexParts::Color)
            glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<GLfloat,4,void,0,void,GLfloat,4>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::TexCoord)
        glInterleavedArrays(GL_T4F_V4F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::TexCoord)
            glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<GLfloat,2,void,0,void,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::TexCoord)
        glInterleavedArrays(GL_T2F_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::TexCoord)
            glTexCoordPointer(sizeof(vertexArray[0]),&vertexArray[0].texCoord);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,GLfloat,4,GLfloat,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::Color|GLVertexParts::Normal)
        glInterleavedArrays(GL_C4F_N3F_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::Color)
            glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
        if(vertexPartsMask&GLVertexParts::Normal)
            glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,void,0,GLfloat,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::Normal)
        glInterleavedArrays(GL_N3F_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::Normal)
            glNormalPointer(sizeof(vertexArray[0]),&vertexArray[0].normal);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,GLfloat,3,void,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::Color)
        glInterleavedArrays(GL_C3F_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::Color)
            glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,GLubyte,4,void,GLfloat,3>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::Color)
        glInterleavedArrays(GL_C4UB_V3F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::Color)
            glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,GLubyte,4,void,GLfloat,2>* vertexArray)
{
    if(vertexPartsMask==GLVertexParts::Color)
        glInterleavedArrays(GL_C4UB_V2F,0,vertexArray);
    else
    {
        /* Set the selected arrays: */
        if(vertexPartsMask&GLVertexParts::Color)
            glColorPointer(sizeof(vertexArray[0]),&vertexArray[0].color);
        glVertexPointer(sizeof(vertexArray[0]),&vertexArray[0].vertex);
    }
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,void,0,void,GLfloat,3>* vertexArray)
{
    glInterleavedArrays(GL_V3F,0,vertexArray);
}

template <>
inline void glVertexPointer(GLenum vertexPartsMask,const GLVertex<void,0,void,0,void,GLfloat,2>* vertexArray)
{
    glInterleavedArrays(GL_V2F,0,vertexArray);
}

#endif
