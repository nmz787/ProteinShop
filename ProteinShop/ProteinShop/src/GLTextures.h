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
GLTextures - Functions to simplify OpenGL programming.
***********************************************************************/


#ifndef GLTEXTURES_INCLUDED
#define GLTEXTURES_INCLUDED

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

inline void glTexImage2D(GLenum target,GLint level,GLint internalFormat,GLsizei width,GLsizei height,GLint border,GLenum format,const GLubyte* pixels)
{
    glTexImage2D(target,level,internalFormat,width,height,border,format,GL_UNSIGNED_BYTE,pixels);
}

inline void glTexImage2D(GLenum target,GLint level,GLint internalFormat,GLsizei width,GLsizei height,GLint border,GLenum format,const GLbyte* pixels)
{
    glTexImage2D(target,level,internalFormat,width,height,border,format,GL_BYTE,pixels);
}

inline void glTexImage2D(GLenum target,GLint level,GLint internalFormat,GLsizei width,GLsizei height,GLint border,GLenum format,const GLushort* pixels)
{
    glTexImage2D(target,level,internalFormat,width,height,border,format,GL_UNSIGNED_SHORT,pixels);
}

inline void glTexImage2D(GLenum target,GLint level,GLint internalFormat,GLsizei width,GLsizei height,GLint border,GLenum format,const GLshort* pixels)
{
    glTexImage2D(target,level,internalFormat,width,height,border,format,GL_SHORT,pixels);
}

inline void glTexImage2D(GLenum target,GLint level,GLint internalFormat,GLsizei width,GLsizei height,GLint border,GLenum format,const GLuint* pixels)
{
    glTexImage2D(target,level,internalFormat,width,height,border,format,GL_UNSIGNED_INT,pixels);
}

inline void glTexImage2D(GLenum target,GLint level,GLint internalFormat,GLsizei width,GLsizei height,GLint border,GLenum format,const GLint* pixels)
{
    glTexImage2D(target,level,internalFormat,width,height,border,format,GL_INT,pixels);
}

inline void glTexImage2D(GLenum target,GLint level,GLint internalFormat,GLsizei width,GLsizei height,GLint border,GLenum format,const GLfloat* pixels)
{
    glTexImage2D(target,level,internalFormat,width,height,border,format,GL_FLOAT,pixels);
}

template <class TexelType>
void glTexSubImage2D(GLenum target,GLint level,GLint xoffset,GLint yoffset,GLsizei width,GLsizei height,GLint columnStride,GLint rowStride,GLenum format,GLenum type,const TexelType* pixels)
{
    /* Set the common pixel pipeline parameters: */
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS,0);

    if ( columnStride == 1 && rowStride == width )
    {
        /* Upload the texture image directly: */
        glPixelStorei(GL_UNPACK_ROW_LENGTH,rowStride);
        glTexSubImage2D(target,level,xoffset,yoffset,width,height,format,type,pixels);
    }
    else
    {
        /* Copy the texture image to a contiguous buffer: */
        TexelType* tempPixels=new TexelType[width*height];
        TexelType* tempPixelPtr=tempPixels;
        const TexelType* rowPtr=pixels;
        for(int y=0;y<height;++y,rowPtr+=rowStride)
        {
            const TexelType* columnPtr=rowPtr;
            for(int x=0;x<width;++x,++tempPixelPtr,columnPtr+=columnStride)
                *tempPixelPtr=*columnPtr;
        }
        /* Upload the temporary texture image: */
        glPixelStorei(GL_UNPACK_ROW_LENGTH,0);
        glTexSubImage2D(target,level,xoffset,yoffset,width,height,format,type,tempPixels);
        
        /* Delete the temporary texture image: */
        delete[] tempPixels;
    }
}

#endif
