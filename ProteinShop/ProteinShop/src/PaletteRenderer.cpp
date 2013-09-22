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
PaletteRenderer - Class for texture-based volume renderers using palette-
based transfer functions.
***********************************************************************/
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#include <GL/glext.h>
#endif

#include <GLTextures.h>

#include "PaletteRenderer.h"


bool PaletteRenderer::isInited = PaletteRenderer::init();


/******************************************
Methods of class PaletteRenderer::DataItem:
******************************************/

PaletteRenderer::DataItem::DataItem(void)
    :cachedColorMapVersion(0),uploadColorMap(false)
{
}

void PaletteRenderer::DataItem::updateTextureCache(const VolumeRenderer* renderer,int majorAxis)
{
    uploadColorMap=false;

    /* Call the base class method: */
    VolumeRenderer::DataItem::updateTextureCache(renderer,majorAxis);

    /* Check if the color map needs to be uploaded: */
    const PaletteRenderer* paletteRenderer=static_cast<const PaletteRenderer*>(renderer);

    /*************************************************************************
    Another bug in SGI's OpenGL implementation: Though color maps are treated
    as a texture object resource, they are not always installed when a texture
    object is bound. Thus, we have to upload them manually everytime we bind;
    painful and sloow.
    *************************************************************************/

    #ifndef __SGI_IRIX__
    if(uploadData||cachedColorMapVersion!=paletteRenderer->colorMapVersion)
    #endif
    {
        cachedColorMapVersion=paletteRenderer->colorMapVersion;
        textureCacheValid=false;
        uploadColorMap=!paletteRenderer->sharePalette;
    }
}

void PaletteRenderer::DataItem::deleteTextureCache(void)
{
    /* Call the base class method: */
    VolumeRenderer::DataItem::deleteTextureCache();

    uploadColorMap=false;
}

/****************************************
Static elements of class PaletteRenderer:
****************************************/

#ifdef _WIN32
PFNGLCOLORTABLEPROC PaletteRenderer::glColorTable;
#endif

/********************************
Methods of class PaletteRenderer:
********************************/

void PaletteRenderer::uploadTexture2D(VolumeRenderer::DataItem* dataItem,int axis,int index) const
{
    if(dataItem->setParameters)
    {
        /* Set the OpenGL texturing parameters: */
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_BASE_LEVEL,0);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAX_LEVEL,0);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,interpolationMode);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,interpolationMode);
    }

    #if 1
    /* Upload a color map only if necessary: */
    if(static_cast<DataItem*>(dataItem)->uploadColorMap)
    {
        /* Set the texture's color map: */
        #ifdef __SGI_IRIX__
        glColorTableSGI(GL_TEXTURE_COLOR_TABLE_SGI,GL_RGBA,256,GL_RGBA,GL_FLOAT,colorMap->getColors());
        #else
        glColorTable(GL_TEXTURE_2D,GL_RGBA,256,GL_RGBA,GL_FLOAT,colorMap->getColors());
        #endif
    }
    #endif

    if(dataItem->uploadData)
    {
        /* Some workaround for SGI's non-compliance of OpenGL (what a shame): */
        #if (__SGI_IRIX__ || __APPLE__)
        const GLenum internalFormat=GL_INTENSITY8;
        const GLenum uploadFormat=GL_LUMINANCE;
        #else
        // use non-ARB extension 78
        const GLenum internalFormat=GL_COLOR_INDEX8_EXT;
        const GLenum uploadFormat=GL_COLOR_INDEX;
        #endif

        /* Upload a texture slice: */
        const Voxel* slicePtr=values+index*increments[axis];
        switch(axis)
        {
            case 0:
                glTexImage2D(GL_TEXTURE_2D,0,internalFormat,textureSize[2],textureSize[1],0,uploadFormat,GL_UNSIGNED_BYTE,0);
                glTexSubImage2D(GL_TEXTURE_2D,0,0,0,size[2],size[1],increments[2],increments[1],uploadFormat,GL_UNSIGNED_BYTE,slicePtr);
                break;
            case 1:
                glTexImage2D(GL_TEXTURE_2D,0,internalFormat,textureSize[2],textureSize[0],0,uploadFormat,GL_UNSIGNED_BYTE,0);
                glTexSubImage2D(GL_TEXTURE_2D,0,0,0,size[2],size[0],increments[2],increments[0],uploadFormat,GL_UNSIGNED_BYTE,slicePtr);
                break;
            case 2:
                glTexImage2D(GL_TEXTURE_2D,0,internalFormat,textureSize[1],textureSize[0],0,uploadFormat,GL_UNSIGNED_BYTE,0);
                glTexSubImage2D(GL_TEXTURE_2D,0,0,0,size[1],size[0],increments[1],increments[0],uploadFormat,GL_UNSIGNED_BYTE,slicePtr);
                break;
        }
    }
}

void PaletteRenderer::uploadTexture3D(VolumeRenderer::DataItem* dataItem) const
{
    if(dataItem->setParameters)
    {
        /* Set the OpenGL texturing parameters: */
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_BASE_LEVEL,0);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAX_LEVEL,0);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,interpolationMode);
        glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,interpolationMode);
    }

    #if 1
    /* Upload a color map only if necessary: */
    if(static_cast<DataItem*>(dataItem)->uploadColorMap)
    {
        /* Set the texture's color map: */
        #ifdef __SGI_IRIX__
        glColorTableSGI(GL_TEXTURE_COLOR_TABLE_SGI,GL_RGBA,256,GL_RGBA,GL_FLOAT,colorMap->getColors());
        #else
        glColorTable(GL_TEXTURE_3D,GL_RGBA,256,GL_RGBA,GL_FLOAT,colorMap->getColors());
        #endif
    }
    #endif
    
    if(dataItem->uploadData)
    {
        /* Some workaround for SGI's non-compliance of OpenGL (what a shame): */
        #if (__SGI_IRIX__ || __APPLE__)
        const GLenum internalFormat=GL_INTENSITY8;
        const GLenum uploadFormat=GL_LUMINANCE;
        #else
        // use non-ARB extension 78
        const GLenum internalFormat=GL_COLOR_INDEX8_EXT;
        const GLenum uploadFormat=GL_COLOR_INDEX;
        #endif

        /* Upload the texture block: */
        glPixelStorei(GL_UNPACK_ALIGNMENT,1);
        glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
        glPixelStorei(GL_UNPACK_ROW_LENGTH,0); // increments[1]); // Seems to be a bug in OpenGL - consistent across SGI/nVidia platforms
        glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
        glPixelStorei(GL_UNPACK_IMAGE_HEIGHT,0); // increments[0]);
        glPixelStorei(GL_UNPACK_SKIP_IMAGES,0);
        #ifdef __SGI_IRIX__
        glTexImage3D(GL_TEXTURE_3D,0,internalFormat,textureSize[2],textureSize[1],textureSize[0],0,uploadFormat,GL_UNSIGNED_BYTE,values);
        #else
        glTexImage3D(GL_TEXTURE_3D,0,internalFormat,textureSize[2],textureSize[1],textureSize[0],0,uploadFormat,GL_UNSIGNED_BYTE,0);
        glTexSubImage3D(GL_TEXTURE_3D,0,0,0,0,size[2],size[1],size[0],uploadFormat,GL_UNSIGNED_BYTE,values);
        #endif
    }
}

void PaletteRenderer::prepareRenderAxisAligned(VolumeRenderer::DataItem* dataItem) const
{
    //DataItem* myDataItem=static_cast<DataItem*>(dataItem);

    #if 1
    /* Manage color palette uploads: */
    if(!textureCachingEnabled)
    {
        /* Sufficient to upload palette right here: */
        #ifdef __SGI_IRIX__
        glColorTableSGI(GL_TEXTURE_COLOR_TABLE_SGI,GL_RGBA,256,GL_RGBA,GL_FLOAT,colorMap->getColors());
        #else
        glColorTable(GL_TEXTURE_2D,GL_RGBA,256,GL_RGBA,GL_FLOAT,colorMap->getColors());
        #endif
    }
    #endif
}

bool PaletteRenderer::init(void)
{
    #ifdef _WIN32
    /* Load OpenGL function pointers from DLL: */
    glColorTable=(PFNGLCOLORTABLEPROC)wglGetProcAddress("glColorTable");
    #endif

    /* Check for applicable OpenGL extensions: */
    /* ... */

    return true;
}

void PaletteRenderer::initContext(GLContextData& contextData) const
{
    contextData.addDataItem(this,new DataItem);
}

void PaletteRenderer::setGLState(GLContextData& contextData) const
{
    VolumeRenderer::setGLState(contextData);

    #ifdef __SGI_IRIX__
    glEnable(GL_TEXTURE_COLOR_TABLE_SGI);
    #endif
    if(sharePalette)
        glEnable(GL_SHARED_TEXTURE_PALETTE_EXT);

    #if 0
    /* Sufficient to upload palette right here: */
    #ifdef __SGI_IRIX__
    glColorTableSGI(GL_TEXTURE_COLOR_TABLE_SGI,GL_RGBA,256,GL_RGBA,GL_FLOAT,colorMap->getColors());
    #else
    glColorTable(GL_TEXTURE_2D,GL_RGBA,256,GL_RGBA,GL_FLOAT,colorMap->getColors());
    #endif
    #endif
}

void PaletteRenderer::resetGLState(GLContextData& contextData) const
{
    #ifdef __SGI_IRIX__
    glDisable(GL_TEXTURE_COLOR_TABLE_SGI);
    #endif
    if(sharePalette)
        glDisable(GL_SHARED_TEXTURE_PALETTE_EXT);

    VolumeRenderer::resetGLState(contextData);
}

void PaletteRenderer::setGlobalColorMap(const GLColorMap* newGlobalColorMap)
{
    if(newGlobalColorMap->getNumEntries()==256)
    {
        /* Sufficient to upload palette right here: */
        #ifdef __SGI_IRIX__
        glColorTableSGI(GL_TEXTURE_COLOR_TABLE_SGI,GL_RGBA,256,GL_RGBA,GL_FLOAT,newGlobalColorMap->getColors());
        #else
        glColorTable(GL_SHARED_TEXTURE_PALETTE_EXT,GL_RGBA,256,GL_RGBA,GL_FLOAT,newGlobalColorMap->getColors());
        #endif
    }
}
