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


#ifndef PALETTERENDERER_INCLUDED
#define PALETTERENDERER_INCLUDED

#ifdef _WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <GLColorMap.h>

#include "VolumeRenderer.h"

class PaletteRenderer:public VolumeRenderer
{
    /* Embedded classes: */
    struct DataItem:public VolumeRenderer::DataItem
    {
        /* Elements: */
        public:
        unsigned int cachedColorMapVersion; // The currently cached color map version

        /* More texture cache update flags: */
        bool uploadColorMap; // Flag to upload color map during texture upload (costly and not always necessary)

        /* Constructors and destructors: */
        DataItem(void);

        /* Methods inherited from VolumeRenderer::ContextData: */
        void updateTextureCache(const VolumeRenderer* renderer,int majorAxis);
        void deleteTextureCache(void);
    };

    static bool init();     // initialize static elements
    static bool isInited;   // flag set by init()

    friend class DataItem;

    /* Static elements: */
    protected:
    #ifdef _WIN32
    static PFNGLCOLORTABLEPROC glColorTable;
    #endif

    /* Elements: */
    unsigned int colorMapVersion; // Version number of the current color map
    const GLColorMap* colorMap; // Colormap containing transfer functions. Must have 256 entries
    bool sharePalette; // Flag if palette renderer uses the global texture palette

    /* Protected methods inherited from VolumeRenderer: */
    void uploadTexture2D(VolumeRenderer::DataItem* dataItem,int axis,int index) const;
    void uploadTexture3D(VolumeRenderer::DataItem* dataItem) const;
    void prepareRenderAxisAligned(VolumeRenderer::DataItem* dataItem) const;

    /* Constructors and destructors: */
    public:
    PaletteRenderer(void) // Creates an uninitialized paletted renderer
        :colorMapVersion(0),colorMap(0),sharePalette(true)
    {
    };
    PaletteRenderer(const char* filename) // Loads a private voxel block from a volume file
        :VolumeRenderer(filename),colorMapVersion(0),colorMap(0),sharePalette(true)
    {
    };
    PaletteRenderer(const Voxel* sValues,const int sSize[3],int sBorderSize,VoxelAlignment sAlignment) // Sets the paletted renderer to a volume block
        :VolumeRenderer(sValues,sSize,sBorderSize,sAlignment),colorMapVersion(0),colorMap(0),sharePalette(true)
    {
    };
    /* Methods inherited from VolumeRenderer: */
    void initContext(GLContextData& contextData) const;
    void setGLState(GLContextData& contextData) const;
    void resetGLState(GLContextData& contextData) const;

    /* New methods: */
    bool isContextInitialized (GLContextData &contextData) const
    {
        return ( contextData.retrieveDataItem<DataItem>(this) );
    }
    const GLColorMap* getColorMap(void) const // Returns a pointer to the used color map
    {
        return colorMap;
    };
    void setColorMap(const GLColorMap* newColorMap) // Sets a new color map
    {
        if(newColorMap->getNumEntries()==256) // Only 8-bit palettes!
        {
            ++colorMapVersion;
            colorMap=newColorMap;
        }
    };
    void setSharePalette(bool newSharePalette) // Sets color palette sharing flag
    {
        sharePalette=newSharePalette;
    };
    static void setGlobalColorMap(const GLColorMap* newGlobalColorMap); // Uploads a global color map to be shared by multiple palette renderers
};

#endif
