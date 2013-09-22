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
VolumeRenderer - Base class for texture-based volume renderers for blocks
of cartesian voxel data.

***********************************************************************/


#include <cassert>
#include <cstdio>
#include <iostream>
using namespace std;

#include <Geometry/HVector.h>
#include <Geometry/ProjectiveTransformation.h>

#include "Buffer.h"
#include "GLTextures.h"
#include "PriorityHeap.h"
#include "VolumeRenderer.h"


bool VolumeRenderer::isInited = VolumeRenderer::init();


static const GLint f_clearBufferTexelFormat = GL_ALPHA;
static const GLenum f_clearBufferPixelType = GL_UNSIGNED_BYTE;

/** A static buffer for the "clearing texture" used in the initial
    glTexImage*D() calls. */
const GLvoid *clearBuffer (GLsizei width, GLsizei height, GLsizei depth = 1)
{
    static Buffer<GLubyte,3> buffer;

    assert ( width && height && depth );
    uint previousCapacity = buffer.capacity();
    buffer.setDimSizes (width, height, depth);
    uint newCapacity = buffer.capacity();
    if ( previousCapacity != newCapacity )
        buffer.initElements (0);
    return (const GLvoid *) buffer;
}


/*****************************************
Methods of class VolumeRenderer::DataItem:
*****************************************/

VolumeRenderer::DataItem::DataItem(void) :
    dataVersion (0),
    settingsVersion (0),
    numTextureObjects (0),
    textureObjectIDs (0),
    textureCacheValid (false),
    cachedAxis (0),
    setParameters (true),
    uploadData (true)
{
}

VolumeRenderer::DataItem::~DataItem(void)
{
    deleteTextureCache();
}

void VolumeRenderer::DataItem::updateTextureCache(const VolumeRenderer* renderer,int majorAxis)
{
    setParameters=false;
    uploadData=false;
    if(dataVersion!=renderer->dataVersion||majorAxis!=cachedAxis)
    {
        /* Calculate the number of required texture objects: */
        int requiredNumTextures=1;
        if(renderer->renderingMode==VolumeRenderer::AXIS_ALIGNED)
        {
            for(int i=0;i<3;++i)
                if(requiredNumTextures<renderer->size[i])
                    requiredNumTextures=renderer->size[i];
        }
        
        /* Reallocate the texture cache if necessary: */
        if(numTextureObjects!=requiredNumTextures)
        {
            if(numTextureObjects!=0)
            {
                glDeleteTextures(numTextureObjects,textureObjectIDs);
                delete[] textureObjectIDs;
            }
            numTextureObjects=requiredNumTextures;
            textureObjectIDs=new GLuint[numTextureObjects];
            glGenTextures(numTextureObjects,textureObjectIDs);
        }
        
        /* Invalidate the texture cache: */
        dataVersion=renderer->dataVersion;
        cachedAxis=majorAxis;
        textureCacheValid=false;
        setParameters=true;
        uploadData=true;
    }
    
    if(settingsVersion!=renderer->settingsVersion)
    {
        /* Invalidate the texture cache: */
        settingsVersion=renderer->settingsVersion;
        textureCacheValid=false;
        setParameters=true;
    }
}

void VolumeRenderer::DataItem::deleteTextureCache(void)
{
    if(numTextureObjects!=0)
    {
        /* Delete all textures: */
        glDeleteTextures(numTextureObjects,textureObjectIDs);

        /* Delete the texture cache: */
        numTextureObjects=0;
        delete[] textureObjectIDs;
        textureObjectIDs=0;
        textureCacheValid=false;
        cachedAxis=-1;
    }
    setParameters=true;
    uploadData=true;
}

/***************************************
Static elements of class VolumeRenderer:
***************************************/

#ifdef _WIN32
PFNGLTEXIMAGE3DPROC VolumeRenderer::glTexImage3D=0;
PFNGLTEXSUBIMAGE3DPROC VolumeRenderer::glTexSubImage3D=0;
#endif

/*******************************
Methods of class VolumeRenderer:
*******************************/

void VolumeRenderer::deletePrivateData(void)
{
    if(privateData)
    {
        /* Calculate the base address of the voxel block: */
        const Voxel* valueBase=values;
        for(int i=0;i<3;++i)
            valueBase-=borderSize*increments[i];
        
        /* Delete the allocated voxel block: */
        delete[] valueBase;
    }
    
    /* Reset elements: */
    privateData=false;
    values=0;
}

void VolumeRenderer::uploadTexture2D(VolumeRenderer::DataItem* dataItem,int axis,int index) const
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

    if(dataItem->uploadData)
    {
        /* Upload a texture slice: */
        const Voxel* slicePtr=values+index*increments[axis];
        switch(axis)
        {
            case 0:
                glTexImage2D (
                    GL_TEXTURE_2D,
                    0,
                    GL_RGBA,
                    textureSize[2],
                    textureSize[1],
                    0,
                    f_clearBufferTexelFormat,
                    f_clearBufferPixelType,
                    clearBuffer(textureSize[2], textureSize[1])
                );
                // note that this call is in GLTextures.h, not OpenGL
                glTexSubImage2D (
                    GL_TEXTURE_2D,
                    0,
                    0,
                    0,
                    size[2],
                    size[1],
                    increments[2],
                    increments[1],
                    GL_RGBA,
                    GL_UNSIGNED_INT_8_8_8_8,
                    slicePtr
                );
                break;
            case 1:
                glTexImage2D (
                    GL_TEXTURE_2D,
                    0,
                    GL_RGBA,
                    textureSize[2],
                    textureSize[0],
                    0,
                    f_clearBufferTexelFormat,
                    f_clearBufferPixelType,
                    clearBuffer(textureSize[2], textureSize[0])
                );
                // note that this call is in GLTextures.h, not OpenGL
                glTexSubImage2D (
                    GL_TEXTURE_2D,
                    0,
                    0,
                    0,
                    size[2],
                    size[0],
                    increments[2],
                    increments[0],
                    GL_RGBA,
                    GL_UNSIGNED_INT_8_8_8_8,
                    slicePtr
                );
                break;
            case 2:
                glTexImage2D (
                    GL_TEXTURE_2D,
                    0,
                    GL_RGBA,
                    textureSize[1],
                    textureSize[0],
                    0,
                    f_clearBufferTexelFormat,
                    f_clearBufferPixelType,
                    clearBuffer(textureSize[1], textureSize[0])
                );
                // note that this call is in GLTextures.h, not OpenGL
                glTexSubImage2D (
                    GL_TEXTURE_2D,
                    0,
                    0,
                    0,
                    size[1],
                    size[0],
                    increments[1],
                    increments[0],
                    GL_RGBA,
                    GL_UNSIGNED_INT_8_8_8_8,
                    slicePtr
                );
                break;
        }
    }
}

void VolumeRenderer::uploadTexture3D(VolumeRenderer::DataItem* dataItem) const
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
    
    if(dataItem->uploadData)
    {
        /* Upload the texture block: */
        glPixelStorei(GL_UNPACK_ALIGNMENT,1);
        glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
        glPixelStorei(GL_UNPACK_ROW_LENGTH,0); // increments[1]); // Seems to be a bug in OpenGL - consistent across SGI/nVidia platforms
        glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
        glPixelStorei(GL_UNPACK_IMAGE_HEIGHT,0); // increments[0]);
        glPixelStorei(GL_UNPACK_SKIP_IMAGES,0);
        glTexImage3D (
            GL_TEXTURE_3D,
            0,
            GL_RGBA,
            textureSize[2],
            textureSize[1],
            textureSize[0],
            0,
            f_clearBufferTexelFormat,
            f_clearBufferPixelType,
            clearBuffer(textureSize[2], textureSize[1], textureSize[0])
        );
        glTexSubImage3D (
            GL_TEXTURE_3D,
            0,
            0,
            0,
            0,
            size[2],
            size[1],
            size[0],
            GL_RGBA,
            GL_UNSIGNED_INT_8_8_8_8,
            values
        );
    }
}

void VolumeRenderer::prepareRenderAxisAligned(VolumeRenderer::DataItem* dataItem) const
{
}

void VolumeRenderer::renderAxisAligned(VolumeRenderer::DataItem* dataItem,const Vector& viewDirection) const
{
    /* Identify the major rendering axis and the stacking direction: */
    int majorAxis;
    int textureAxis[2];
    if(fabs(viewDirection[0])>=fabs(viewDirection[1])&&fabs(viewDirection[0])>=fabs(viewDirection[2]))
    {
        /* Major axis is x: */
        majorAxis=0;
        textureAxis[0]=2;
        textureAxis[1]=1;
    }
    else if(fabs(viewDirection[1])>=fabs(viewDirection[2]))
    {
        /* Major axis is y: */
        majorAxis=1;
        textureAxis[0]=2;
        textureAxis[1]=0;
    }
    else
    {
        /* Major axis is z: */
        majorAxis=2;
        textureAxis[0]=1;
        textureAxis[1]=0;
    }

    /* Determine stacking order and calculate the slices' corner positions and texture coordinates: */
    static int cornerIndices[3][2][4] = {
        {{0,2,6,4},{1,5,7,3}},
        {{0,4,5,1},{2,3,7,6}},
        {{0,1,3,2},{4,6,7,5}}
    };
    int stackingOrder;
    int sliceIndex,sliceIncrement,lastSlice;
    double quadCornerIncrement;
    Vertex quadCorner[4];
    double quadTexCoord[4][2];
    if(viewDirection[majorAxis]<0.0)
    {
        /* Stacking order is upwards: */
        stackingOrder=0;
        sliceIndex=subOrigin[majorAxis];
        lastSlice=subOrigin[majorAxis]+subSize[majorAxis];
        if(alignment==VERTEX_CENTERED)
            ++lastSlice;
        sliceIncrement=1;
        quadCornerIncrement=corners[7].position[majorAxis]-corners[0].position[majorAxis];
    }
    else
    {
        /* Stacking order is downwards: */
        stackingOrder=1;
        sliceIndex=subOrigin[majorAxis]+subSize[majorAxis]-1;
        if(alignment==VERTEX_CENTERED)
            ++sliceIndex;
        lastSlice=subOrigin[majorAxis]-1;
        sliceIncrement=-1;
        quadCornerIncrement=corners[0].position[majorAxis]-corners[7].position[majorAxis];
    }
    quadCornerIncrement/=double(subSize[majorAxis]);

    /* Copy positions and texture coordinates from the box structure: */
    for(int i=0;i<4;++i)
    {
        const BoxCorner& c=corners[cornerIndices[majorAxis][stackingOrder][i]];
        quadCorner[i]=c.position;
        for(int j=0;j<2;++j)
            quadTexCoord[i][j]=c.texture[2-textureAxis[j]];
    }

    /* Adjust for cell-centered voxels (texture slices are aligned with cell centers): */
    if(alignment==CELL_CENTERED)
        for(int i=0;i<4;++i)
            quadCorner[i][majorAxis]+=quadCornerIncrement*0.5;

    /* Create/delete the texture cache if necessary: */
    if(textureCachingEnabled)
        dataItem->updateTextureCache(this,majorAxis);
    else
        dataItem->deleteTextureCache();

    /* Prepare rendering: */
    prepareRenderAxisAligned(dataItem);

    /* Render each slice as a textured quadrilateral: */
    for(;sliceIndex!=lastSlice;sliceIndex+=sliceIncrement)
    {
        #if 1
        /* Upload the slice texture: */
        if(textureCachingEnabled)
        {
            glBindTexture(GL_TEXTURE_2D,dataItem->textureObjectIDs[sliceIndex]);
            if(!dataItem->textureCacheValid)
                uploadTexture2D(dataItem,majorAxis,sliceIndex);
        }
        else
            uploadTexture2D(dataItem,majorAxis,sliceIndex);
        #endif

        /* Render a quadrilateral: */
        glBegin(GL_QUADS);
        for(int i=0;i<4;++i)
        {
            glTexCoord2dv(quadTexCoord[i]);
            glVertex3dv(quadCorner[i].getComponents());
            quadCorner[i][majorAxis]+=quadCornerIncrement;
        }
        glEnd();
    }
    
    #if 1
    if(textureCachingEnabled)
    {
        /* Unbind the last texture to prevent someone else from tampering with it: */
        glBindTexture(GL_TEXTURE_2D,0);
        
        /* Validate the texture cache: */
        dataItem->textureCacheValid=true;
    }
    #endif
}

void VolumeRenderer::prepareRenderViewPerpendicular(VolumeRenderer::DataItem* dataItem) const
{
}

void VolumeRenderer::renderViewPerpendicular(VolumeRenderer::DataItem* dataItem,const Vector& viewDirection) const
{
    /* Calculate the corners' parameters along the viewing direction: */
    double cornerD[8];
    for(int i=0;i<8;++i)
        cornerD[i]=corners[i].position*viewDirection;
    
    /* Find the box's distance range and the farthest away corner: */
    int maxCorner=0;
    double minD=cornerD[0];
    double maxD=cornerD[0];
    for(int i=1;i<8;++i)
    {
        if(minD>cornerD[i])
            minD=cornerD[i];
        else if(maxD<cornerD[i])
        {
            maxD=cornerD[i];
            maxCorner=i;
        }
    }
    
    /* Calculate the distance of the farthest slice: */
    double sliceOffset=sliceCenter*viewDirection;
    double sliceD=floor((maxD-sliceOffset)/sliceDistance)*sliceDistance+sliceOffset;
    
    /* Initialize the list of active edges: */
    ActiveEdge edges[12];
    ActiveEdge* firstEdge=&edges[0];
    ActiveEdge* nextEdge=&edges[0];
    PriorityHeap<EdgeExpiration,EdgeExpiration> expirations(6);
    for(int i=0;i<3;++i,++nextEdge)
    {
        /* Initialize the edge: */
        nextEdge->expired=false;
        nextEdge->startIndex=maxCorner;
        int endCorner=corners[maxCorner].neighbours[i];
        nextEdge->endIndex=endCorner;
        double rangeD=cornerD[endCorner]-cornerD[maxCorner];
        if(rangeD!=0.0)
        {
            nextEdge->dPoint=(corners[endCorner].position-corners[maxCorner].position)/rangeD;
            nextEdge->point=corners[maxCorner].position+nextEdge->dPoint*(sliceD-cornerD[maxCorner]);
            nextEdge->dPoint*=sliceDistance;
            nextEdge->dTexture=(corners[endCorner].texture-corners[maxCorner].texture)/rangeD;
            nextEdge->texture=corners[maxCorner].texture+nextEdge->dTexture*(sliceD-cornerD[maxCorner]);
            nextEdge->dTexture*=sliceDistance;
        }
        nextEdge->pred=&edges[(i+2)%3];
        nextEdge->succ=&edges[(i+1)%3];
        
        /* Store its expiration distance: */
        expirations.insert(EdgeExpiration(cornerD[endCorner],nextEdge));
    }
    
    /* Create/delete the texture cache if necessary: */
    if(textureCachingEnabled)
        dataItem->updateTextureCache(this,-1);
    else
        dataItem->deleteTextureCache();
    
    /* Set up OpenGL texturing parameters: */
    prepareRenderViewPerpendicular(dataItem);
    
    #if 1
    /* Upload the block texture: */
    if(textureCachingEnabled)
    {
        glBindTexture(GL_TEXTURE_3D,dataItem->textureObjectIDs[0]);
        if(!dataItem->textureCacheValid)
            uploadTexture3D(dataItem);
    }
    else
        uploadTexture3D(dataItem);
    #endif
    
    /* Generate slices while updating the active edge list: */
    while(sliceD>minD)
    {
        /* Process all expired edges: */
        while(expirations.getSmallest().endD>=sliceD)
        {
            /* Distinguish the four expiration cases: */
            ActiveEdge* edge=expirations.getSmallest().edge;
            int startIndex=edge->endIndex;
            if(edge->expired)
            {
                /* Edge has already expired; just remove it from the expiration queue: */
                expirations.removeSmallest();
            }
            else if(startIndex!=edge->pred->endIndex&&startIndex!=edge->succ->endIndex)
            {
                /* Split the edge: */
                edge->expired=true;
                
                /* Create the two new edges: */
                nextEdge->expired=false;
                nextEdge->startIndex=startIndex;
                int endIndex1=corners[startIndex].incomingEdgeSuccessors[edge->startIndex];
                nextEdge->endIndex=endIndex1;
                double rangeD1=cornerD[endIndex1]-cornerD[startIndex];
                if(rangeD1!=0.0)
                {
                    nextEdge->dPoint=(corners[endIndex1].position-corners[startIndex].position)/rangeD1;
                    nextEdge->point=corners[startIndex].position+nextEdge->dPoint*(sliceD-cornerD[startIndex]);
                    nextEdge->dPoint*=sliceDistance;
                    nextEdge->dTexture=(corners[endIndex1].texture-corners[startIndex].texture)/rangeD1;
                    nextEdge->texture=corners[startIndex].texture+nextEdge->dTexture*(sliceD-cornerD[startIndex]);
                    nextEdge->dTexture*=sliceDistance;
                }
                nextEdge->pred=edge->pred;
                nextEdge->pred->succ=nextEdge;
                nextEdge->succ=nextEdge+1;
                expirations.getSmallest().endD=cornerD[endIndex1];
                expirations.getSmallest().edge=nextEdge;
                expirations.reinsertSmallest();
                ++nextEdge;
                nextEdge->expired=false;
                nextEdge->startIndex=startIndex;
                int endIndex2=corners[startIndex].incomingEdgeSuccessors[endIndex1];
                nextEdge->endIndex=endIndex2;
                double rangeD2=cornerD[endIndex2]-cornerD[startIndex];
                if(rangeD2!=0.0)
                {
                    nextEdge->dPoint=(corners[endIndex2].position-corners[startIndex].position)/rangeD2;
                    nextEdge->point=corners[startIndex].position+nextEdge->dPoint*(sliceD-cornerD[startIndex]);
                    nextEdge->dPoint*=sliceDistance;
                    nextEdge->dTexture=(corners[endIndex2].texture-corners[startIndex].texture)/rangeD2;
                    nextEdge->texture=corners[startIndex].texture+nextEdge->dTexture*(sliceD-cornerD[startIndex]);
                    nextEdge->dTexture*=sliceDistance;
                }
                nextEdge->pred=nextEdge-1;
                nextEdge->succ=edge->succ;
                nextEdge->succ->pred=nextEdge;
                firstEdge=nextEdge;
                expirations.insert(EdgeExpiration(cornerD[endIndex2],nextEdge));
                ++nextEdge;
            }
            else
            {
                /* Merge the edge with one of its neighbours: */
                ActiveEdge* pred;
                ActiveEdge* succ;
                if(startIndex==edge->pred->endIndex)
                {
                    /* Merge with the clockwise neighbour: */
                    pred=edge->pred;
                    succ=edge;
                }
                else
                {
                    /* Merge with the counter-clockwise neighbour: */
                    pred=edge;
                    succ=edge->succ;
                }
                pred->expired=true;
                succ->expired=true;

                /* Create the new edge: */
                nextEdge->expired=false;
                nextEdge->startIndex=startIndex;
                int endIndex=corners[startIndex].incomingEdgeSuccessors[pred->startIndex];
                nextEdge->endIndex=endIndex;
                double rangeD=cornerD[endIndex]-cornerD[startIndex];
                if(rangeD!=0.0)
                {
                    nextEdge->dPoint=(corners[endIndex].position-corners[startIndex].position)/rangeD;
                    nextEdge->point=corners[startIndex].position+nextEdge->dPoint*(sliceD-cornerD[startIndex]);
                    nextEdge->dPoint*=sliceDistance;
                    nextEdge->dTexture=(corners[endIndex].texture-corners[startIndex].texture)/rangeD;
                    nextEdge->texture=corners[startIndex].texture+nextEdge->dTexture*(sliceD-cornerD[startIndex]);
                    nextEdge->dTexture*=sliceDistance;
                }
                nextEdge->pred=pred->pred;
                nextEdge->pred->succ=nextEdge;
                nextEdge->succ=succ->succ;
                nextEdge->succ->pred=nextEdge;
                firstEdge=nextEdge;
                expirations.getSmallest().endD=cornerD[endIndex];
                expirations.getSmallest().edge=nextEdge;
                expirations.reinsertSmallest();
                ++nextEdge;
            }
        }

        /* Generate the current polygon: */
        #if 1
        glBegin(GL_POLYGON);
        ActiveEdge* ePtr=firstEdge;
        do
        {
            glTexCoord3dv(ePtr->texture.getComponents());
            ePtr->texture-=ePtr->dTexture;
            glVertex3dv(ePtr->point.getComponents());
            ePtr->point-=ePtr->dPoint;
            ePtr=ePtr->succ;
        }
        while(ePtr!=firstEdge);
        glEnd();
        #else
        ActiveEdge* ePtr=firstEdge;
        do
        {
            ePtr->texture-=ePtr->dTexture;
            ePtr->point-=ePtr->dPoint;
            ePtr=ePtr->succ;
        }
        while(ePtr!=firstEdge);
        #endif
        
        /* Go to the next slice: */
        sliceD-=sliceDistance;
    }
    
    #if 1
    if(textureCachingEnabled)
    {
        /* Unbind the texture to prevent someone else from tampering with it: */
        glBindTexture(GL_TEXTURE_3D,0);
        
        /* Validate the texture cache: */
        dataItem->textureCacheValid=true;
    }
    #endif
}

void VolumeRenderer::calcIncrements(void)
{
    increments[2]=1;
    if(rowLength==0)
        increments[1]=increments[2]*(size[2]+borderSize+borderSize);
    else
        increments[1]=rowLength;
    if(imageHeight==0)
        increments[0]=increments[1]*(size[1]+borderSize+borderSize);
    else
        increments[0]=imageHeight;
}

VolumeRenderer::Voxel* VolumeRenderer::createPrivateMemoryBlock(const int newSize[3],int newBorderSize)
{
    /* Calculate the private block's specification: */
    privateData=true;
    for(int i=0;i<3;++i)
        size[i]=newSize[i];
    borderSize=newBorderSize;
    rowLength=0;
    imageHeight=0;
    
    int numVoxels=1;
    #ifdef __SGI_IRIX__
    /**************************************************************
    The most fucking ugly workaround ever committed by mankind:
    Calculate the texture size OpenGL will be using ahead of time,
    and allocate a user memory array of that size. Only fill in the
    used part of the array though, and set up "fake" increments to
    fool the other parts of this class into working correctly. If
    this ain't gonna win the golden raspberry, I don't know what.
    --
    In retrospect, this actually works incredibly well.
    **************************************************************/
    
    /* Calculate the fake image size: */
    int imageSize[3];
    for(int i=0;i<3;++i)
    {
        int bSize=size[i]+borderSize+borderSize;
        for(imageSize[i]=1;imageSize[i]<bSize;imageSize[i]+=imageSize[i])
            ;
        numVoxels*=imageSize[i];
    }
    increments[2]=1;
    increments[1]=imageSize[2];
    increments[0]=imageSize[1]*imageSize[2];
    #else
    calcIncrements();
    for(int i=0;i<3;++i)
    {
        int bSize=size[i]+borderSize+borderSize;
        numVoxels*=bSize;
    }
    #endif
    Voxel* result=new Voxel[numVoxels];
    values=result;
    return result;
}

void VolumeRenderer::initBoxStructure(void)
{
    /* Construct the box's connectivity: */
    static const int cornerNeighbours[8][3]={{1,2,4},{0,5,3},{0,3,6},{1,7,2},{0,6,5},{1,4,7},{2,7,4},{3,5,6}};
    for(int i=0;i<8;++i)
    {
        for(int j=0;j<3;++j)
        {
            corners[i].neighbours[j]=cornerNeighbours[i][j];
            corners[i].incomingEdgeSuccessors[cornerNeighbours[i][j]]=cornerNeighbours[i][(j+1)%3];
        }
    }
}

void VolumeRenderer::updateVoxelBlock(void)
{
    for(int i=0;i<3;++i)
    {
        /* Calculate the number of cells: */
        numCells[i]=size[i];
        if(alignment==VERTEX_CENTERED)
            --numCells[i];
        
        /* Reset the subblock selection: */
        subOrigin[i]=0;
        subSize[i]=numCells[i];
        
        /* Calculate the texture image size (always a power of two): */
        for(textureSize[i]=1;textureSize[i]<size[i];textureSize[i]+=textureSize[i])
            ;
    }
    
    /* Update other settings depending on the voxel block size: */
    calcBoxTexCoords();
    calcBoxGeometry();
    calcSlicingParameters();
    
    /* Update the data version counter: */
    ++dataVersion;
}

void VolumeRenderer::calcBoxTexCoords(void)
{
    int iMask=0x1;
    for(int i=0;i<3;++i,iMask+=iMask)
    {
        double texMin,texMax;
        if(alignment==CELL_CENTERED)
        {
            texMin=double(subOrigin[i])/double(textureSize[i]);
            texMax=double(subOrigin[i]+subSize[i])/double(textureSize[i]);
        }
        else
        {
            texMin=(double(subOrigin[i])+0.5)/double(textureSize[i]);
            texMax=(double(subOrigin[i]+subSize[i])+0.5)/double(textureSize[i]);
        }
        
        /* Update the box's texture coordinates: */
        for(int j=0;j<8;++j)
            corners[j].texture[2-i]=j&iMask?texMax:texMin;
    }
}

void VolumeRenderer::calcBoxGeometry(void)
{
    /* Calculate the corner positions in model coordinates: */
    int iMask=0x1;
    for(int i=0;i<3;++i,iMask+=iMask)
    {
        double coordMin=origin[i]+double(subOrigin[i])*extent[i]/double(numCells[i]);
        double coordMax=origin[i]+double(subOrigin[i]+subSize[i])*extent[i]/double(numCells[i]);
        
        /* Update the box's corner coordinates: */
        for(int j=0;j<8;++j)
            corners[j].position[i]=j&iMask?coordMax:coordMin;
    }
}

void VolumeRenderer::calcSlicingParameters(void)
{
    /* Calculate the minimal cell side length: */
    minCellSize=fabs(extent[0]/double(numCells[0]));
    for(int i=1;i<3;++i)
    {
        double cellSize=fabs(extent[i]/double(numCells[i]));
        if(minCellSize>cellSize)
            minCellSize=cellSize;
    }
    
    /* Calculate the slicing distance for view-perpendicular rendering: */
    sliceDistance=minCellSize*sliceFactor;
}

bool VolumeRenderer::init(void)
{
    #ifdef _WIN32
    /* Load OpenGL function pointers from DLL: */
    glTexImage3D=(PFNGLTEXIMAGE3DPROC)wglGetProcAddress("glTexImage3D");
    glTexSubImage3D=(PFNGLTEXSUBIMAGE3DPROC)wglGetProcAddress("glTexSubImage3D");
    #endif

    /* Determine the default slicing mode: */
    /* ... */

    /* Check for applicable OpenGL extensions: */
    /* ... */

    return true;
}

VolumeRenderer::VolumeRenderer(void) :
    privateData (false),
    values (0),
    borderSize (0),
    borderValue (0),
    alignment (CELL_CENTERED),
    rowLength (0),
    imageHeight (0),
    origin (0.0, 0.0, 0.0),
    minCellSize (0.0),
    renderingMode (AXIS_ALIGNED),
    interpolationMode (GL_NEAREST),
    textureFunction (GL_REPLACE),
    sliceFactor (0.5),
    sliceDistance (1.0),
    autosaveGLState (true),
    textureCachingEnabled (false),
    dataVersion (0),
    settingsVersion (0)
{
    for ( uint i = 0; i < 3; ++i )
    {
        // give everything a defined initial value
        size[i] = numCells[i] = increments[i] = textureSize[i] = 0;
        subOrigin[i] = subSize[i] = 0;
        extent[i] = 0.0;
    }
    initBoxStructure();
}

VolumeRenderer::~VolumeRenderer(void)
{
    deletePrivateData();
}

void VolumeRenderer::clearVoxelBlock(void)
{
    deletePrivateData();
}

VolumeRenderer::Voxel* VolumeRenderer::createVoxelBlock(const int newSize[3],int newBorderSize,VoxelAlignment newAlignment,int blockIncrements[3])
{
    deletePrivateData();
    
    /* Create the private memory block: */
    Voxel* result=createPrivateMemoryBlock(newSize,newBorderSize);
    
    /* Set the private block's specification: */
    alignment=newAlignment;
    
    /* Return the resulting block: */
    for(int i=0;i<3;++i)
        blockIncrements[i]=increments[i];
    return result;
}

void VolumeRenderer::finishVoxelBlock(void)
{
    /* Update other data depending on the block specification: */
    updateVoxelBlock();
}

void VolumeRenderer::setVoxelBlock(const VolumeRenderer::Voxel* newValues,const int newSize[3],int newBorderSize,VolumeRenderer::VoxelAlignment newAlignment)
{
    deletePrivateData();
    
    /* Use the given array as non-private data: */
    privateData=false;
    values=newValues;
    
    /* Copy the given specifications: */
    for(int i=0;i<3;++i)
        size[i]=newSize[i];
    borderSize=newBorderSize;
    alignment=newAlignment;
    calcIncrements();
    
    /* Update other data depending on the block specification: */
    updateVoxelBlock();
}

void VolumeRenderer::setVoxelBlock(const VolumeRenderer::Voxel* newValues,const int newSize[3],int newBorderSize,const int newIncrements[3],VolumeRenderer::VoxelAlignment newAlignment)
{
    deletePrivateData();
    
    /* Create the private memory block: */
    Voxel* ownValues=createPrivateMemoryBlock(newSize,newBorderSize);
    
    /* Set the private block's specification: */
    alignment=newAlignment;
    
    /* Copy all source values: */
    int bSize[3];
    for(int i=0;i<3;++i)
        bSize[i]=size[i]+borderSize+borderSize;
    const Voxel* svPlane=newValues;
    for(int x=0;x<bSize[0];++x,svPlane+=newIncrements[0])
    {
        const Voxel* svRow=svPlane;
        for(int y=0;y<bSize[1];++y,svRow+=newIncrements[1])
        {
            const Voxel* svPtr=svRow;
            for(int z=0;z<bSize[2];++z,svPtr+=newIncrements[2])
                ownValues[x*increments[0]+y*increments[1]+z*increments[2]]=*svPtr;
        }
    }
    
    /* Update other data depending on the block specification: */
    updateVoxelBlock();
}

void VolumeRenderer::setRowLength(int newRowLength)
{
    rowLength=newRowLength;
    calcIncrements();
}

void VolumeRenderer::setImageHeight(int newImageHeight)
{
    imageHeight=newImageHeight;
    calcIncrements();
}

void VolumeRenderer::initContext(GLContextData& contextData) const
{
    contextData.addDataItem(this,new DataItem);
}

void VolumeRenderer::clearContext(GLContextData& contextData)
{
    /* Get a pointer to the context data: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);

    /* At this point we should also remove the context data entry itself from the context! */
    contextData.removeDataItem(this);

    /* Delete all associated texture objects: */
    delete dataItem;
}

void VolumeRenderer::setGLState(GLContextData& contextData) const
{
    /* Set up the OpenGL state: */
    glPushAttrib(GL_COLOR_BUFFER_BIT|GL_CURRENT_BIT|GL_DEPTH_BUFFER_BIT|GL_ENABLE_BIT|GL_POLYGON_BIT|GL_TEXTURE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE,GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    switch(renderingMode)
    {
        case AXIS_ALIGNED:
            glEnable(GL_TEXTURE_2D);
            break;

        case VIEW_PERPENDICULAR:
            glEnable(GL_TEXTURE_3D);
            break;
    }
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,textureFunction);
}

void VolumeRenderer::resetGLState(GLContextData& contextData) const
{
    /* Restore the OpenGL state: */
    glPopAttrib();
}

VolumeRenderer::Vector VolumeRenderer::calcViewDirection(void)
{
    /* Retrieve the viewing direction in model coordinates: */
    GLfloat tempArray[4][4];
    glGetFloatv(GL_PROJECTION_MATRIX,(GLfloat*)tempArray);
    Transformation pmv (Matrix::fromColumnMajor((GLfloat*)tempArray));
//    Transformation pmv(tempArray);
    glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat*)tempArray);
    pmv *= Transformation(Matrix::fromColumnMajor((GLfloat*)tempArray));
//    pmv*=Transformation(tempArray);
    #if 1
    HVector x=pmv.inverseTransform(HVector(1.0,0.0,0.0,0.0));
    HVector y=pmv.inverseTransform(HVector(0.0,1.0,0.0,0.0));
    Vector viewDirection=cross(y.toVector(),x.toVector());
    #else
    HVector eye=pmv.inverseTransform(HVector(0.0,0.0,1.0,0.0));
    Vector viewDirection;
    if(eye.w!=0.0)
    {
        /* Perspective projection, the eye is located at an affine point: */
        viewDirection=(origin+Vector(extent)*0.5)-eye.toVertex();
    }
    else
    {
        /* Parallel projection, the eye is at infinity: */
        viewDirection=eye.toVector();
    }
    #endif
    viewDirection.normalize();

    return viewDirection;
}

void VolumeRenderer::renderBlock(GLContextData& contextData) const
{
    /* Render the block with automatically calculated viewing direction: */
    renderBlock(contextData,calcViewDirection());
}

void VolumeRenderer::renderBlock(GLContextData& contextData,const Vector& viewDirection) const
{
    /* Get a pointer to the context data: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    
    /* Render the voxel block using the current rendering mode: */
    if(autosaveGLState)
        setGLState(contextData);
    switch(renderingMode)
    {
        case AXIS_ALIGNED:
            renderAxisAligned(dataItem,viewDirection);
            break;
        
        case VIEW_PERPENDICULAR:
            renderViewPerpendicular(dataItem,viewDirection);
            break;
    }
    if(autosaveGLState)
        resetGLState(contextData);
}
