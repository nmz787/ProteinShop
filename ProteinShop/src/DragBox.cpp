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
DragBox - Class for 6-DOF rigid body movement interaction.
***********************************************************************/

#include <Math/Math.h>
#include <Math/Constants.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLTypes.h>
#include <GLGeometry.h>

#include "DragBox.h"

/************************
Methods of class DragBox:
************************/

double DragBox::intersectVertex(int vertexIndex,const DragBox::Point& start,const DragBox::Vector& direction) const
{
    /* Construct the vertex' position: */
    Point p=center;
    for(int i=0;i<3;++i)
    {
        if(vertexIndex&(1<<i))
            p+=axis[i]*(size[i]*0.5);
        else
            p-=axis[i]*(size[i]*0.5);
    }
    
    /* Intersect the ray with the vertex sphere: */
    double denominator=Geometry::sqr(direction);
    if(denominator!=0.0)
    {
        Vector sp=start-p;
        double spd=sp*direction;
        double det=Math::sqr(spd)-(Geometry::sqr(sp)-Math::sqr(vertexRadius))*denominator;
        if(det>=0.0)
        {
            det=Math::sqrt(det);
            double lambda=(-spd-det)/denominator;
            if(lambda<0.0)
                lambda=(-spd+det)/denominator;
            if(lambda>=0.0)
                return lambda;
            else
                return Math::Constants<double>::max;
        }
        else
            return Math::Constants<double>::max;
    }
    else
        return Math::Constants<double>::max;
}

double DragBox::intersectEdge(int axisIndex,int vertexMask,const DragBox::Point& start,const DragBox::Vector& direction) const
{
    /* Construct the edge's start point: */
    Point p=center;
    vertexMask&=~(1<<axisIndex);
    for(int i=0;i<3;++i)
    {
        if(vertexMask&(1<<i))
            p+=axis[i]*(size[i]*0.5);
        else
            p-=axis[i]*(size[i]*0.5);
    }
    
    /* Intersect the ray with the edge cylinder: */
    double da=direction*axis[axisIndex];
    double denominator=Geometry::sqr(direction)-Math::sqr(da);
    if(denominator!=0.0)
    {
        Vector sp=start-p;
        double spa=sp*axis[axisIndex];
        double ph=sp*direction-spa*da;
        double q=Geometry::sqr(sp)-Math::sqr(spa)-Math::sqr(edgeRadius);
        double determinant=Math::sqr(ph)-q*denominator;
        if(determinant>=0.0)
        {
            /* Calculate the first intersection: */
            determinant=Math::sqrt(determinant);
            double lambda=(-ph-determinant)/denominator;
            double beta=spa+da*lambda;
            if(lambda<0.0||beta<0.0||beta>size[axisIndex])
            {
                /* Calculate the second intersection: */
                lambda=(-ph+determinant)/denominator;
                beta=spa+da*lambda;
            }
            if(lambda>=0.0&&beta>=0.0&&beta<=size[axisIndex])
                return lambda; // We've got a valid intersection!
            else
                return Math::Constants<double>::max; // Intersection is not inside edge
        }
        else
            return Math::Constants<double>::max; // Ray does not intersect edge's cylinder
    }
    else
        return Math::Constants<double>::max; // Ray is parallel to edge's cylinder axis (ignore cylinder caps here)
}

double DragBox::intersectFace(int axisIndex,int faceSign,const DragBox::Point& start,const DragBox::Vector& direction) const
{
    /* Intersect the ray and the face's plane: */
    double denominator=direction*axis[axisIndex];
    if(denominator!=0.0)
    {
        /* Calculate the face's plane offset: */
        double dist=axis[axisIndex]*center;
        if(faceSign>0)
            dist+=size[axisIndex]*0.5;
        else
            dist-=size[axisIndex]*0.5;
        
        /* Calculate the intersection point: */
        double lambda=(dist-start*axis[axisIndex])/denominator;
        Point p=start+direction*lambda;
        
        /* Check the intersection point for validity: */
        Vector pc=p-center;
        for(int i=0;i<3;++i)
        {
            if(i!=axisIndex)
            {
                double offset=pc*axis[i];
                if(offset<-size[i]*0.5||offset>size[i]*0.5)
                    return Math::Constants<double>::max; // Point is outside face
            }
        }
        return lambda; // We've got a valid intersection!
    }
    else
        return Math::Constants<double>::max; // Ray is parallel to face's plane
}

DragBox::DragBox(void)
    :center(0.0,0.0,0.0),edgeRadius(0.0),vertexRadius(0.0),
     currentPickMode(TRANSPARENT),currentDragMode(NONE),
     rotateCenter(center)
{
    axis[0]=Vector(1.0,0.0,0.0);
    axis[1]=Vector(0.0,1.0,0.0);
    axis[2]=Vector(0.0,0.0,1.0);
    size[0]=1.0;
    size[1]=1.0;
    size[2]=1.0;
    for(int i=0;i<12;++i)
        edgeHighlightFlags[i]=false;
}

void DragBox::setCenter(const DragBox::Point& newCenter)
{
    center=newCenter;
    rotateCenter=newCenter;
}

void DragBox::setAxis(int axisIndex,const DragBox::Vector& newAxis)
{
    axis[axisIndex]=newAxis;
    axis[axisIndex].normalize();
}

void DragBox::setSize(int axisIndex,double newSize)
{
    size[axisIndex]=newSize;
}

void DragBox::setEdgeRadius(double newEdgeRadius)
{
    edgeRadius=newEdgeRadius;
    vertexRadius=1.5*newEdgeRadius;
}

void DragBox::setRotateCenter(const DragBox::Point& newRotateCenter)
{
    rotateCenter=newRotateCenter;
    if(currentDragMode==ROTATING_AXIS)
    {
        /* Have to update cylinder radius: */
        Vector lc=lastIntersection-rotateCenter;
        rotateCylinderRadius=Math::sqrt(Geometry::sqr(lc)-Math::sqr(lc*rotateAxis));
    }
}

bool DragBox::pickIncremental(void)
{
    /* Prepare for incremental 6-DOF dragging: */
    currentDragMode=SIXDOF;
    
    /* Highlight all edges: */
    for(int i=0;i<12;++i)
        edgeHighlightFlags[i]=true;

    /* Initialize the dragging transformation: */
    dragTransformation=Transformation::identity;
    
    return true;
}

bool DragBox::pick(const DragBox::Transformation& transformation)
{
    /* Initialize intersection result: */
    currentDragMode=NONE;
    
    /* Check if the transformation's origin is inside the box: */
    Vector pc=transformation.getOrigin()-center;
    bool inside=true;
    for(int i=0;i<3&&inside;++i)
    {
        double offset=pc*axis[i];
        inside=offset>=-size[i]*0.5&&offset<=size[i]*0.5;
    }
    
    if(inside)
    {
        currentDragMode=SIXDOF;
        
        /* Highlight all edges: */
        for(int i=0;i<12;++i)
            edgeHighlightFlags[i]=true;
        
        /* Calculate the pre-transformation for subsequent dragging: */
        initialTransformation=Geometry::invert(transformation);
        dragTransformation=Transformation::identity;
    }
    
    return inside;
}

bool DragBox::pick(const DragBox::Point& start,const DragBox::Vector& direction)
{
    /* Initialize intersection result: */
    double lambdaMin=Math::Constants<double>::max;
    currentDragMode=NONE;
    
    /* Intersect ray with all edges: */
    int rotateAxisIndex=-1;
    static const int edgeIndicators[12][2]={{0,0x0},{0,0x2},{0,0x6},{0,0x4},
                                            {1,0x0},{1,0x4},{1,0x5},{1,0x1},
                                            {2,0x0},{2,0x1},{2,0x3},{2,0x2}};
    for(int edge=0;edge<12;++edge)
    {
        double lambda=intersectEdge(edgeIndicators[edge][0],edgeIndicators[edge][1],start,direction);
        if(lambdaMin>lambda)
        {
            lambdaMin=lambda;
            currentDragMode=ROTATING_AXIS;
            lastIntersection=start+direction*lambda;
            rotateAxisIndex=edgeIndicators[edge][0];
            rotateAxis=axis[rotateAxisIndex];
            Vector lc=lastIntersection-rotateCenter;
            rotateCylinderRadius=Math::sqrt(Geometry::sqr(lc)-Math::sqr(lc*rotateAxis));
        }
    }
    
    /* Intersect ray with all faces: */
    int translateFaceIndex=-1;
    static const int faceIndicators[6][2]={{0,-1},{0,1},{1,-1},{1,1},{2,-1},{2,1}};
    if(currentPickMode==OPAQUE||currentDragMode==NONE)
    {
        for(int face=0;face<6;++face)
        {
            double lambda=intersectFace(faceIndicators[face][0],faceIndicators[face][1],start,direction);
            if(lambdaMin>lambda)
            {
                lambdaMin=lambda;
                currentDragMode=TRANSLATING;
                lastIntersection=start+direction*lambda;
                translateFaceIndex=face;
                translateFaceNormal=axis[faceIndicators[face][0]];
                if(faceIndicators[face][1]<0)
                    translateFaceNormal*=-1.0;
                translateFaceOffset=lastIntersection*translateFaceNormal;
            }
        }
    }
    
    /* Hightlight picked edge(s): */
    if(currentDragMode==TRANSLATING)
    {
        static const int faceEdgeIndices[6][4]={{4,6,8,10},{5,7,9,11},{0,2,8,9},{1,3,10,11},{0,1,4,5},{2,3,6,7}};
        for(int i=0;i<4;++i)
            edgeHighlightFlags[faceEdgeIndices[translateFaceIndex][i]]=true;
    }
    else if(currentDragMode==ROTATING_AXIS)
    {
        static const int axisEdgeIndices[3][4]={{0,1,2,3},{4,5,6,7},{8,9,10,11}};
        for(int i=0;i<4;++i)
            edgeHighlightFlags[axisEdgeIndices[rotateAxisIndex][i]]=true;
    }
    
    /* Return intersection result: */
    if(currentDragMode!=NONE)
        dragTransformation=Transformation::identity;
    return currentDragMode!=NONE;
}

bool DragBox::pick(const DragBox::HTransformation& modelView,const DragBox::HTransformation& projection,const DragBox::Point& mouseClip)
{
    /* Create ray by reprojection mouse position into model coordinates: */
    Point start=modelView.inverseTransform(projection.inverseTransform(Point(mouseClip[0],mouseClip[1],-1.0)));
    Vector direction=modelView.inverseTransform(projection.inverseTransform(Point(mouseClip[0],mouseClip[1],1.0)))-start;
    
    /* Initialize intersection result: */
    double lambdaMin=Math::Constants<double>::max;
    currentDragMode=NONE;
    
    /* Intersect ray with all vertices: */
    int rotateVertexIndex=-1;
    for(int vertex=0;vertex<8;++vertex)
    {
        double lambda=intersectVertex(vertex,start,direction);
        if(lambdaMin>lambda)
        {
            lambdaMin=lambda;
            currentDragMode=ROTATING_VERTEX;
            lastMouse=mouseClip;
	    rotateVertexIndex=vertex;
        }
    }
    
    /* Intersect ray with all edges: */
    int rotateAxisIndex=-1;
    static const int edgeIndicators[12][2]={{0,0x0},{0,0x2},{0,0x6},{0,0x4},
                                            {1,0x0},{1,0x4},{1,0x5},{1,0x1},
                                            {2,0x0},{2,0x1},{2,0x3},{2,0x2}};
    if(currentPickMode==OPAQUE||currentDragMode==NONE)
    {
        for(int edge=0;edge<12;++edge)
        {
            double lambda=intersectEdge(edgeIndicators[edge][0],edgeIndicators[edge][1],start,direction);
            if(lambdaMin>lambda)
            {
                lambdaMin=lambda;
                currentDragMode=ROTATING_AXIS;
                lastMouse=mouseClip;
                rotateAxisIndex=edgeIndicators[edge][0];
                rotateAxis=axis[rotateAxisIndex];
            }
        }
    }
    
    /* Intersect ray with all faces: */
    int translateFaceIndex=-1;
    static const int faceIndicators[6][2]={{0,-1},{0,1},{1,-1},{1,1},{2,-1},{2,1}};
    if(currentPickMode==OPAQUE||currentDragMode==NONE)
    {
        for(int face=0;face<6;++face)
        {
            double lambda=intersectFace(faceIndicators[face][0],faceIndicators[face][1],start,direction);
            if(lambdaMin>lambda)
            {
                lambdaMin=lambda;
                currentDragMode=TRANSLATING;
                lastIntersection=start+direction*lambda;
                translateFaceIndex=face;
                translateFaceNormal=axis[faceIndicators[face][0]];
                if(faceIndicators[face][1]<0)
                    translateFaceNormal*=-1.0;
                translateFaceOffset=lastIntersection*translateFaceNormal;
            }
        }
    }
    
    /* Save z distance between picked point and center of rotation: */
    if(currentDragMode==ROTATING_VERTEX||currentDragMode==ROTATING_AXIS)
    {
        Point rotateCenterClip=projection.transform(modelView.transform(rotateCenter));
        Point intersectionClip=projection.transform(modelView.transform(start+lambdaMin*direction));
        rotateDepth=rotateCenterClip[2]-intersectionClip[2];
    }
    
    /* Hightlight picked edge(s): */
    if(currentDragMode==TRANSLATING)
    {
        static const int faceEdgeIndices[6][4]={{4,6,8,10},{5,7,9,11},{0,2,8,9},{1,3,10,11},{0,1,4,5},{2,3,6,7}};
        for(int i=0;i<4;++i)
            edgeHighlightFlags[faceEdgeIndices[translateFaceIndex][i]]=true;
    }
    else if(currentDragMode==ROTATING_AXIS)
    {
        static const int axisEdgeIndices[3][4]={{0,1,2,3},{4,5,6,7},{8,9,10,11}};
        for(int i=0;i<4;++i)
            edgeHighlightFlags[axisEdgeIndices[rotateAxisIndex][i]]=true;
    }
    
    /* Return intersection result: */
    if(currentDragMode!=NONE)
        dragTransformation=Transformation::identity;
    return currentDragMode!=NONE;
}

void DragBox::dragIncremental(const DragBox::Transformation& transformation)
{
    dragTransformation.leftMultiply(Transformation::translate(Point::origin-rotateCenter));
    dragTransformation.leftMultiply(transformation);
    dragTransformation.leftMultiply(Transformation::translate(rotateCenter-Point::origin));
}

void DragBox::drag(const DragBox::Transformation& transformation)
{
    dragTransformation=transformation*initialTransformation;
}

void DragBox::drag(const DragBox::Point& start,const DragBox::Vector& direction)
{
    Transformation deltaT=Transformation::identity;
    switch(currentDragMode)
    {
        case TRANSLATING:
        {
            /* Intersect the new ray with the translating plane: */
            double denominator=direction*translateFaceNormal;
            if(denominator!=0.0)
            {
                double lambda=(translateFaceOffset-start*translateFaceNormal)/denominator;
                if(lambda>=0.0)
                {
                    /* Calculate the intersection point: */
                    Point newIntersection=start+direction*lambda;

                    /* Calculate the result translation: */
                    deltaT=Transformation::translate(newIntersection-lastIntersection);
                    lastIntersection=newIntersection;
                }
            }
            break;
        }
        
        case ROTATING_AXIS:
        {
            /* Intersect the new ray with the rotating cylinder: */
            double da=direction*rotateAxis;
            double denominator=Geometry::sqr(direction)-Math::sqr(da);
            if(denominator!=0.0)
            {
                Vector sp=start-rotateCenter;
                double spa=sp*rotateAxis;
                double ph=sp*direction-spa*da;
                double q=Geometry::sqr(sp)-Math::sqr(spa)-Math::sqr(rotateCylinderRadius);
                double determinant=Math::sqr(ph)-q*denominator;
                if(determinant>=0.0)
                {
                    /* Calculate the first intersection: */
                    determinant=Math::sqrt(determinant);
                    double lambda=(-ph-determinant)/denominator;
                    if(lambda<0.0)
                    {
                        /* Calculate the second intersection: */
                        lambda=(-ph+determinant)/denominator;
                    }
                    if(lambda>=0.0)
                    {
                        /* Calculate the intersection point: */
                        Point newIntersection=start+direction*lambda;

                        /* Calculate the result rotation: */
                        Vector lc=lastIntersection-rotateCenter;
                        Vector nc=newIntersection-rotateCenter;
                        Vector normal1=Geometry::cross(rotateAxis,lc);
                        Vector normal2=Geometry::cross(rotateAxis,nc);
                        double angle=Math::acos((normal1*normal2)/Math::sqr(rotateCylinderRadius));
                        if(Geometry::cross(lc,nc)*rotateAxis<0.0)
                            angle*=-1.0;
                        deltaT=Transformation::translate(rotateCenter-Point::origin);
                        deltaT*=Transformation::rotate(Transformation::Rotation::rotateAxis(rotateAxis,angle));
                        deltaT*=Transformation::translate(Point::origin-rotateCenter);
                        lastIntersection=newIntersection;
                    }
                }
            }
            break;
        }

        case NONE:
            break;
        case ROTATING_VERTEX:
            break;
        case SIXDOF:
            break;
    }
    
    dragTransformation.leftMultiply(deltaT);
}

void DragBox::drag(const DragBox::HTransformation& modelView,const DragBox::HTransformation& projection,const DragBox::Point& mouseClip)
{
    Transformation deltaT=Transformation::identity;
    switch(currentDragMode)
    {
        case TRANSLATING:
        {
            /* Create ray by reprojection mouse position into model coordinates: */
            Point start=modelView.inverseTransform(projection.inverseTransform(Point(mouseClip[0],mouseClip[1],-1.0)));
            Vector direction=modelView.inverseTransform(projection.inverseTransform(Point(mouseClip[0],mouseClip[1],1.0)))-start;
            
            /* Intersect the new ray with the translating plane: */
            double denominator=direction*translateFaceNormal;
            if(denominator!=0.0)
            {
                double lambda=(translateFaceOffset-start*translateFaceNormal)/denominator;
                if(lambda>=0.0)
                {
                    /* Calculate the intersection point: */
                    Point newIntersection=start+direction*lambda;

                    /* Calculate the result translation: */
                    deltaT=Transformation::translate(newIntersection-lastIntersection);
                    lastIntersection=newIntersection;
                }
            }
            break;
        }
        
        case ROTATING_AXIS:
        {
            /* Project rotation axis to screen coordinates: */
            Point axisStartClip=projection.transform(modelView.transform(rotateCenter));
            Point axisEndClip=projection.transform(modelView.transform(rotateCenter+rotateAxis));
            Vector axisClip=axisEndClip-axisStartClip;
            double axisClipLen=Math::sqrt(Math::sqr(axisClip[0])+Math::sqr(axisClip[1]));
            
            /* Switch between disc and cylinder interaction depending on the axis orientation: */
            double angle=0.0;
            if(axisClipLen<Math::abs(axisClip[2])*0.001)
            {
                /* Disc interaction: */
                angle=(mouseClip[0]-lastMouse[0])*(lastMouse[1]-axisStartClip[1])-(mouseClip[1]-lastMouse[1])*(lastMouse[0]-axisStartClip[0]);
            }
            else
            {
                /* Cylinder interaction: */
                angle=((mouseClip[0]-lastMouse[0])*axisClip[1]-(mouseClip[1]-lastMouse[1])*axisClip[0])/axisClipLen;
                angle*=50.0;
                angle=0.33*angle+angle*Math::abs(angle);
                angle*=0.02;
            }
            
            /* Calculate the result rotation: */
            deltaT=Transformation::translate(rotateCenter-Point::origin);
            deltaT*=Transformation::rotate(Transformation::Rotation::rotateAxis(rotateAxis,angle));
            deltaT*=Transformation::translate(Point::origin-rotateCenter);
            lastMouse=mouseClip;
            break;
        }
        
        case ROTATING_VERTEX:
        {
            /* Project rotation center to screen coordinates: */
            Point centerClip=projection.transform(modelView.transform(rotateCenter));
            
            /* Calculate rotation axis: */
            Vector off=lastMouse-centerClip;
            off[2]=rotateDepth*0.25;
            Vector mouseDelta=mouseClip-lastMouse;
            Vector axis=Geometry::cross(off,mouseDelta);
            deltaT=Transformation::translate(rotateCenter-Point::origin);
            if(Geometry::sqr(axis)>1.0e-8)
            {
                axis=modelView.inverseTransform(projection.inverseTransform(axis));
                deltaT*=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,1.5*Geometry::mag(mouseDelta)));
            }
            deltaT*=Transformation::translate(Point::origin-rotateCenter);
            lastMouse=mouseClip;
            break;
        }

        case NONE:
            break;

        case SIXDOF:
            break;
	
    }
    
    dragTransformation.leftMultiply(deltaT);
}

void DragBox::release(void)
{
    /* Apply the current dragging transformation: */
    center=dragTransformation.transform(center);
    for(int i=0;i<3;++i)
        axis[i]=dragTransformation.transform(axis[i]);
    rotateCenter=dragTransformation.transform(rotateCenter);
    
    currentDragMode=NONE;
    
    /* Un-highlight all edges: */
    for(int i=0;i<12;++i)
        edgeHighlightFlags[i]=false;
}

void DragBox::draw(void) const
{
    Point c=center;
    Vector sx=axis[0]*size[0];
    Vector sy=axis[1]*size[1];
    Vector sz=axis[2]*size[2];
    if(currentDragMode!=NONE)
    {
        c=dragTransformation.transform(c);
        sx=dragTransformation.transform(sx);
        sy=dragTransformation.transform(sy);
        sz=dragTransformation.transform(sz);
    }
    
    Point corners[8];
    corners[0]=c-(sx+sy+sz)*0.5;
    corners[1]=corners[0]+sx;
    corners[3]=corners[1]+sy;
    corners[2]=corners[3]-sx;
    corners[6]=corners[2]+sz;
    corners[7]=corners[6]+sx;
    corners[5]=corners[7]-sy;
    corners[4]=corners[5]-sx;
    
    GLColor<GLfloat,4> baseColor=glGetColor<GLfloat>(GL_CURRENT_COLOR);
    GLColor<GLfloat,4> highlightColor(1.0,1.0,0.0);
    
    /* Draw the box's edges: */
    static const int edgeVertexIndices[12][2]={{0,1},{2,3},{4,5},{6,7},{0,2},{1,3},{4,6},{5,7},{0,4},{1,5},{2,6},{3,7}};
    glBegin(GL_LINES);
    for(int i=0;i<12;++i)
    {
        if(edgeHighlightFlags[i])
            glColor(highlightColor);
        else
            glColor(baseColor);
        glVertex(corners[edgeVertexIndices[i][0]]);
        glVertex(corners[edgeVertexIndices[i][1]]);
    }
    glEnd();
    
    /* Draw the box's vertices: */
    if(currentDragMode==ROTATING_VERTEX)
        glColor(highlightColor);
    else
        glColor(baseColor);
    glBegin(GL_POINTS);
    for(int i=0;i<8;++i)
        glVertex(corners[i]);
    glEnd();
    
    /* Render the box's faces semi-transparently: */
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    glColor(baseColor);
    glBegin(GL_TRIANGLE_FAN);
    glVertex(corners[0]);
    glVertex(corners[1]);
    glVertex(corners[5]);
    glVertex(corners[4]);
    glVertex(corners[6]);
    glVertex(corners[2]);
    glVertex(corners[3]);
    glVertex(corners[1]);
    glEnd();
    glBegin(GL_TRIANGLE_FAN);
    glVertex(corners[7]);
    glVertex(corners[6]);
    glVertex(corners[4]);
    glVertex(corners[5]);
    glVertex(corners[1]);
    glVertex(corners[3]);
    glVertex(corners[2]);
    glVertex(corners[6]);
    glEnd();
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
}
