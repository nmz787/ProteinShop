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

#ifndef DRAGBOX_INCLUDED
#define DRAGBOX_INCLUDED

#include <Geometry/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/ProjectiveTransformation.h>

class DragBox
{
    /* Embedded classes: */
    public:
    typedef Geometry::Vector<double,3> Vector;
    typedef Geometry::Point<double,3> Point;
    typedef Geometry::OrthonormalTransformation<double,3> Transformation;
    typedef Geometry::ProjectiveTransformation<double,3> HTransformation;
    
    enum PickMode // Enumerated type for picking modes
    {
        TRANSPARENT,OPAQUE
    };
    
    enum DragMode // Enumerated type for dragging modes
    {
        NONE,TRANSLATING,ROTATING_AXIS,ROTATING_VERTEX,SIXDOF
    };
    
    /* Elements: */
    private:
    Point center; // The drag box's center
    Vector axis[3]; // The drag box's three (normalized) axis vectors
    double size[3]; // The drag box's sizes along the axis vectors
    double edgeRadius; // The radius of an edge cylinder
    double vertexRadius; // The radius of a vertex sphere
    
    PickMode currentPickMode; // The current picking mode (determines how box is picked with 2D devices)
    DragMode currentDragMode; // The current dragging mode (gets set after a successful pick operation)
    bool edgeHighlightFlags[12]; // Highlight flags for each edge for rendering
    Point lastIntersection; // The last intersection point for dragging
    double rotateDepth; // Z distance between picked point and rotation center in clip coordinates
    Point lastMouse; // Last picked/dragged mouse position in clip coordinates
    Vector translateFaceNormal; // Normal vector of translating plane
    double translateFaceOffset; // Offset of translating plane
    Vector rotateAxis; // Axis of rotation
    Point rotateCenter; // Center of rotation
    double rotateCylinderRadius; // Radius of rotation cylinder
    Transformation initialTransformation; // Transformation matrix to pre-apply when dragging with a 6-DOF dragger
    Transformation dragTransformation; // Complete transformation matrix during dragging
    
    /* Private methods */
    double intersectVertex(int vertexIndex,const Point& start,const Vector& direction) const;
    double intersectEdge(int axisIndex,int vertexMask,const Point& start,const Vector& direction) const;
    double intersectFace(int axisIndex,int faceSign,const Point& start,const Vector& direction) const;
    
    /* Constructors and destructors: */
    public:
    DragBox(void); // Constructs a default drag box
    
    /* Methods: */
    const Point& getCenter(void) // Returns box's current center
    {
        return center;
    };
    void setCenter(const Point& newCenter); // Sets a new center point
    void setAxis(int axisIndex,const Vector& newAxis); // Sets one axis
    void setSize(int axisIndex,double newSize); // Sets one size
    void setEdgeRadius(double newEdgeRadius); // Sets a new radius for ray/edge intersection tests
    const Point& getRotateCenter(void) const // Returns current center of rotation
    {
        return rotateCenter;
    };
    void setRotateCenter(const Point& newRotateCenter); // Sets a new center of rotation
    PickMode getPickMode(void) const // Returns current picking mode
    {
        return currentPickMode;
    };
    void setPickMode(PickMode newPickMode) // Sets current picking mode
    {
        currentPickMode=newPickMode;
    };
    bool pickIncremental(void); // Prepares box for subsequent incremental 6-DOF dragging
    bool pick(const Transformation& transformation); // Returns true if the transformation's origin is inside the box; prepares for subsequent 6-DOF box dragging
    bool pick(const Point& start,const Vector& direction); // Returns true if ray intersects box; sets up internal state for subsequent dragging
    bool pick(const HTransformation& modelView,const HTransformation& projection,const Point& mouseClip); // Returns true if mouse picks box; sets up internal state for dragging
    DragMode getDragMode(void) const // Returns current dragging mode
    {
        return currentDragMode;
    };
    void dragIncremental(const Transformation& transformation); // Drags the box with an incremental 6-DOF dragger
    void drag(const Transformation& transformation); // Drags the box with a 6-DOF dragger
    void drag(const Point& start,const Vector& direction); // Drags the box via a ray intersecting the 3D widget
    void drag(const HTransformation& modelView,const HTransformation& projection,const Point& mouseClip); // Drags the box via screen-space position
    const Transformation& getDragTransformation(void) const // Returns the current dragging transformation
    {
        return dragTransformation;
    };
    void release(void); // Resets the dragging mode to NONE
    void draw(void) const; // Draws the box
};

#endif
