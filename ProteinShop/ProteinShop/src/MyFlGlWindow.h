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
MyFlGlWindow - Class for rendering protein structures using OpenGL.
Implements a trackball navigation interface and hotkeys for render
toggles.

***********************************************************************/

#ifndef MYFLGLWINDOW_INCLUDED
#define MYFLGLWINDOW_INCLUDED

#include <ctime>
#include <list>

#include <Geometry/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/Ray.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/ProjectiveTransformation.h>

#include <FL/gl.h>
#include <FL/Fl_Gl_Window.H>

#include "GLTypes.h"
#include "GLLightSource.h"
#include "GLContextData.h"
#include "Globals.h"

/* Forward declarations: */
class MySpaceBall;

class MyFlGlWindow:public Fl_Gl_Window
{
    /* Embedded classes: */
    public:
    typedef Geometry::Vector<double,3> Vector;
    typedef Geometry::Point<double,3> Point;
    typedef Geometry::Ray<double,3> Ray;
    typedef Geometry::Rotation<double,3> Rotation;
    typedef Geometry::OrthonormalTransformation<double,3> MVTransformation;
    typedef Geometry::ProjectiveTransformation<double,3> Transformation;

    private:
    enum ZipMode // Direction flag for automatic beta-strand bonding
    {
        ANTIPARALLEL,PARALLEL
    };
    
    struct ViewSpec // Structure to store view specifications
    {
        /* Elements: */
        public:
        double zoom; // Zoom factor
        double dolly; // Distance from viewpoint to scene center
        MVTransformation modelView; // Modelview transformation
        
        /* Constructors and destructors: */
        ViewSpec(double sZoom,double sDolly,const MVTransformation& sModelView)
            :zoom(sZoom),dolly(sDolly),modelView(sModelView)
        {
        };
    };
    
    /* Elements: */
    private:
    
    /* GL parameters: */
    GLContextData contextData; // A GL data context for rendering
    GLColor<GLfloat,4> backgroundColor; // Background color for the GL window
    GLLightSource headlight; // The default headlight
    bool contextInitialized;
    bool projectionMatrixChanged;
    bool settime;
	time_t oldtime;
	double energyValue;
	double rmsd;
    
    /* Viewing specification: */
    double windowAspect; // Aspect ratio of rendering window in pixels
    double frontPlaneDistance,backPlaneDistance; // Distance from eye point to front/back plane
    double zoom; // Current zoom factor
    double dolly; // Current distance from viewpoint to scene center
    MVTransformation modelView; // Current modelview transformation
    std::list<ViewSpec> viewSpecStack; // Stack of stored viewing specifications
    
    /* OpenGL viewing parameters: */
    Transformation projectionMatrix,modelViewMatrix; // Local copies of OpenGL transformations
    
    /* Interaction parameters: */
    int rotateButtonMask,panButtonMask,zoomButtonMask; // Mouse button masks for viewpoint manipulation
    double winRadius; // Trackball radius
    int dragX,dragY,dragButton; // State to keep track of mouse movements inside the window
    MySpaceBall* spaceBall; // Pointer to SpaceBall input device, if one is connected
    int numZoomingSteps; // Number of zooming steps for smooth transition; doubles as flag for zooming if >0
    Vector zoomingTranslation; // Displacement vector for each zooming step
    double zoomingDollyStep; // Dolly increment/decrement for each zooming step
    bool isSpinning; // Flag if the display is spinning
    MVTransformation spinTransformation; // Last transformation increment before mouse button was released
    bool isManipulating; // Flag if a drag box is being dragged
    bool isZipping; // Flag if user wants to align two residues automatically
    ZipMode zipMode; // Current zipping mode
    
    /* Parameters for zipping animation: */
    int zipStep,numZipSteps;
    Point zipRotateCenterStart,zipRotateCenterEnd;
    Vector zipTranslation;
    Vector zipRotationAxis;
    double zipRotationAngle;
    
    /* Methods inherited from Fl_Gl_Window: */
    protected:
    int pick(int mouseX,int mouseY,int which); // Picks a secondary structure (which==0) or residue (which==1) by mouse coordinates
    virtual void draw(void);
    public:
    virtual int handle(int eventType);
    
    /* New protected methods: */
    protected:
    Ray calcMouseRay(void) const;
    MVTransformation calcRotation(int dX,int dY) const;
    static void zoomCB(void* cbData);
    static void spinCB(void* cbData);
    static void spaceBallCB(void* cbData);
    static void waitForIkCB(void* cbData);
    static void zipCB(void* cbData);
    
    /* Constructors and destructors: */
    public:
    MyFlGlWindow(int w,int h,const char* label =0);
    MyFlGlWindow(int x,int y,int w,int h,const char* label =0);
    ~MyFlGlWindow(void);
    
    /* New methods: */
    void clearContext(ProteinState *state); // Clears GLContextData for a protein prior to destruction
    void changeRenderer(void); // Notifies the window that a new protein renderer has been created
    void centerView(void); // Resets the viewing transformation
    void pushViewSpec(void); // Pushes current viewing specification onto stack
    void popViewSpec(void); // Pops topmost viewing specification from stack
    void zipIt(const MVTransformation& zipTransformation);
    void redrawNow(void); // Redraw window contents synchronously
    bool saveScreenshot(const char* fileName) const; // Saves a screenshot of the main window as a binary PPM file
	void updateEnergyVaule(void);
	void setRMSDValue(double value);
};

#endif
