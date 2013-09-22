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

#include <stdio.h>
#include <Math/Math.h>
#include <GLTemplates.h>
#include <GLGeometry.h>
#include <GLTransformations.h>
#include <FL/Fl.H>

#include "MySpaceBall.h"

#include "MyFlGlWindow.h"

extern void updateRamachandranWindow( );

/*****************************
Methods of class MyFlGlWindow:
*****************************/

void MyFlGlWindow::draw(void)
{
    using namespace MD;
    typedef std::pair<double, double> Range;

    if ( !numProteins() )
    {
        // clear any afterimage of previous buffer contents
        glClearColor (backgroundColor);
        glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        return;
    }
    if(!valid())
    {
        if(!contextInitialized)
        {
            /* Set OpenGL state: */
            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LEQUAL);
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
            glFrontFace(GL_CCW);
            glEnable(GL_LIGHTING);
            glEnable(GL_NORMALIZE);
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
            for ( uint i = 0; i < numProteins(); ++i )
            {
                // initialize uninitialized ProteinRenderer(s)
                ProteinState *state = getProtein (i);
                if ( state )
                {
                    ProteinRenderer *pr = state->proteinRenderer;
                    if ( pr && !pr->isContextInitialized(contextData) )
                        pr->initContext (contextData);
                    EnergyRenderer *er = state->energyRenderer;
                    if ( er && !er->isContextInitialized(contextData) )
                        er->initContext (contextData);
                }
            }
            contextInitialized=true;
		}

        /* Calculate the window's aspect ratio: */
        windowAspect=double(w())/double(h());

        /* Calculate the window radius: */
        if(w()>=h())
            winRadius=double(w())*0.5;
        else
            winRadius=double(h())*0.5;

        /* Set the viewport: */
        glViewport(0,0,w(),h());
    }

    /* Clear the screen: */
    glClearColor(backgroundColor);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    double frontDist = dolly, backDist = dolly;
    for ( uint i = 0; i < numProteins(); ++i )
    {
        // take the maximum depth range of all proteins
        ProteinState *state = getProtein (i);
        if ( !state ) continue;
        Range depthRange (state->proteinRenderer->calcDepthRange(modelView));
        depthRange.first += dolly;
        depthRange.second += dolly;
        if ( frontDist > depthRange.first )
            frontDist = depthRange.first;
        if ( backDist < depthRange.second )
            backDist = depthRange.second;
    }
    // calculate front and back plane distances
    if ( frontDist < frontPlaneDistance )
        frontDist = frontPlaneDistance;
    if ( backDist < (frontDist + 1.0) )
        backDist = frontDist + 1.0;

    /* Set the projection matrix: */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double frontPlaneSize=frontDist*zoom;
    if(windowAspect>=1.0)
        glFrustum(-frontPlaneSize*windowAspect,frontPlaneSize*windowAspect,-frontPlaneSize,frontPlaneSize,frontDist,backDist);
    else
        glFrustum(-frontPlaneSize,frontPlaneSize,-frontPlaneSize/windowAspect,frontPlaneSize/windowAspect,frontDist,backDist);
    projectionMatrix=glGetMatrix<double>(GL_PROJECTION_MATRIX);

    /* Reset the modelview matrix: */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    /* Set the light sources: */
    glEnable(GL_LIGHT0);
    glLight(GL_LIGHT0,headlight);

    /* Set the modelview matrix: */
    glTranslate(0.0,0.0,-dolly);
    glMultMatrix(modelView);
    modelViewMatrix=glGetMatrix<double>(GL_MODELVIEW_MATRIX);

    /* Lock the protein's state for rendering: */
    pthread_mutex_lock(&proteinMutex);

    for ( uint i = 0; i < numProteins(); ++i )
    {
        ProteinState *state = getProtein (i);
        if ( !state ) continue;
        glPushMatrix();
        state->interactor->glTransformProtein (contextData);

        /* Draw the protein: */
        state->proteinRenderer->glRenderAction (contextData);
        if ( state->selectedResidue )
        {
		    state->proteinRenderer->highlightResidue (
                contextData,
                state->selectedResidue
            );
        }
		if ( state->energyRenderer )
            state->energyRenderer->glRenderAction (contextData);
        glDisable (GL_LIGHTING);
		
		updateRamachandranWindow();
		updateEnergyVaule();

		if (state->proteinRenderer->getGrayOutCartoon())
			state->proteinRenderer->showRMSD(rmsd);

		/* Draw the structure selection box: */
        glPopMatrix();
        state->interactor->glRenderAction (contextData);
    }

    /* Release the protein lock: */
    pthread_mutex_unlock(&proteinMutex);
	focus(this);
}

int MyFlGlWindow::pick(int mouseX,int mouseY,int which)
{
    ProteinState *state = curProtein();
    if ( !state ) return 0;

    /* Activate this window's OpenGL context: */
    make_current();

    /* Calculate front- and backplane distances: */
    std::pair<double,double> depthRange=state->proteinRenderer->calcDepthRange(modelView);
    double frontDist=depthRange.first+dolly;
    if(frontDist<frontPlaneDistance)
        frontDist=frontPlaneDistance;
    double backDist=depthRange.second+dolly;
    if(backDist<frontDist+1.0)
        backDist=frontDist+1.0;

    /* Set the picking matrix: */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    GLdouble pickMatrix[4*4];
    for(int i=0;i<16;++i)
        pickMatrix[i]=0.0;
    double ws=double(w())/4.0;
    double hs=double(h())/4.0;
    pickMatrix[0]=ws;
    pickMatrix[5]=hs;
    pickMatrix[10]=1.0;
    pickMatrix[12]=(1.0-2.0*double(mouseX)/double(w()))*ws;
    pickMatrix[13]=(2.0*double(mouseY)/double(h())-1.0)*hs;
    pickMatrix[15]=1.0;
    glMultMatrixd(pickMatrix);

    /* Set the projection matrix: */
    double frontPlaneSize=frontDist*zoom;
    if(windowAspect>=1.0)
        glFrustum(-frontPlaneSize*windowAspect,frontPlaneSize*windowAspect,-frontPlaneSize,frontPlaneSize,frontDist,backDist);
    else
        glFrustum(-frontPlaneSize,frontPlaneSize,-frontPlaneSize/windowAspect,frontPlaneSize/windowAspect,frontDist,backDist);
    projectionMatrix=glGetMatrix<double>(GL_PROJECTION_MATRIX);

    /* Reset the modelview matrix: */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    /* Set the modelview matrix: */
    glTranslate(0.0,0.0,-dolly);
    glMultMatrix(modelView);

    /* Lock the protein's state for rendering: */
    pthread_mutex_lock(&proteinMutex);

    glPushMatrix();
    state->interactor->glTransformProtein(contextData);

    /* Draw the protein: */
    int structureIndex=state->proteinRenderer->glPick(contextData,which);

    glPopMatrix();

    /* Release the protein lock: */
    pthread_mutex_unlock(&proteinMutex);

    return structureIndex;
}

int MyFlGlWindow::handle(int eventType)
{
    ProteinState *state = curProtein();
    if ( !state )
        return Fl_Gl_Window::handle (eventType);

    int result=0;
    float plane[3], center[3];
    bool mustRedraw=false;
    const int modMask=FL_CAPS_LOCK|FL_SHIFT|FL_CTRL|FL_ALT;
    switch(eventType)
    {
        case FL_PUSH:
            isSpinning=false;
            dragX=Fl::event_x();
            dragY=Fl::event_y();
            if(dragButton==0&&Fl::event_button()==1)
            {
                if(numZoomingSteps>0)
                {
                    /* Try picking an atom: */
                    MD::Protein::ConstAtomIterator atom=state->protein->pickAtom(calcMouseRay());
                    if(atom!=state->protein->atomsEnd())
                    {
                        /* Move the point of interest to the picked atom: */
                        Point newViewCenter=atom->getPosition();
                        Point currentViewCenter=modelView.inverseTransform(Point::origin);
                        zoomingTranslation=(currentViewCenter-newViewCenter)/double(numZoomingSteps);
                        double newViewCenterDist=dolly-modelView.transform(newViewCenter)[2];
                        zoomingDollyStep=(newViewCenterDist*0.75-dolly)/double(numZoomingSteps);
                        Fl::add_timeout(0.01,zoomCB,this);
                    }
                    else
                        numZoomingSteps=0;
                }
                if((Fl::event_state()&modMask)==FL_CTRL)
                {
                    /* Try picking a secondary structure: */
                    // MD::Protein::StructureSelector picked=protein->pickStructure(calcMouseRay());
                    MD::Protein::StructureSelector picked = state->protein->pickStructure (
                        pick(Fl::event_x(), Fl::event_y(), 0)
                    );
                    state->interactor->selectStructure(picked);
                    updateGui();
                    mustRedraw=true;
                }
                if((Fl::event_state()&modMask)==(FL_CTRL|FL_SHIFT))
                {
                    /* Toggle a coil region: */
                    // MD::Protein::StructureSelector picked=protein->pickStructure(calcMouseRay());
                    MD::Protein::StructureSelector picked = state->protein->pickStructure (
                        pick(Fl::event_x(), Fl::event_y(), 0)
                    );
                    state->interactor->toggleCoil(picked);
                    updateGui();
                    mustRedraw=true;
                }
                if((Fl::event_state()&modMask)==FL_ALT)
                {
                    /* Try picking a residue: */
                    // MD::Protein::Residue* res=protein->pickResidue(calcMouseRay());
                    MD::Protein::Residue* res = state->protein->pickResidue (
                        pick(Fl::event_x(), Fl::event_y(), 1)
                    );
                    if(isZipping)
                    {
                        if(res!=0)
                        {
                            if(zipMode==ANTIPARALLEL)
                                zipAntiParallel(state->selectedResidue,res, plane, center);
                            else
                                zipParallel(state->selectedResidue,res, plane, center);
                            isZipping=false;
                        }
                    }
                    else
                        state->selectedResidue=res;
                    state->interactor->resetDragBox();
                    updateGui();
                    mustRedraw=true;
                }
                if((Fl::event_state()&modMask)==FL_SHIFT)
                {
                    if ( state->interactor->pickDragBox(
                            modelViewMatrix,
                            projectionMatrix,
                            DragBox::Point(2.0*double(dragX)/double(w())-1.0,
                            1.0-2.0*double(dragY)/double(h()),0.0)
                         ) )
                    {
                        isManipulating=true;

                        /* Initialize communication with the IK update thread: */
                        pthread_mutex_lock(&proteinMutex);
                        ikUpdateRequested=false;
                        ikUpdatePosted=false;
                        pthread_mutex_unlock(&proteinMutex);
                        Fl::add_timeout(0.01f,waitForIkCB,this);

                        mustRedraw=true;
                    }
                }
                else
                    dragButton=0x1;
            }
            else
                dragButton|=1<<(Fl::event_button()-1);
            result=1;
            break;

        case FL_RELEASE:
        {
            isSpinning=false;
            if(isManipulating)
            {
                isManipulating=false;

                /* Wait for the IK thread to finish: */
                pthread_mutex_lock(&proteinMutex);
                ikUpdateRequested=false;
                if(!ikUpdateDone)
                    pthread_cond_wait(&ikUpdateDoneCond,&proteinMutex);
                pthread_mutex_unlock(&proteinMutex);

                /* Snap the drag box back to its true position: */
                state->interactor->releaseDragBox();
                state->interactor->finishInteraction();

                updateGui();
                mustRedraw=true;
            }
            else if((Fl::event_state()&modMask)==0x0&&Fl::event_button()==1&&dragButton==1)
            {
                /* Button 1 has just been released. Check for spinning: */
                int dX=Fl::event_x()-dragX;
                int dY=dragY-Fl::event_y();

                if(dX!=0||dY!=0)
                {
                    /* Calculate the spin rotation: */
                    spinTransformation=calcRotation(dX,dY);
                    isSpinning=true;
                    Fl::add_timeout(0.01f,spinCB,this);
                }
            }
            // dragButton&=~(1<<(Fl::event_button()-1));
            dragButton=0; // Account for that Fltk bug!
            result=1;
            break;
        }

        case FL_DRAG:
        {
            isSpinning=false;
            int dX=Fl::event_x()-dragX;
            int dY=dragY-Fl::event_y();
            if(dX!=0||dY!=0)
            {
                if(isManipulating)
                {
                    /* Drag the box: */
                    pthread_mutex_lock(&proteinMutex);
                    state->interactor->dragBox (
                        modelViewMatrix,
                        projectionMatrix,
                        DragBox::Point(2.0*double(dragX)/double(w())-1.0,
                                       1.0-2.0*double(dragY)/double(h()),
                                       0.0)
                    );

                    /* Post an IK update request: */
                    ikUpdateRequested=true;
                    if(ikUpdateDone)
                        pthread_cond_signal(&ikUpdateRequestedCond);
                    pthread_mutex_unlock(&proteinMutex);
                }
                else if((Fl::event_state()&modMask)==0x0)
                {
                    if(dragButton==rotateButtonMask)
                        modelView.leftMultiply(calcRotation(dX,dY));
                    else if(dragButton==panButtonMask)
                        modelView.leftMultiply(MVTransformation::translate(Vector(dX*dolly*zoom/winRadius,dY*dolly*zoom/winRadius,0.0)));
                    else if(dragButton==zoomButtonMask)
                        dolly*=Math::exp(2.0*dY/winRadius);
                }
            }

            dragX=Fl::event_x();
            dragY=Fl::event_y();
            mustRedraw=true;
            result=1;
            break;
        }
        case FL_MOVE:
            result=1;
            break;
        case FL_ENTER:
            make_current();
            result=1;
            break;

        case FL_FOCUS:
            make_current();
            result=1;
            break;

        case FL_UNFOCUS:
            result=1;
            break;

        case FL_KEYBOARD:
        {
            bool eventTaken=false;
            MD::Protein::StructureSelector selectedStructure=state->interactor->getStructure();
            switch(Fl::event_key())
            {
                case 'a':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawAtoms (
                            selectedStructure,
                            !state->proteinRenderer->getDrawAtoms(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawAtoms (
                            !state->proteinRenderer->getDrawAtoms()
                        );
                        checkDrawToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 'b':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawBonds (
                            selectedStructure,
                            !state->proteinRenderer->getDrawBonds(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawBonds (
                            !state->proteinRenderer->getDrawBonds()
                        );
                        checkDrawToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 'c':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawHydrogenCages (
                            selectedStructure,
                            !state->proteinRenderer->getDrawHydrogenCages(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawHydrogenCages (
                            !state->proteinRenderer->getDrawHydrogenCages()
                        );
                        checkShowToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 'd':
                    if((Fl::event_state()&modMask)==0x0)
                    {
                        MD::BuildBeta();
                        eventTaken=true;
                    }
		    else
		        printf("Fl::event_state()&modMask %x\n", Fl::event_state()&modMask);
                    break;

                case 'e':
                    if((Fl::event_state()&modMask)==0x0)
                    {
                        if(state->energyCalculator!=0)
                            recalculateEnergy();
                        eventTaken=true;
                    }
                    break;

                case 'h':
                    if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawHydrogenBonds (
                            !state->proteinRenderer->getDrawHydrogenBonds()
                        );
                        checkShowToggles(state);
						eventTaken=true;
                    }
                    break;
                
                case 'k':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawCPK (
                            selectedStructure,
                            !state->proteinRenderer->getDrawCPK(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawCPK (
                            !state->proteinRenderer->getDrawCPK()
                        );
                        checkDrawToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 'l':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawLargeHydrogenCages (
                            selectedStructure,
                            !state->proteinRenderer->getDrawLargeHydrogenCages(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    break;

                case 'n':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawCartoon (
                            selectedStructure,
                            !state->proteinRenderer->getDrawCartoon(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawCartoon (
                            !state->proteinRenderer->getDrawCartoon()
                        );
                        checkDrawToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 'p':
                    if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawCollisions (
                            !state->proteinRenderer->getDrawCollisions()
                        );
                        checkShowToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 'r':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawBackboneRibbon (
                            selectedStructure,
                            !state->proteinRenderer->getDrawBackboneRibbon(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawBackboneRibbon (
                            !state->proteinRenderer->getDrawBackboneRibbon()
                        );
                        checkDrawToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 's':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawHydrogenBondSites (
                            selectedStructure,
                            !state->proteinRenderer->getDrawHydrogenBondSites(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawHydrogenBondSites (
                            !state->proteinRenderer->getDrawHydrogenBondSites()
                        );
                        checkShowToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 't':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawTube (
                            selectedStructure,
                            !state->proteinRenderer->getDrawTube(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawTube (
                            !state->proteinRenderer->getDrawTube()
                        );
                        checkDrawToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 'u':
                    if((Fl::event_state()&modMask)==0x0)
                    {
                        state->interactor->toggleUpdateDirection();
                        eventTaken=true;
                    }
                    break;

                case 'w':
                    if((Fl::event_state()&modMask)==FL_CTRL)
                    {
                        state->proteinRenderer->setDrawLine (
                            selectedStructure,
                            !state->proteinRenderer->getDrawLine(selectedStructure)
                        );
                        eventTaken=true;
                    }
                    else if((Fl::event_state()&modMask)==0x0)
                    {
                        state->proteinRenderer->setDrawLine (
                            !state->proteinRenderer->getDrawLine()
                        );
                        checkDrawToggles(state);
                        eventTaken=true;
                    }
                    break;

                case 'x':
                    if((Fl::event_state()&modMask)==0x0)
                    {
                        /* The currently selected residue and the next selected one will be zipped to form an anti-parallel bond: */
                        if(state->selectedResidue!=0)
                        {
                            isZipping=true;
                            zipMode=ANTIPARALLEL;
                        }
                        eventTaken=true;
                    }
                    break;

                case 'y':
                    if((Fl::event_state()&modMask)==0x0)
                    {
                        /* The currently selected residue and the next selected one will be zipped to form a parallel bond: */
                        if(state->selectedResidue!=0)
                        {
                            isZipping=true;
                            zipMode=PARALLEL;
                        }
                        eventTaken=true;
                    }
                    break;

                case 'z':
                    if((Fl::event_state()&modMask)==0x0)
                    {
                        isSpinning=false;
                        numZoomingSteps=50;
                        result=1;
                    }
                    break;
            }

            if(eventTaken)
            {
                updateGui();
                mustRedraw=true;
                result=1;
            }
            break;
        }

        default:
            result=Fl_Gl_Window::handle(eventType);
    }

    /* Redraw window contents only when necessary: */
    if(mustRedraw)
        redraw();

    return result;
}

MyFlGlWindow::Ray MyFlGlWindow::calcMouseRay(void) const
{
    double mouseX=2.0*double(dragX)/double(w())-1.0;
    double mouseY=1.0-2.0*double(dragY)/double(h());
    Point mouseStart(mouseX,mouseY,-1.0);
    Point mouseEnd(mouseX,mouseY,1.0);
    mouseStart=projectionMatrix.inverseTransform(mouseStart);
    mouseEnd=projectionMatrix.inverseTransform(mouseEnd);
    mouseStart=modelViewMatrix.inverseTransform(mouseStart);
    mouseEnd=modelViewMatrix.inverseTransform(mouseEnd);
    return Ray(mouseStart,mouseEnd-mouseStart);
}

MyFlGlWindow::MVTransformation MyFlGlWindow::calcRotation(int dX,int dY) const
{
    double offX=2.0*double(dragX)/double(w())-1.0;
    double offY=1.0-2.0*double(dragY)/double(h());
    double offZ=0.25;
    Geometry::Vector<double,3> v1(offX,offY,offZ);
    Geometry::Vector<double,3> v2(dX,dY,0.0);
    Geometry::Vector<double,3> axis=Geometry::cross(v1,v2);
    return MVTransformation::rotate(Rotation::rotateAxis(axis,1.5*Geometry::mag(v2)/winRadius));
}

void MyFlGlWindow::zoomCB(void* cbData)
{
    MyFlGlWindow* thisPtr=static_cast<MyFlGlWindow*>(cbData);
    if(thisPtr->numZoomingSteps>0)
    {
        --thisPtr->numZoomingSteps;
        if(thisPtr->numZoomingSteps>0)
            Fl::repeat_timeout(0.01f,zoomCB,thisPtr);
        thisPtr->modelView*=MVTransformation::translate(thisPtr->zoomingTranslation);
        thisPtr->dolly+=thisPtr->zoomingDollyStep;
        thisPtr->redraw();
    }
}

void MyFlGlWindow::spinCB(void* cbData)
{
    MyFlGlWindow* thisPtr=static_cast<MyFlGlWindow*>(cbData);
    if(thisPtr->isSpinning)
    {
        Fl::repeat_timeout(0.01f,spinCB,thisPtr);
        thisPtr->modelView.leftMultiply(thisPtr->spinTransformation);
        thisPtr->redraw();
    }
}

void MyFlGlWindow::spaceBallCB(void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    MyFlGlWindow* thisPtr=static_cast<MyFlGlWindow*>(cbData);
    MySpaceBall* spaceBall=thisPtr->spaceBall;

    bool mustRedraw=false;

    if(thisPtr->isManipulating&&ikUpdatePosted)
    {
        mustRedraw=true;
        ikUpdatePosted=false;
    }

    if(spaceBall->isButtonEventPending())
    {
        if(spaceBall->getEventButtonIndex()==spaceBall->getDragButtonIndex())
        {
            if(spaceBall->getDragState())
            {
                /* Button just went down; start dragging: */
                if(state->interactor->isValid()&&state->interactor->startInteraction())
                {
                    state->interactor->getBox().pickIncremental();
                    thisPtr->isManipulating=true;

                    /* Initialize communication with the IK update thread: */
                    pthread_mutex_lock(&proteinMutex);
                    ikUpdateRequested=false;
                    ikUpdatePosted=false;
                    pthread_mutex_unlock(&proteinMutex);

                    mustRedraw=true;
                }
            }
            else
            {
                /* Button just went up; stop dragging: */
                if(thisPtr->isManipulating)
                {
                    thisPtr->isManipulating=false;

                    /* Wait for the IK thread to finish: */
                    pthread_mutex_lock(&proteinMutex);
                    ikUpdateRequested=false;
                    if(!ikUpdateDone)
                        pthread_cond_wait(&ikUpdateDoneCond,&proteinMutex);
                    pthread_mutex_unlock(&proteinMutex);

                    /* Snap the drag box back to its true position: */
                    state->interactor->releaseDragBox();
                    state->interactor->finishInteraction();

                    updateGui();

                    mustRedraw=true;
                }
            }
        }
    }

    if(spaceBall->isAxisEventPending())
    {
        /* Get last translation/rotation vectors: */
        Vector translation,scaledRotationAxis;
        for(int i=0;i<3;++i)
        {
            translation[i]=spaceBall->getModifiedAxisState(i);
            scaledRotationAxis[i]=spaceBall->getModifiedAxisState(3+i);
        }

        if(thisPtr->isManipulating)
        {
            /* Transform motion into model coordinates: */
            translation=thisPtr->modelView.inverseTransform(translation);
            scaledRotationAxis=thisPtr->modelView.inverseTransform(scaledRotationAxis);
            translation*=thisPtr->dolly*thisPtr->zoom;

            /* Apply motion to drag box: */
            pthread_mutex_lock(&proteinMutex);
            state->interactor->getBox().dragIncremental (
                MVTransformation(translation,Rotation(scaledRotationAxis))
            );

            /* Post an IK update request: */
            ikUpdateRequested=true;
            if(ikUpdateDone)
                pthread_cond_signal(&ikUpdateRequestedCond);
            pthread_mutex_unlock(&proteinMutex);
        }
        else
        {
            /* Adjust linear motion to window size and zoom: */
            translation*=thisPtr->dolly*thisPtr->zoom;

            /* Switch between dollying and z-movement: */
            double dollyFactor=1.0;
            if(spaceBall->getZTranslationMode()==MySpaceBall::DOLLY)
            {
                dollyFactor=Math::exp(-translation[2]*spaceBall->getDollyGain());
                translation[2]=0.0;
            }

            /* Apply motion to current modelview transformation: */
            thisPtr->modelView.leftMultiply(MVTransformation(translation,Rotation(scaledRotationAxis)));
            thisPtr->dolly*=dollyFactor;
        }

        mustRedraw=true;
    }

    if(mustRedraw)
    {
        /* Redraw window: */
        redrawRenderWindows();
    }

    Fl::repeat_timeout(0.01f,spaceBallCB,thisPtr);
}

void MyFlGlWindow::waitForIkCB(void* cbData)
{
    MyFlGlWindow* thisPtr=static_cast<MyFlGlWindow*>(cbData);
    if(thisPtr->isManipulating)
    {
        Fl::repeat_timeout(0.1f,waitForIkCB,thisPtr);
        if(ikUpdatePosted)
        {
            ikUpdatePosted=false;
            redrawRenderWindows();
        }
    }
}

void MyFlGlWindow::zipCB(void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    MyFlGlWindow* thisPtr=static_cast<MyFlGlWindow*>(cbData);
    if(thisPtr->zipStep<thisPtr->numZipSteps)
    {
        Fl::repeat_timeout(0.01f,zipCB,thisPtr);
        ++thisPtr->zipStep;
        double s=double(thisPtr->zipStep)/double(thisPtr->numZipSteps);
        Rotation r=Rotation::rotateAxis(thisPtr->zipRotationAxis,thisPtr->zipRotationAngle*s);
        Point rc=Geometry::affineCombination(thisPtr->zipRotateCenterStart,thisPtr->zipRotateCenterEnd,s);
        Vector t=rc-r.transform(thisPtr->zipRotateCenterStart);
        MVTransformation stepTransformation(t,r);
        // MVTransformation stepTransformation(thisPtr->zipTranslation*s,Rotation::rotateAxis(thisPtr->zipRotationAxis,thisPtr->zipRotationAngle*s));
        state->interactor->drag(stepTransformation);
        state->interactor->applyChanges();
        state->proteinRenderer->updateProtein();
        if(state->energyCalculator!=0)
            state->energyCalculator->updateProtein();
        thisPtr->redraw();
        if(thisPtr->zipStep==thisPtr->numZipSteps)
        {
            state->interactor->finishInteraction();
            updateDihedralAngles();
        }
    }
}


void MyFlGlWindow::updateEnergyVaule(void)
{
    ProteinState *state = curProtein();
    if ( !state || !state->energyCalculator
		|| !state->proteinRenderer->getEnergyUpdate()) return;
	if (!settime) {
		time(&oldtime);
		settime = true;
		//energyValue = getEnergyOutput();
	}

	time_t newtime;
	time(&newtime);
	if (difftime(newtime, oldtime) > state->engUpdateRate)
	{
        energyValue = getEnergyOutput();
    	updateAtomEnergies();
		settime = false;
	}
	state->proteinRenderer->updateEnergy(energyValue);
}


void MyFlGlWindow::setRMSDValue(double value)
{
	rmsd = value;
}

MyFlGlWindow::MyFlGlWindow(int w,int h,const char* label)
    :Fl_Gl_Window(w,h,label),
     contextData(101),
     backgroundColor(retrieveValue(*configFile,"./RenderWindow/backgroundColor",GLColor<GLfloat,4>(0.0f,0.0f,0.0f))),
     headlight(GLLightSource::Color(1.0f,1.0f,1.0f),GLLightSource::Vector(0.0f,0.0f,1.0f,0.0f)),
     contextInitialized(false),projectionMatrixChanged(true),
     rotateButtonMask(retrieveValue(*configFile,"./RenderWindow/rotateButtonMask",0x1)),
     panButtonMask(retrieveValue(*configFile,"./RenderWindow/panButtonMask",0x2)),
     zoomButtonMask(retrieveValue(*configFile,"./RenderWindow/zoomButtonMask",0x3)),
     dragX(0),dragY(0),dragButton(0),
     spaceBall(0),settime(false),energyValue(0),rmsd(-1),
     numZoomingSteps(0),isSpinning(false),isManipulating(false),
     isZipping(false),zipStep(0),numZipSteps(0)
{
    /* Try creating a SpaceBall object: */
    try
    {
        std::string spaceBallSection=std::string("./RenderWindow/")+retrieveString(*configFile,"./RenderWindow/spaceBall");
        spaceBall=new MySpaceBall(configFile->getSectionIterator(spaceBallSection.c_str()));
        Fl::add_timeout(0.01f,spaceBallCB,this);
    }
    catch(...)
    {
        /* Ignore any errors, but disable SpaceBall: */
        spaceBall=0;
    }
}

MyFlGlWindow::MyFlGlWindow(int x,int y,int w,int h,const char* label)
    :Fl_Gl_Window(x,y,w,h,label),
     contextData(101),
     backgroundColor(retrieveValue(*configFile,"./RenderWindow/backgroundColor",GLColor<GLfloat,4>(0.0f,0.0f,0.0f))),
     headlight(GLLightSource::Color(1.0f,1.0f,1.0f),GLLightSource::Vector(0.0f,0.0f,1.0f,0.0f)),
     contextInitialized(false),projectionMatrixChanged(true),
     rotateButtonMask(retrieveValue(*configFile,"./RenderWindow/rotateButtonMask",0x1)),
     panButtonMask(retrieveValue(*configFile,"./RenderWindow/panButtonMask",0x2)),
     zoomButtonMask(retrieveValue(*configFile,"./RenderWindow/zoomButtonMask",0x3)),
     dragX(0),dragY(0),dragButton(0),
     spaceBall(0),settime(false),energyValue(0),rmsd(-1),
     numZoomingSteps(0),isSpinning(false),isManipulating(false),
     isZipping(false),zipStep(0),numZipSteps(0)
{
    /* Try creating a SpaceBall object: */
    try
    {
        std::string spaceBallSection=std::string("./RenderWindow/")+retrieveString(*configFile,"./RenderWindow/spaceBall");
        spaceBall=new MySpaceBall(configFile->getSectionIterator(spaceBallSection.c_str()));
        Fl::add_timeout(0.01f,spaceBallCB,this);
    }
    catch(...)
    {
        /* Ignore any errors, but disable SpaceBall: */
        spaceBall=0;
    }
}

MyFlGlWindow::~MyFlGlWindow(void)
{
    delete spaceBall;
}

void MyFlGlWindow::clearContext(ProteinState *state)
{
    if ( state->proteinRenderer )
        state->proteinRenderer->clearContext (contextData);
    if ( state->energyRenderer )
        state->energyRenderer->clearContext (contextData);
}

void MyFlGlWindow::changeRenderer(void)
{
    contextInitialized=false;
    invalidate();
}

void MyFlGlWindow::centerView(void)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    Point proteinCenter=state->protein->calcCentroid();
    double proteinRadius=state->protein->calcRadius();
    frontPlaneDistance=0.1;
    zoom=0.5;
    backPlaneDistance=30.0*proteinRadius;
    dolly=2.0*proteinRadius;
    modelView=MVTransformation::translate(Point::origin-proteinCenter);
    projectionMatrixChanged=true;
}

void MyFlGlWindow::pushViewSpec(void)
{
    /* Save current viewing specification in stack: */
    viewSpecStack.push_back(ViewSpec(zoom,dolly,modelView));
}

void MyFlGlWindow::popViewSpec(void)
{
    /* Retrieve topmost viewing specification from stack: */
    ViewSpec vs=viewSpecStack.front();
    viewSpecStack.pop_back();

    /* Set viewing specification: */
    zoom=vs.zoom;
    dolly=vs.dolly;
    modelView=vs.modelView;
}

void MyFlGlWindow::zipIt(const MyFlGlWindow::MVTransformation& zipTransformation)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    /* Decompose the given transformation: */
    zipTranslation=zipTransformation.getTranslation();
    zipRotationAxis=zipTransformation.getRotation().getAxis();
    zipRotationAngle=zipTransformation.getRotation().getAngle();
    zipRotateCenterStart=state->interactor->getBoxRotateCenter();
    zipRotateCenterEnd=zipTransformation.transform(zipRotateCenterStart);
    int numSteps=int(Math::ceil(Geometry::dist(zipRotateCenterStart,zipRotateCenterEnd)/2.0)); // At most 5 Angstrom per step
    int numRotateSteps=int(Math::ceil(Math::abs(zipRotationAngle)/0.2)); // At most 0.5 radians per step
    if(numSteps<numRotateSteps)
        numSteps=numRotateSteps;
    if(state->interactor->startInteraction())
    {
        zipStep=0;
        numZipSteps=numSteps;
        Fl::add_timeout(0.01f,zipCB,this);
    }
}

void MyFlGlWindow::redrawNow(void)
{
    make_current();
    draw();
    swap_buffers();
}

bool MyFlGlWindow::saveScreenshot(const char* fileName) const
{
    /* Force immediate redraw of window: */
    // redrawNow();

    /* Create a buffer for window contents: */
    GLColor<GLubyte,3>* buffer=new GLColor<GLubyte,3>[w()*h()];

    /* Set up OpenGL pixel transfer pipeline: */
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glPixelStorei(GL_PACK_SKIP_PIXELS,0);
    glPixelStorei(GL_PACK_ROW_LENGTH,0);
    glPixelStorei(GL_PACK_SKIP_ROWS,0);

    /* Retrieve window contents from OpenGL: */
    glReadPixels(0,0,w(),h(),GL_RGB,GL_UNSIGNED_BYTE,buffer);

    /* Open PPM file: */
    FILE* ppmFile=fopen(fileName,"wb");
    if(ppmFile==0) // Error during opening file?
    {
        /* Return: */
        delete[] buffer;
        return false;
    }

    /* Write binary PPM file header: */
    fprintf(ppmFile,"P6\n");
    fprintf(ppmFile,"# Created by ProtoShop\n");
    fprintf(ppmFile,"%d %d\n",w(),h());
    fprintf(ppmFile,"255\n");

    /* Write buffer to a PPM file one line at a time to reverse orientation: */
    for(int y=h()-1;y>=0;--y)
    {
        /* Write one row worth of pixel values: */
        fwrite(buffer+y*w(),sizeof(GLColor<GLubyte,3>),w(),ppmFile);
    }

    /* Clean up and return: */
    fclose(ppmFile);
    delete[] buffer;
    return true;
}
