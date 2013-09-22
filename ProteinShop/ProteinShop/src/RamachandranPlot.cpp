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
RamachandranPlot - Class for Fltk windows displaying Ramachandran plots.
***********************************************************************/

#include <Math/Constants.h>
#include <GLTemplates.h>
#include <FL/Fl.H>

#include "ProteinInteractor.h"
#include "RamachandranPlot.h"
#include "GLText.h"
#include "Globals.h"

extern void updateVisualization (bool atomPositionsChanged = true);

/*********************************
Methods of class RamachandranPlot:
*********************************/

void RamachandranPlot::draw(void)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
	 
	 GLText text;
    if(!valid())
    {
        /* Set the viewport: */
        glViewport(0,0,w(),h());
        
        /* Set the projection and modelview matrices: */
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        float factor=1.0/Math::Constants<double>::pi;
        glScale(factor,factor,factor);
        glDisable(GL_DEPTH_TEST);
    }
    
    /* Clear window: */
    glClearColor(backgroundColor);
    glClear(GL_COLOR_BUFFER_BIT);
    
    /* Draw ramachandran texture background: */
    glEnable(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_BASE_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAX_LEVEL,0);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS,0);
    glPixelStorei(GL_UNPACK_ROW_LENGTH,0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS,0);
    glPixelStorei(GL_UNPACK_IMAGE_HEIGHT,0);
    glPixelStorei(GL_UNPACK_SKIP_IMAGES,0);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGB8,ramachandranTextureHeight,ramachandranTextureWidth,0,GL_RGB,GL_UNSIGNED_BYTE,ramachandranTexture);
    glColor3f(1.0f,0.0f,1.0f);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f,0.0f);
    glVertex2f(-Math::Constants<float>::pi,-Math::Constants<float>::pi);
    glTexCoord2f(1.0f,0.0f);
    glVertex2f(Math::Constants<float>::pi,-Math::Constants<float>::pi);
    glTexCoord2f(1.0f,1.0f);
    glVertex2f(Math::Constants<float>::pi,Math::Constants<float>::pi);
    glTexCoord2f(0.0f,1.0f);
    glVertex2f(-Math::Constants<float>::pi,Math::Constants<float>::pi);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    
    /* Draw coordinate axes: */
    glColor3f(1.0f-backgroundColor[0],1.0f-backgroundColor[1],1.0f-backgroundColor[2]);
    glBegin(GL_LINES);
    glVertex2f(0.0f,-4.0f);
    glVertex2f(0.0f,4.0f);
    glVertex2f(-4.0f,0.0f);
    glVertex2f(4.0f,0.0f);
    glEnd();

	/* Draw axes data range */
	text.setFontColor(0.0, 0.0, 0.0);
	/* h, v */
	/*	+  */
	/*-   +*/
	/*	-  */	
	text.drawFreeString("pi",2.3f,-2.7f);
	text.drawFreeString("-pi",-2.8f,-2.7f);
	text.drawFreeString("pi",-3.0f,2.6f);
	text.drawFreeString("0",-3.0f,-0.3f);
	text.drawFreeString("0",-0.3f,-2.7f);
	text.drawFreeString("Phi",0.15f,2.5f);
	text.drawFreeString("Psi",2.5f,0.15f);

	/* Draw axes label */
	text.setFontScale(0.003, 0.003);
	text.setFontColor(0.0, 0.0, 1.0);
    
    if(structureValid)
    {
        /* Visualize dihedral angles inside a selected structure: */
        glPointSize(3.0f);
        glBegin(GL_POINTS);
        glColor(selectedStructureColor);
        for(int i=0;i<numStructureResidues;++i)
            glVertex(structurePhis[i],structurePsis[i]);
        glEnd();
    }
    
    /* Visualize dihedral angles inside all active coil regions: */
    if(activeCoilsValid)
    {
        /* Visualize dihedral angles inside all active coil regions: */
        glPointSize(3.0f);
        glBegin(GL_POINTS);
        glColor(activeCoilRegionColor);
        for(int i=0;i<numActiveCoilsResidues;++i)
            glVertex(activeCoilsPhis[i],activeCoilsPsis[i]);
        glEnd();
    }

	if ( state->selectedResidue != 0 )
	{
    	/* Highlight dihedral angles of selected residue: */
		glEnable(GL_BLEND);
		glEnable(GL_POINT_SMOOTH);
        glPointSize(6.0f);
		glBegin(GL_POINTS);
        
		if(boundary)
			glColor3f(0.0f, 0.0f, 5.0f);
		else
			glColor3f(1.0f, 1.0f, 1.0f);
			//printf("draw phi %f psi %f\n", oldPhiAngle,oldPsiAngle);
			glVertex2f(oldPhiAngle,oldPsiAngle);
        
		glEnd();
		}
}

void RamachandranPlot::loadRamachandranTexture(void)
{
    /* Open texture PPM file: */
    FILE* ppmFile=fopen("Ramachandran64.ppm","rb");
    
    /* Read image file header: */
    char line[80];
    fgets(line,80,ppmFile); // Skip "P6"
    fgets(line,80,ppmFile); // Skip creator comment
    fgets(line,80,ppmFile); // Read width and height
    sscanf(line,"%d %d",&ramachandranTextureWidth,&ramachandranTextureHeight);
    fgets(line,80,ppmFile); // Skip max value line
    
    /* Read texture image: */
    ramachandranTexture=new GLubyte[ramachandranTextureWidth*ramachandranTextureHeight*3];
    for(int y=ramachandranTextureHeight-1;y>=0;--y)
        fread(ramachandranTexture+y*ramachandranTextureWidth*3,sizeof(GLubyte)*3,ramachandranTextureWidth,ppmFile);
    fclose(ppmFile);
}

RamachandranPlot::RamachandranPlot(int w,int h,const char* label)
    :Fl_Gl_Window(w,h,label),
     backgroundColor(retrieveValue(*configFile,"./RamachandranPlot/backgroundColor",GLColor<GLfloat,4>(0.0f,0.0f,0.0f))),
     selectedStructureColor(retrieveValue(*configFile,"./RamachandranPlot/selectedStructureColor",GLColor<GLfloat,4>(0.0f,1.0f,0.0f))),
     activeCoilRegionColor(retrieveValue(*configFile,"./RamachandranPlot/activeCoilRegionColor",GLColor<GLfloat,4>(1.0f,1.0f,0.0f))),
     protein(0),proteinInteractor(0),
     structureValid(false),numStructureResidues(0),structurePhis(0),structurePsis(0),
     activeCoilsValid(false),numActiveCoilsResidues(0),activeCoilsPhis(0),activeCoilsPsis(0),
     ramachandranTextureWidth(0),ramachandranTextureHeight(0),ramachandranTexture(0)
{
    loadRamachandranTexture();
}

RamachandranPlot::RamachandranPlot(int x,int y,int w,int h,const char* label)
    :Fl_Gl_Window(x,y,w,h,label),
     backgroundColor(retrieveValue(*configFile,"./RamachandranPlot/backgroundColor",GLColor<GLfloat,4>(0.0f,0.0f,0.0f))),
     selectedStructureColor(retrieveValue(*configFile,"./RamachandranPlot/selectedStructureColor",GLColor<GLfloat,4>(0.0f,1.0f,0.0f))),
     activeCoilRegionColor(retrieveValue(*configFile,"./RamachandranPlot/activeCoilRegionColor",GLColor<GLfloat,4>(1.0f,1.0f,0.0f))),
     protein(0),proteinInteractor(0),
     structureValid(false),numStructureResidues(0),structurePhis(0),structurePsis(0),
     activeCoilsValid(false),numActiveCoilsResidues(0),activeCoilsPhis(0),activeCoilsPsis(0),
     ramachandranTextureWidth(0),ramachandranTextureHeight(0),ramachandranTexture(0)
{
    loadRamachandranTexture();
}

RamachandranPlot::~RamachandranPlot(void)
{
    delete[] structurePhis;
    delete[] structurePsis;
    delete[] activeCoilsPhis;
    delete[] activeCoilsPsis;
    delete[] ramachandranTexture;
}

void RamachandranPlot::setProtein(MD::Protein* newProtein)
{
    protein=newProtein;
    proteinInteractor=0;
}

void RamachandranPlot::setProteinInteractor(ProteinInteractor* newProteinInteractor)
{
    proteinInteractor=newProteinInteractor;
}

void RamachandranPlot::updateProtein(void)
{
    /* Retrieve dihedral angles from selected structure, if any: */
    if(protein!=0&&proteinInteractor!=0)
    {
        MD::Protein::StructureSelector selectedStructure=proteinInteractor->getStructure();
        structureValid=selectedStructure.isValid();
        if(structureValid)
        {
            /* Retrieve dihedral angles from the selected structure: */
            if(numStructureResidues!=selectedStructure.getNumResidues())
            {
                /* Reallocate dihedral angle arrays: */
                delete[] structurePhis;
                delete[] structurePsis;
                numStructureResidues=selectedStructure.getNumResidues();
                structurePhis=new MD::Scalar[numStructureResidues];
                structurePsis=new MD::Scalar[numStructureResidues];
            }
            selectedStructure.getDihedralAngles(structurePhis,structurePsis);
        }
    }
    else
        structureValid=false;
    
    /* Retrieve dihedral angles from active coil regions: */
    int newNumActiveCoilsResidues=0;
    for(std::list<MD::Protein::StructureSelector>::const_iterator cIt=proteinInteractor->getLeftCoils().begin();cIt!=proteinInteractor->getLeftCoils().end();++cIt)
        newNumActiveCoilsResidues+=cIt->getNumResidues();
    for(std::list<MD::Protein::StructureSelector>::const_iterator cIt=proteinInteractor->getRightCoils().begin();cIt!=proteinInteractor->getRightCoils().end();++cIt)
        newNumActiveCoilsResidues+=cIt->getNumResidues();
    if(newNumActiveCoilsResidues!=0)
    {
        /* Retrieve dihedral angles from each active coil region: */
        if(numActiveCoilsResidues!=newNumActiveCoilsResidues)
        {
            /* Reallocate dihedral angle arrays: */
            delete[] activeCoilsPhis;
            delete[] activeCoilsPsis;
            numActiveCoilsResidues=newNumActiveCoilsResidues;
            activeCoilsPhis=new MD::Scalar[numActiveCoilsResidues];
            activeCoilsPsis=new MD::Scalar[numActiveCoilsResidues];
        }
        int baseIndex=0;
        for(std::list<MD::Protein::StructureSelector>::const_iterator cIt=proteinInteractor->getLeftCoils().begin();cIt!=proteinInteractor->getLeftCoils().end();++cIt)
        {
            cIt->getDihedralAngles(activeCoilsPhis+baseIndex,activeCoilsPsis+baseIndex);
            baseIndex+=cIt->getNumResidues();
        }
        for(std::list<MD::Protein::StructureSelector>::const_iterator cIt=proteinInteractor->getRightCoils().begin();cIt!=proteinInteractor->getRightCoils().end();++cIt)
        {
            cIt->getDihedralAngles(activeCoilsPhis+baseIndex,activeCoilsPsis+baseIndex);
            baseIndex+=cIt->getNumResidues();
        }
        activeCoilsValid=true;
    }
    else
        activeCoilsValid=false;
}
void RamachandranPlot::resetSelectedResidueAngles()
{
    ProteinState *state = curProtein();
    if ( !state ) return;
#ifdef COIL
	if ( state->selectedResidue != 0 &&
         state->selectedResidue->getType() != MD::Protein::Residue::UNK &&
         state->selectedResidue->getSecondaryStructure()->getStructureType() ==
            MD::Protein::SecondaryStructure::COIL )
#endif
	if ( state->selectedResidue != 0 &&
         state->selectedResidue->getType() != MD::Protein::Residue::UNK  )
    {
        /* reset  dihedral angles of selected residue: */
        MD::Scalar newPhi = MD::Scalar( phiAngle );
        MD::Scalar newPsi = MD::Scalar( psiAngle );

        /* Apply the change: */
        undoBuffer.startInteraction (
            state->protein,
            state->selectedResidue,
            state->interactor->getUpdateDirection()
        );
        state->protein->setDihedralAngles (
            state->protein->getResidueIndex(state->selectedResidue),
            1,&newPhi,&newPsi,state->interactor->getUpdateDirection()
        );

        undoBuffer.finishInteraction();
        /* Update visualization state: */
		updateVisualization();
	}
}


void RamachandranPlot::updateSelectedResidueAngles( double phi, double psi)
{
	/* update dihedral angles for Ramachandran Plot: */
	phiAngle = phi;
	psiAngle = psi;
	oldPhiAngle=phi; 
	oldPsiAngle=psi; 
	redraw();
}

int RamachandranPlot::handle(int eventType)
{
    int result=0;
    //ProteinState *state = curProtein();
    //if ( !state ) return result;
    float deltaX=0.0f;
    float deltaY=0.0f;

    bool mustRedraw=false;
    const int modMask=FL_ALT;
    switch(eventType)
    {
        case FL_PUSH:
            make_current();

			dragX=Fl::event_x();
            dragY=Fl::event_y();
            if(dragButton==0&&Fl::event_button()==1)
            {
                if((Fl::event_state()&modMask)==FL_ALT)
                {
                    //pick(Fl::event_x(), Fl::event_y());
                    //mustRedraw=true;
                }
            }
            else
                dragButton|=1<<(Fl::event_button()-1);
            	result=1;
            break;

        case FL_RELEASE:
        {
            // dragButton&=~(1<<(Fl::event_button()-1));
            dragButton=0; // Account for that Fltk bug!
            result=1;
            break;
        }
        case FL_ENTER:
        {
			boundary = true;
			mustRedraw=true;
            dragButton=0; // Account for that Fltk bug!
            result=1;
            break;
        }
        case FL_LEAVE:
        {
			boundary = false;
			mustRedraw=true;
            dragButton=0; // Account for that Fltk bug!
            result=1;
            break;
        }
        case FL_DRAG:
        {
            if((float)Fl::event_x()<0 || (float)Fl::event_y()<0 )
			break;
			
			deltaX=((float)Fl::event_x()-dragX)/10;
            deltaY=((float)dragY-Fl::event_y())/10;
			if(deltaX!=0||deltaY!=0)
            {
		    	phiAngle+=deltaX*2.0f*Math::Constants<double>::pi/w(); 
				psiAngle+=deltaY*2.0f*Math::Constants<double>::pi/h();
 				/* Update the dihedral angles: */
				resetSelectedResidueAngles();
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

        case FL_FOCUS:
            result=1;
            break;

        case FL_UNFOCUS:
            result=1;
            break;

        default:
            result=Fl_Gl_Window::handle(eventType);
    }

    /* Redraw window contents only when necessary: */
    if(mustRedraw)
        redraw();

    return result;
}

void RamachandranPlot::redrawNow(void)
{
    make_current();
    draw();
    //swap_buffers();
}
