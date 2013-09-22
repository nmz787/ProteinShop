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

#ifndef RAMACHANDRANPLOT_INCLUDED
#define RAMACHANDRANPLOT_INCLUDED

#include <FL/gl.h>
#include <FL/Fl_Gl_Window.H>
#include <GLTypes.h>

#include "MDGeometry.h"

/* Forward declarations: */
namespace MD {
class Protein;
}
class ProteinInteractor;

class RamachandranPlot:public Fl_Gl_Window
{
    /* Elements: */
    private:
    GLColor<GLfloat,4> backgroundColor; // Background color for the GL window
    GLColor<GLfloat,4> selectedStructureColor; // Rendering color for the selected structure
    GLColor<GLfloat,4> activeCoilRegionColor; // Rendering color for active coil regions
    MD::Protein* protein; // Pointer to protein whose dihedral angles are being visualized
    ProteinInteractor* proteinInteractor; // Pointer to protein interactor containing selected structure and active coil regions
    bool structureValid; // Flag if the interactor currently contains a selected secondary structure
    int numStructureResidues; // Number of dihedral angles in currently selected structure
    MD::Scalar* structurePhis; // Array of phi dihedral angles for selected structure
    MD::Scalar* structurePsis; // Array of psi dihedral angles for selected structure
    bool activeCoilsValid; // Flag if the interactor currently contains active coil regions
    int numActiveCoilsResidues; // Total number of dihedral angles in currently active coil regions
    MD::Scalar* activeCoilsPhis; // Array of phi dihedral angles for all active coil regions
    MD::Scalar* activeCoilsPsis; // Array of psi dihedral angles for all active coil regions
    int ramachandranTextureWidth,ramachandranTextureHeight;
    GLubyte* ramachandranTexture; // Color texture containing Ramachandran probabilities
	bool activeResidueValid; // Flag if the interactor currently contains active coil regions
	double phiAngle, psiAngle;
	double oldPhiAngle, oldPsiAngle;
	int dragX,dragY,dragButton; // State to keep track of mouse movements inside the window
	bool boundary;
    MD::Scalar* phis;
    MD::Scalar* psis;
    
    /* Methods inherited from Fl_Gl_Window: */
    protected:
    virtual void draw(void);
    int handle(int eventType);
    
    /* New protected methods: */
    void loadRamachandranTexture(void);
    
    /* Constructors and destructors: */
    public:
    RamachandranPlot(int w,int h,const char* label =0);
    RamachandranPlot(int x,int y,int w,int h,const char* label =0);
    ~RamachandranPlot(void);
    
    /* Methods: */
    void setProtein(MD::Protein* newProtein); // Sets a new protein to visualize
    void setProteinInteractor(ProteinInteractor* newProteinInteractor); // Sets a protein interactor to visualize dihedral angles
    void updateProtein(void); // Notifies the Ramachandran plot window that the protein has been changed
    void redrawNow(void); // Redraws the Ramachandran plot immediately
    void updateSelectedResidueAngles(double phi, double psi); // Redraws the Ramachandran plot immediately
    void updateSelectedResidueAngles(void); // Redraws the Ramachandran plot immediately
    void resetSelectedResidueAngles(); // Redraws the Ramachandran plot immediately
//	void updateWindow( ) const;
};

#endif
