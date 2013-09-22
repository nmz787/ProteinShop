/*****************************************************************************
Copyright (c) 2004, The Regents of the University of California, 
through Lawrence Berkeley National Laboratory, Univ. of Calif. at Davis, 
and Lawrence Livermore National Laboratory 
(subject to receipt of any required approvals from U.S. Dept. of Energy).  
All rights reserved.

This program is free software; you can redistribute it and/or 
modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation; either version 2 
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License 
along with this program; if not, write to the Free Software 
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 
******************************************************************************/ 

/***********************************************************************
ProteinFltk - GUI callbacks, thread callbacks, and entry point for
ProteinShop (main() at bottom of file).

***********************************************************************/

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <Geometry/MathTemplates.h>
#include <Geometry/MathConstants.h>
#include <Geometry/Vector.h>

#include <FL/fl_ask.H>
#ifdef FLTK_NEW
#include <FL/Fl_File_Chooser.H>
#else
#include <FL/fl_file_chooser.H>
#endif

#include "ColorFunction.h"
#include "CreateProtein.h"
#include "fitsq.h"
#include "ProteinGui.h"
#include "ParsePdbFile.h"
#include "Globals.h"
#include "Stride.h"
#include "ProteinFactory.h"

extern int offlineBuildBeta;
extern float strand_to_coil_penalty, flatfrac, flatPhi, flflatPhi, raflatPhi,
   flatPsi, flflatPsi, raflatPsi;

static int* energyDialogComponentToggleCBDatas=0;
bool  validPdbForEnergy = false;
int alignBase = 2;
void energyDialogComponentToggleCB(Fl_Check_Button* button,void* cbData);
const size_t inputfilename_size = 200;
char inputfilename[inputfilename_size];
Fl_File_Chooser	*pdbfc, *predfc;
Fl_Input *pdbFileInput, *predFileInput; 
Fl_Value_Input *copyBeginIndexInput, *copyEndIndexInput;

// variables used to control the energy optimizer
static const char runLabel[] = "Run Local Optimization";
static const char stopLabel[] = "Stop Local Optimization";
static ProteinState *optState = 0;
static uint updateNumber, updateInterval;

// buffer used for snprintf() calls
static const uint buflen = 256;
static char buf[buflen];

// variables for the volume dialog UI
static uint volumeDialogClassNumber = 0;
static char volumeDialogClassifierGroupLabel[buflen];
static char volumeDialogColorChoiceLabel[buflen];

// variables for the record dialog UI
static ProteinRecord *proteinRecord = 0;
static ProteinRecord *backupRecord = 0;

typedef ProteinInteractor::Transformation Transformation;
typedef Transformation::Vector Vector;

/*****************
Utility functions:
*****************/

float testCross(Vector d1, Vector d2, Vector c) {
    int i;
    float den;
    float s1 = 0., s2 = 0., sc = 0.;
    for (i = 0; i < 3; ++i)  {
       s1 = s1 + d1[i]*d1[i];
       s2 = s2 + d2[i]*d2[i];
       sc = sc + c[i]*c[i];
       }
    den = s1*s2;
    if(den == 0.) return 0.;
    return sc/den;
    }

void updateAtomEnergies(void)
{
    ProteinState *state = curProtein();
    if ( !state || !state->energyCalculator || !state->visualizeEnergy ) return;

    /* Recalculate atom energy values according to current visualization switches: */
    for ( uint i = 0; i < state->protein->getNumAtoms(); ++i )
        state->proteinRenderer->setAtomValue (i, state->energyCalculator->getAtomEnergy(i));
}

void recalculateEnergy(void)
{
    ProteinState *state = curProtein();
    if ( !state || !state->energyCalculator ) return;

    energyDialogEnergyValueOutput->value (state->energyCalculator->calcEnergy());
    updateAtomEnergies();

    renderWindow->redraw();
    updateGui();
}

void updateEnergyOutput(double value)
{
    if ( !isnan(value) )
        energyDialogEnergyValueOutput->value(value);
}

double getEnergyOutput()
{
    return energyDialogEnergyValueOutput->value();
}

void initializeVolumeDialogChoices()
{
    // load color choice with menu items
    volumeDialogColorChoice->clear();
    for ( uint i = 0; i < ColorFunction::numFunctions(); ++i )
    {
        ColorFunction *function = ColorFunction::get (i);
        volumeDialogColorChoice->add (function->name());
    }
    // load classifier choice with menu items
    volumeDialogClassifierChoice->clear();
    for ( uint i = 0; i < AtomClassifier::numClassifiers(); ++i )
    {
        AtomClassifier *classifier = AtomClassifier::get (i);
        volumeDialogClassifierChoice->add (classifier->name());
    }
}

void loadEnergyLibrary(const char* energyLibraryName)
{
    /* Delete current energy library: */
    delete energyLibrary;
    energyLibrary=0;

    /* Remove all energy component toggles from energy component group: */
    int dialogHeight=energyDialog->h();
    int groupHeight=energyDialogComponentGroup->h();
    for(int i=energyDialogComponentGroup->children()-1;i>=0;--i)
    {
        Fl_Widget* toggle=energyDialogComponentGroup->child(i);
        dialogHeight-=toggle->h();
        groupHeight-=toggle->h();
        energyDialogComponentGroup->remove(toggle);
        delete toggle;
    }
    energyDialogComponentGroup->size(energyDialogComponentGroup->w(),groupHeight);
    energyDialog->size(energyDialog->w(),dialogHeight);
    delete[] energyDialogComponentToggleCBDatas;

    try
    {
        /* Load new energy library: */
        energyLibrary=EnergyLibrary::load(energyLibraryName);
    }
    catch(std::runtime_error error)
    {
        /* Disable the energy calculation module: */
        energyLibrary=0;

        /* Print error message: */
        msg.Error(ENGID,"",error.what());

        return;
    }

    /* Create toggle buttons for all energy components in energy component group: */
    energyDialogComponentToggleCBDatas=new int[energyLibrary->getNumEnergyComponents()];
    groupHeight+=energyLibrary->getNumEnergyComponents()*25;
    energyDialogComponentGroup->size(energyDialogComponentGroup->w(),groupHeight);
    dialogHeight+=energyLibrary->getNumEnergyComponents()*25;
    energyDialog->size(energyDialog->w(),dialogHeight);

    // energyDialogComponentGroup->begin();
    Fl_Group::current(0);
    int x=energyDialogComponentGroup->x()+5;
    int y=energyDialogComponentGroup->y();
    int w=energyDialogComponentGroup->w()-10;
    for(int i=0;i<energyLibrary->getNumEnergyComponents();++i)
    {
        y+=25;
        Fl_Check_Button* newToggle=new Fl_Check_Button(x,y,w,25,energyLibrary->getEnergyComponentName(i));
        newToggle->down_box(FL_DIAMOND_DOWN_BOX);
        energyDialogComponentToggleCBDatas[i]=i;
        newToggle->callback((Fl_Callback*)energyDialogComponentToggleCB,&energyDialogComponentToggleCBDatas[i]);
        energyDialogComponentGroup->add(newToggle);
    }
    // energyDialogComponentGroup->end();

    initializeVolumeDialogChoices();
}


void hideRecordDialog()
{
    recordDialog->hide();
    delete proteinRecord;
    proteinRecord = 0;
    ProteinState *state = curProtein();
    if ( state && backupRecord )
        backupRecord->apply (*state);
    delete backupRecord;
    backupRecord = 0;
    energyDialogLoadRecordButton->activate();
}


void selectProtein (ProteinState *state)
{
    // called to update the workspace when the current protein changes
    if ( !state ) return;
    setCurProtein (state);
    sequenceDialogSequenceWidget->setProtein (state->protein);
    sequenceDialogSequenceWidget->setProteinRenderer (state->proteinRenderer);
    ramachandranPlot->setProtein (state->protein);
    ramachandranPlot->setProteinInteractor (state->interactor);
    state->proteinRenderer->updateProtein();
    state->interactor->resetDragBox();
    renderWindow->changeRenderer();
    renderWindow->redraw();
    if ( recordDialog->visible_r() )
        hideRecordDialog();
    updateGui();
}


void hideEnergyDialog()
{
    if ( optState && !strcmp(energyDialogOptimizeButton->label(), stopLabel) )
    {
        EnergyOptimizer *optimizer = optState->energyCalculator->optimizer();
        if ( optimizer )
        {
            optimizer->cancel();
            energyDialogOptimizeButton->label (runLabel);
        }
        optState = 0;
    }
    energyDialog->hide();
    volumeDialog->hide();
    hideRecordDialog();
}


void addProtein (MD::Protein* newProtein, int Id)
{
    // called to add a new protein to the workspace
    using namespace MD;

    if ( !newProtein ) return;
    ProteinState *state = createProtein (inputfilename);
    if ( !state ) return;

    state->protein = newProtein;
    configFile->setCurrentSection ("/ProteinRenderer");
    state->proteinRenderer = new ProteinRenderer (
        configFile->getSectionIterator("/ProteinRenderer"),
        state->protein
    );
    configFile->setCurrentSection("/");
    state->proteinId = Id;
//	printf("pro id %d\n",Id);
	if ( validPdbForEnergy )
    {
        showEnergyDialogToggle->activate();
        msg.Debug (PDBID, "This PDB file has the correct number of atoms.");
        if ( energyLibrary )
        {
            state->energyCalculator = energyLibrary->createEnergyCalculator (state->protein);
            if ( state->energyCalculator )
            {
                state->visualizeEnergy = energyDialogVisualizeEnergyToggle->value();
                state->visualizeEnergyMinRange = energyDialogRangeMinInput->value();
                state->visualizeEnergyMaxRange = energyDialogRangeMaxInput->value();
                state->energyRenderer = new EnergyRenderer (state);
            }
            else
            {
                showEnergyDialogToggle->deactivate();
                msg.Warn (PDBID, "There is no Energy Calculator.");
            }
        }
    }
    else
    {
        // some hydrogen atoms may be missing from this molecule
        showEnergyDialogToggle->deactivate();
        hideEnergyDialog();
    }
    state->proteinRenderer->setMapAtomValues (state->visualizeEnergy);
    if ( state->visualizeEnergy )
    {
        state->proteinRenderer->setMapAtomValueRange (
            state->visualizeEnergyMinRange,
            state->visualizeEnergyMaxRange
        );
    }
    state->interactor = new ProteinInteractor (
        state->protein,
        state->proteinRenderer,
        state->energyCalculator,
        &undoBuffer
    );
    // update the user interface
    selectionDialogBrowser->add (state->name, state);
    selectionDialogBrowser->deselect();
    selectionDialogBrowser->select (selectionDialogBrowser->size());
    selectionDialogBrowser->bottomline (selectionDialogBrowser->size());
    selectionDialogAlignButton->deactivate();
    selectProtein (state);
    updateAtomEnergies();
    renderWindow->centerView();
}


//void checkDrawToggles(void)
void checkDrawToggles(ProteinState *state)
{
    /* Need to draw backbone? */
    bool dontDrawBackbone=state->proteinRenderer->getDrawAtoms();
    dontDrawBackbone=dontDrawBackbone||state->proteinRenderer->getDrawBonds();
    dontDrawBackbone=dontDrawBackbone||state->proteinRenderer->getDrawCartoon();
    dontDrawBackbone=dontDrawBackbone||state->proteinRenderer->getDrawCPK();
    dontDrawBackbone=dontDrawBackbone||state->proteinRenderer->getDrawTube();
    dontDrawBackbone=dontDrawBackbone||state->proteinRenderer->getDrawLine();
    state->proteinRenderer->setDrawBackbone(!dontDrawBackbone);
	
	// update all proteins with the same id
	for ( uint i = 0; i < numProteins(); ++i )
    {
        ProteinState *s2 = getProtein (i);
        if ( !s2 || s2 == state) continue;
		if(state->proteinId == s2->proteinId)
    	{
			s2->proteinRenderer->setDrawBonds(state->proteinRenderer->getDrawBonds());
			s2->proteinRenderer->setDrawCartoon(state->proteinRenderer->getDrawCartoon());
    		s2->proteinRenderer->setDrawCPK(state->proteinRenderer->getDrawCPK());
    		s2->proteinRenderer->setDrawTube(state->proteinRenderer->getDrawTube());
    		s2->proteinRenderer->setDrawLine(state->proteinRenderer->getDrawLine());
    		s2->proteinRenderer->setDrawAtoms(state->proteinRenderer->getDrawAtoms());
    		s2->proteinRenderer->setDrawBackboneRibbon(state->proteinRenderer->getDrawBackboneRibbon());
			s2->proteinRenderer->setDrawBackbone(!dontDrawBackbone);
		}
	}
}
void checkShowToggles(ProteinState *state)
{
	// update all proteins with the same id
	for ( uint i = 0; i < numProteins(); ++i )
    {
        ProteinState *s2 = getProtein (i);
        if ( !s2 || s2 == state) continue;
		if(state->proteinId == s2->proteinId)
    	{
    		s2->proteinRenderer->setDrawHydrogenBonds(state->proteinRenderer->getDrawHydrogenBonds());
    		s2->proteinRenderer->setDrawHydrogenBondSites(state->proteinRenderer->getDrawHydrogenBondSites());
    		s2->proteinRenderer->setDrawHydrogenCages(state->proteinRenderer->getDrawHydrogenCages());
    		s2->proteinRenderer->setDrawCollisions(state->proteinRenderer->getDrawCollisions());
    		s2->proteinRenderer->setShowHydrophobic(state->proteinRenderer->getShowHydrophobic());
    		s2->proteinRenderer->setShowHydrophilic(state->proteinRenderer->getShowHydrophilic());
   			s2->proteinRenderer->setShowDisulfide(state->proteinRenderer->getShowDisulfide());
		}
	}
}
/**********************************
Graphical user inferface callbacks:
**********************************/

void updateVisualization (bool atomPositionsChanged = true)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
    state->proteinRenderer->updateProtein();
    state->interactor->resetDragBox();
    ramachandranPlot->updateProtein();
    if ( state->energyCalculator && atomPositionsChanged )
        state->energyCalculator->updateProtein();
    renderWindow->redraw();
    ramachandranPlot->redraw();
}

void updateDihedralAngles(void)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->selectedResidue!=0&&residueDialog->visible_r())
    {
        double phis[1],psis[1];
        state->protein->getDihedralAngles (
            state->protein->getResidueIndex(state->selectedResidue),1,phis,psis
        );
        residueDialogPhiOutput->value(Math::deg(phis[0]));
        residueDialogPsiOutput->value(Math::deg(psis[0]));
    }
}

void updateRamachandranWindow( )
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->selectedResidue!=0)
    {
    	MD::Scalar* phis=new MD::Scalar[1];
    	MD::Scalar* psis=new MD::Scalar[1];
    	state->protein->getDihedralAngles (
	    	state->protein->getResidueIndex(state->selectedResidue),1,phis,psis );

		/* ramachandranPlot redarw after the angles update */
		ramachandranPlot->updateSelectedResidueAngles(phis[0], psis[0]);	
	
		delete[] phis;
		delete[] psis;
		updateDihedralAngles();
	}
}

void redrawRenderWindows(void)
{
    renderWindow->redraw();
    ramachandranPlot->redraw();
}

void updateProteinNow(void)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->updateProtein();
    ramachandranPlot->updateProtein();
    if(state->energyCalculator!=0)
        state->energyCalculator->updateProtein();
    // interactor->resetDragBox();
    if(renderWindow->visible_r())
        renderWindow->redrawNow();
    if(ramachandranPlotWindow->visible_r())
        ramachandranPlot->redrawNow();
}


void updateVolumeDialogGui()
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
    {
        uint classifierNumber = state->energyRenderer->classifier();
        volumeDialogClassifierChoice->value (classifierNumber);
        AtomClassifier *classifier = AtomClassifier::get (classifierNumber);
        if ( classifier )
        {
            snprintf (
                volumeDialogClassifierGroupLabel, buflen,
                "%s Classifier (%u of %u class%s)",
                classifier->name(),
                (volumeDialogClassNumber + 1),
                classifier->numClasses(),
                ( (classifier->numClasses() == 1) ? "" : "es" )
            );
            snprintf (
                volumeDialogColorChoiceLabel, buflen,
                "Color function for the %s class",
                classifier->className(volumeDialogClassNumber)
            );
            volumeDialogClassifierGroup->label (
                volumeDialogClassifierGroupLabel
            );
            volumeDialogColorChoice->label (volumeDialogColorChoiceLabel);
            volumeDialogColorChoice->value (
                state->energyRenderer->colorFunction(volumeDialogClassNumber)
            );
        }
        volumeDialogRadiusSlider->value (
            state->energyRenderer->radiusMultiplier()
        );
        volumeDialogResolutionSlider->value (
            state->energyRenderer->texelsPerAngstrom()
        );
        if ( state->energyRenderer->isGradient() )
        {
            volumeDialogSubsetSumRadioButton->value (0);
            volumeDialogGradientRadioButton->value (1);
        }
        else
        {
            volumeDialogSubsetSumRadioButton->value (1);
            volumeDialogGradientRadioButton->value (0);
        }
        volumeDialogUniformRadioButton->value (0);
        volumeDialogRadiusRadioButton->value (0);
        volumeDialogVDWRadioButton->value (0);
        switch ( state->energyRenderer->radiusType() )
        {
            case EnergyRenderer::UNIFORM_RADIUS:
                volumeDialogUniformRadioButton->value (1);
                break;
            case EnergyRenderer::ATOM_RADIUS:
                volumeDialogRadiusRadioButton->value (1);
                break;
            case EnergyRenderer::VAN_DER_WAALS_RADIUS:
                volumeDialogVDWRadioButton->value (1);
                break;
            default:  break;
        }
        volumeDialogRangeMinInput->value (state->energyRenderer->minNormal());
        volumeDialogRangeMaxInput->value (state->energyRenderer->maxNormal());
        if ( state->energyRenderer->isAutoNormalizing() )
        {
            volumeDialogAutoNormalizeToggle->value (1);
            volumeDialogRangeMinInput->deactivate();
            volumeDialogRangeMaxInput->deactivate();
        }
        else
        {
            volumeDialogAutoNormalizeToggle->value (0);
            volumeDialogRangeMinInput->activate();
            volumeDialogRangeMaxInput->activate();
        }
    }
}


void updateGui(void)
{
    ProteinState *state = curProtein();
    if ( !state )
    {
        // no protein --> disable some items in the menu
        loadAngleFileButton->deactivate();
        saveAngleFileButton->deactivate();
        savePdbFileButton->deactivate();
        saveScreenshotButton->deactivate();
        serverMenu->deactivate();
        editMenu->deactivate();
        viewMenu->deactivate();
        windowsMenu->deactivate();
    }
    else
    {
        // enable those items instead
        loadAngleFileButton->activate();
        saveAngleFileButton->activate();
        savePdbFileButton->activate();
        saveScreenshotButton->activate();
        serverMenu->activate();
        editMenu->activate();
        viewMenu->activate();
        windowsMenu->activate();
    }
    // activate and deactivate don't trigger a redraw in this case
    mainMenuBar->redraw();
    if ( !state ) return;

    /* Update main menu toggles: */
    if(state->client!=0)
    {
        disconnectButton->activate();
        getOptimizationTreeButton->activate();
        getConfigButton->activate();
        getConfigByIdButton->activate();
        getBestConfigButton->activate();
        addConfigButton->activate();
    }
    else
    {
        disconnectButton->deactivate();
        getOptimizationTreeButton->deactivate();
        getConfigButton->deactivate();
        getConfigByIdButton->deactivate();
        getBestConfigButton->deactivate();
        addConfigButton->deactivate();
    }
    if(undoBuffer.canUndo())
        undoButton->activate();
    else
        undoButton->deactivate();
    if(undoBuffer.canRedo())
        redoButton->activate();
    else
        redoButton->deactivate();

    if(structureDialog->visible_r())
        showStructureDialogToggle->set();
    else
        showStructureDialogToggle->clear();
    if(residueDialog->visible_r())
        showResidueDialogToggle->set();
    else
        showResidueDialogToggle->clear();
    if(sequenceDialog->visible_r())
        showSequenceDialogToggle->set();
    else
        showSequenceDialogToggle->clear();
    if(ramachandranPlotWindow->visible_r())
        showRamachandranPlotWindowToggle->set();
    else
        showRamachandranPlotWindowToggle->clear();
    if(renderingDialog->visible_r())
        showRenderingDialogToggle->set();
    else
        showRenderingDialogToggle->clear();
    if ( selectionDialog->visible_r() )
        showSelectionDialogToggle->set();
    else
        showSelectionDialogToggle->clear();
    if ( state->energyCalculator )
    {
        // allow the energy dialog to be shown
        showEnergyDialogToggle->activate();
        if ( energyDialog->visible_r() )
            showEnergyDialogToggle->set();
        else
            showEnergyDialogToggle->clear();
    }
    else
    {
        // energy dialog will not show if protein has no calculator
        showEnergyDialogToggle->clear();
        showEnergyDialogToggle->deactivate();
        if ( energyDialog->visible_r() )
            hideEnergyDialog();
    }

    MD::Protein::StructureSelector selectedStructure=state->interactor->getStructure();

    /* Update residue dialog: */
    if(residueDialog->visible_r())
    {
        std::pair<int,int> residueIndexRange=state->protein->getResidueIndexRange();
        residueDialogResidueIndexSlider->range(residueIndexRange.first-1,residueIndexRange.second);
        residueDialogResidueTypeChoice->set_output();
        if(state->selectedResidue!=0)
        {
            residueDialogResidueIndexSlider->value(state->selectedResidue->getPdbResidueIndex());
            residueDialogResidueTypeChoice->activate();
            residueDialogResidueTypeChoice->value(state->selectedResidue->getType());
            residueDialogPdbResidueNameOutput->activate();
            residueDialogPdbResidueNameOutput->value(state->selectedResidue->getPdbResidueName());
            residueDialogStructureTypeChoice->activate();
            switch(state->selectedResidue->getSecondaryStructure()->getStructureType())
            {
                case MD::Protein::SecondaryStructure::COIL:
                    residueDialogStructureTypeChoice->value(0);
                    break;

                case MD::Protein::SecondaryStructure::ALPHA_HELIX:
                    residueDialogStructureTypeChoice->value(1);
                    break;

                case MD::Protein::SecondaryStructure::BETA_STRAND:
                    residueDialogStructureTypeChoice->value(2);
                    break;

				default:
					residueDialogStructureTypeChoice->value(3);
		    		break;
            }
            residueDialogPhiOutput->activate();
            residueDialogPsiOutput->activate();
            if ( state->selectedResidue->getType() != MD::Protein::Residue::UNK &&
                 state->selectedResidue->getSecondaryStructure()->getStructureType() ==
                    MD::Protein::SecondaryStructure::COIL )
                residueDialogRandomizeAnglesButton->activate();
            else
                residueDialogRandomizeAnglesButton->deactivate();
            updateDihedralAngles();
			labelResidueNameToggle->value(
				state->proteinRenderer->getDrawResidueName());
			showResidueAnglesToggle->value(
				state->proteinRenderer->getDrawResidueAngles());
        	labelAtomNamesToggle->value(
				state->proteinRenderer->getDrawAtomNames());
			showPhobicButton->value(state->proteinRenderer->getShowHydrophobic());
			showPhilicButton->value(state->proteinRenderer->getShowHydrophilic());	 
			showDisulfideButton->value(state->proteinRenderer->getShowDisulfide());	 
		
		}
        else
        {
            residueDialogResidueIndexSlider->value(residueIndexRange.first-1);
			residueDialogResidueTypeChoice->value(23);
			residueDialogResidueTypeChoice->deactivate();
            residueDialogPdbResidueNameOutput->value(0);
			residueDialogPdbResidueNameOutput->deactivate();
            residueDialogStructureTypeChoice->value(3);
            residueDialogStructureTypeChoice->deactivate();
            residueDialogPhiOutput->deactivate();
            residueDialogPsiOutput->deactivate();
            residueDialogRandomizeAnglesButton->deactivate();
			labelResidueNameToggle->value(0);
			showResidueAnglesToggle->value(0);
        	labelAtomNamesToggle->value(0);
			showPhobicButton->value(0);	 
			showPhilicButton->value(0);	 
			showDisulfideButton->value(0);	 
			
        }
    }

    /* Update structure dialog: */
    if(structureDialog->visible_r())
    {
        structureDialogStructureIndexSlider->range(-1,state->protein->getNumStructures()-1);
        structureDialogStructureTypeChoice->set_output();
        structureDialogFirstResidueIndexOutput->range(0,state->protein->getNumResidues()-1);
        structureDialogLastResidueIndexOutput->range(0,state->protein->getNumResidues()-1);
        if(selectedStructure.isValid())
        {
            /* Set widgets in dialog header: */
            structureDialogStructureIndexSlider->value(selectedStructure.getStructureIndex());
            structureDialogStructureTypeChoice->activate();
            switch(selectedStructure.getStructureType())
            {
                case MD::Protein::SecondaryStructure::COIL:
                    structureDialogStructureTypeChoice->value(0);
                    break;

                case MD::Protein::SecondaryStructure::ALPHA_HELIX:
                    structureDialogStructureTypeChoice->value(1);
                    break;

                case MD::Protein::SecondaryStructure::BETA_STRAND:
                    structureDialogStructureTypeChoice->value(2);
                    break;

                default:
                    break;
            }

            /* Set widgets in residue range group: */
            structureDialogFirstResidueIndexOutput->activate();
            structureDialogFirstResidueIndexOutput->value(selectedStructure.getFirstResidueIndex());
            structureDialogLastResidueIndexOutput->activate();
            structureDialogLastResidueIndexOutput->value(selectedStructure.getFirstResidueIndex()+selectedStructure.getNumResidues()-1);

            /* Set widgets in beta strand shape adjustment group: */
            if(selectedStructure.getStructureType()==MD::Protein::SecondaryStructure::BETA_STRAND)
            {
                structureDialogCurlRoller->activate();
                structureDialogTwistRoller->activate();
                structureDialogPleatRoller->activate();
                structureDialogBraidRoller->activate();
                structureDialogFlattenBetaStrandButton->activate();
            }
            else
            {
                structureDialogCurlRoller->deactivate();
                structureDialogTwistRoller->deactivate();
                structureDialogPleatRoller->deactivate();
                structureDialogBraidRoller->deactivate();
                structureDialogFlattenBetaStrandButton->deactivate();
            }
            if(selectedStructure.getStructureType()!=MD::Protein::SecondaryStructure::NONE)
                structureDialogResetBetaStrandButton->activate();
            else
                structureDialogResetBetaStrandButton->deactivate();
        }
        else
        {
            structureDialogStructureIndexSlider->value(-1);
            structureDialogStructureTypeChoice->deactivate();
            structureDialogFirstResidueIndexOutput->deactivate();
            structureDialogLastResidueIndexOutput->deactivate();
            structureDialogCurlRoller->deactivate();
            structureDialogTwistRoller->deactivate();
            structureDialogPleatRoller->deactivate();
            structureDialogBraidRoller->deactivate();
            structureDialogFlattenBetaStrandButton->deactivate();
            structureDialogResetBetaStrandButton->deactivate();
        }
    }

    /* Update rendering dialog: */
    if(renderingDialog->visible_r())
    {
        renderingDialogDrawAtomsToggle->value(state->proteinRenderer->getDrawAtoms());
        renderingDialogDrawBondsToggle->value(state->proteinRenderer->getDrawBonds());
        renderingDialogDrawBackboneToggle->value(state->proteinRenderer->getDrawBackboneRibbon());
        renderingDialogDrawCartoonToggle->value(state->proteinRenderer->getDrawCartoon());
        renderingDialogDrawCPKToggle->value(state->proteinRenderer->getDrawCPK());
        renderingDialogDrawTubeToggle->value(state->proteinRenderer->getDrawTube());
        renderingDialogDrawLineToggle->value(state->proteinRenderer->getDrawLine());
        renderingDialogDrawHydrogenBondsToggle->value (
            state->proteinRenderer->getDrawHydrogenBonds()
        );
        renderingDialogDrawHydrogenBondSitesToggle->value (
            state->proteinRenderer->getDrawHydrogenBondSites()
        );
        renderingDialogDrawHydrogenCagesToggle->value (
            state->proteinRenderer->getDrawHydrogenCages()
        );
        renderingDialogDrawCollisionsToggle->value (
            state->proteinRenderer->getDrawCollisions()
        );

        if(selectedStructure.isValid())
        {
            renderingDialogStructureDrawAtomsToggle->activate();
            renderingDialogStructureDrawAtomsToggle->value (
                state->proteinRenderer->getDrawAtoms(selectedStructure)
            );
            renderingDialogStructureDrawBondsToggle->activate();
            renderingDialogStructureDrawBondsToggle->value (
                state->proteinRenderer->getDrawBonds(selectedStructure)
            );
            renderingDialogStructureDrawBackboneToggle->activate();
            renderingDialogStructureDrawBackboneToggle->value (
                state->proteinRenderer->getDrawBackboneRibbon(selectedStructure)
            );
            renderingDialogStructureDrawCartoonToggle->activate();
            renderingDialogStructureDrawCartoonToggle->value (
                state->proteinRenderer->getDrawCartoon(selectedStructure)
            );
            renderingDialogStructureDrawCPKToggle->activate();
            renderingDialogStructureDrawCPKToggle->value (
                state->proteinRenderer->getDrawCPK(selectedStructure)
            );
            renderingDialogStructureDrawTubeToggle->activate();
            renderingDialogStructureDrawTubeToggle->value (
                state->proteinRenderer->getDrawTube(selectedStructure)
            );
            renderingDialogStructureDrawLineToggle->activate();
            renderingDialogStructureDrawLineToggle->value (
                state->proteinRenderer->getDrawLine(selectedStructure)
            );
            renderingDialogStructureDrawHydrogenBondSitesToggle->activate();
            renderingDialogStructureDrawHydrogenBondSitesToggle->value (
                state->proteinRenderer->getDrawHydrogenBondSites(selectedStructure)
            );
            renderingDialogStructureDrawHydrogenCagesToggle->activate();
            renderingDialogStructureDrawHydrogenCagesToggle->value (
                state->proteinRenderer->getDrawHydrogenCages(selectedStructure)
            );
            if(state->proteinRenderer->getDrawHydrogenCages(selectedStructure))
            {
                renderingDialogStructureDrawLargeHydrogenCagesToggle->activate();
                renderingDialogStructureDrawLargeHydrogenCagesToggle->value (
                    state->proteinRenderer->getDrawLargeHydrogenCages(selectedStructure)
                );
            }
            else
                renderingDialogStructureDrawLargeHydrogenCagesToggle->deactivate();
        }
        else
        {
            renderingDialogStructureDrawAtomsToggle->deactivate();
            renderingDialogStructureDrawBondsToggle->deactivate();
            renderingDialogStructureDrawBackboneToggle->deactivate();
            renderingDialogStructureDrawCartoonToggle->deactivate();
            renderingDialogStructureDrawCPKToggle->deactivate();
            renderingDialogStructureDrawTubeToggle->deactivate();
            renderingDialogStructureDrawLineToggle->deactivate();
            renderingDialogStructureDrawHydrogenBondSitesToggle->deactivate();
            renderingDialogStructureDrawHydrogenCagesToggle->deactivate();
            renderingDialogStructureDrawLargeHydrogenCagesToggle->deactivate();
        }
    }

    /* Update sequence dialog: */
    if(sequenceDialog->visible_r())
    {
        sequenceDialog->redraw();
    }

    /* Update Ramachandran plot window: */
    if(ramachandranPlotWindow->visible_r())
    {
        ramachandranPlot->updateProtein();
        ramachandranPlot->redraw();
    }

    /* Update energy visualization dialog: */
    if ( state->energyCalculator && energyDialog->visible_r() )
    {
        energyDialogEnergyValueOutput->value (
            state->energyCalculator->calcEnergy()
        );
        energyDialogDisplayEnergyToggle->value (
            state->proteinRenderer->getEnergyUpdate()
        );
        energyUpdateRatInput->value (state->engUpdateRate);
        energyDialogVisualizeEnergyToggle->value (state->visualizeEnergy);
        energyDialogRangeMinInput->value (state->visualizeEnergyMinRange);
        energyDialogRangeMaxInput->value (state->visualizeEnergyMaxRange);
        for ( int i = 0; i < energyLibrary->getNumEnergyComponents(); ++i )
            static_cast<Fl_Check_Button*>(
                energyDialogComponentGroup->child(i)
            )->value (
                state->energyCalculator->getEnergyComponentState(i)
            );
        if ( volumeDialog->visible_r() && state->energyRenderer )
            updateVolumeDialogGui();
    }
}


void hideAllDialogs()
{
    residueDialog->hide();
    structureDialog->hide();
    sequenceDialog->hide();
    ramachandranPlotWindow->hide();
    renderingDialog->hide();
    selectionDialog->hide();
    if ( energyLibrary )
        hideEnergyDialog();
}


#ifdef FLTK_DOUBLE
void windowCB(Fl_Double_Window* window,void* cbData)
#else
void windowCB(Fl_Window* window,void* cbData)
#endif
{
    if(window==mainWindow)
    {
        if(Fl::event_key()!=FL_Escape) // Don't close main window on escape
        {
            if(undoBuffer.isSaved()||fl_ask("Do you really want to exit the program and discard all unsaved changes?"))
            {
                /* Hide all windows (thereby exiting the program quite cleanly): */
                hideAllDialogs();
                mainWindow->hide();
            }
        }
    }
    else if ( window == recordDialog )
        hideRecordDialog();
    else if ( window == energyDialog )
        hideEnergyDialog();
    else
        window->hide();

    updateGui();
}

static MD::Protein *loadPrediction (const char *filename, int setprotein)
{
    using namespace MD;
    Protein *protein = 0;
    if ( filename )
    {
        strncpy (inputfilename, filename, inputfilename_size);
        inputfilename[inputfilename_size-1] = 0;
        try
        {
            /* Build a protein model: */
            ReadStandards ("Standards/");
            protein = ReadPredictionFile (filename, setprotein);
			// Silently assumed all pred file are ready for energy calculation
			validPdbForEnergy = true;
        }
        catch ( std::runtime_error error )
        {
            fl_alert ("Error loading prediction file:\n%s", error.what());
            msg.Error (PREDID, error.what());
        }
    }
    return protein;
}


void loadPredictionFileCB(Fl_Menu_* menuBar,void* cbData)
{
    using namespace MD;

    /* Show a file select dialog to choose a prediction file: */
    char* predictionFilename=fl_file_chooser("Load Prediction File...","*.pred",0);
	if(predictionFilename!=0)
    {
		validPdbForEnergy = true;
        Protein *protein = loadPrediction (predictionFilename, 1);
        if ( protein ) {
	   addProtein (protein,selectionDialogBrowser->size());
           initBuild (protein);
	   printf("returned from initBuild\n");
           if ( offlineBuildBeta == 1 )
           {
              BuildBeta();
              exit(-1);
           }
        }
    	else inputfilename[0] = 0;
		
    }
    updateGui();
}

void loadAngleFileCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    /* Show a file select dialog to choose an angle file: */
    char* angleFilename=fl_file_chooser("Load Angle File...","*.dih",0);
    if(angleFilename!=0)
        state->interactor->readCoilRegionAngles(angleFilename);

    updateGui();
    renderWindow->redraw();
}

void saveAngleFileCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    /* Show a file select dialog to choose an angle file: */
    char* angleFilename=fl_file_chooser("Save Angle File...","*.dih",0);
    if(angleFilename!=0)
        state->interactor->writeCoilRegionAngles(angleFilename);
}

void loadProteins (const char *filename)
{
    using namespace MD;
    Protein *protein = 0;
    Stride* pred = new Stride;
	char* chain = new char[1];
	int id = selectionDialogBrowser->size();
	int models = parseModels(filename);
	const char* chainList = parseChains(filename);
	char* selectedModel = new char[models];
	char* selectedChain = new char[strlen(chainList)];

    if ( filename )
    {
		if(!pred->stridePrediction(filename))
		{
			fl_alert ("Not able to predict the second structure!");
			return;
		}
		delete pred;
		
		strncpy (inputfilename, filename, inputfilename_size);
        inputfilename[inputfilename_size-1] = 0;
		try
        {
            /* Check the requirements of energy comuutation for every PDB file: */
            validPdbForEnergy = checkPdbFile (filename);
        }
        catch ( std::runtime_error error )
        {
            /* Silently disable the EnergyDialog */
            validPdbForEnergy = false;
            fl_message ("PDB file not valid for energy computation.");
        }
        
		if ( energyLibrary && !validPdbForEnergy)
            msg.Warn (PDBID, "This PDB file is not ready for Energy Calculation.");

		try
        {
   			// chain selection
			if(strlen(chainList)>1)
			{
			msg.Info(PDBID, "This is a X-ray diffraction solved structure. ", filename);
			strcpy(selectedChain,configFile->retrieveString("/ProteinCreator/chainSelection","").c_str());	
			if (!strlen(selectedChain)) 
			{
				// all the chains will be loaded
				for(int j =0; j< strlen(chainList); j++)
				{
					strncpy(chain, &chainList[j], 1);
					chain[1]='\0';
					// Load a protein model: Disable ACE/NME cap if this pdb file is not valid
           			if(validPdbForEnergy)
						protein = parsePdbFile (filename, chain);
					else 
						protein = parsePdbFile ( filename, false, chain);
					if ( protein ) addProtein (protein, id);
				}
			}
			else
			{
				selectedChain = strtok( selectedChain, " " );
				while( selectedChain != NULL )
				{
					if (sscanf( selectedChain, "%c", chain ))
						chain[1]='\0';

					selectedChain = strtok( NULL, " " );
			    	bool found = false;
					// load selected chain
					for(int j =0; j< strlen(chainList); j++)
        			{
						// make sure the existness of selected chain
						if(strncmp(chain, &chainList[j], 1)==0)
						{
            				if(validPdbForEnergy)
								protein = parsePdbFile (filename, chain);
							else 
								protein = parsePdbFile ( filename, false, chain);
							if ( protein ) addProtein (protein, id);
							found = true;
						}
					}
					if(!found)
						msg.Warn(PDBID, "This chain (", chain, ") was not found! Check configuration file.");
				}
			}
			}
			else // model selection
			{
				msg.Info(PDBID, "The structure was determined using NMR spectroscopy. ", filename);
				strcpy(selectedModel,configFile->retrieveString("/ProteinCreator/multipleModel","").c_str());	
				if ( !strlen(selectedModel) )
           		{
					// all the models will be loaded
					for (int i = 1; i< models+1; i++)
					{
						if(validPdbForEnergy)
							protein = parsePdbFile (filename, i);
						else 
							protein = parsePdbFile ( filename, false, 0, i);
						if ( protein ) addProtein (protein, id);
 					}
				}
				else
				{
 					int number;
					selectedModel = strtok( selectedModel, " " );
					while( selectedModel != NULL )
					{
						if (!sscanf(selectedModel , "%d", &number ))
							msg.Warn(PDBID, "No number found!");

            	  		selectedModel = strtok( NULL, " " );
						
						if(number <0 || number > models)
						{
							sprintf(msg.buf, "This model (%d) was not found! Check configuration file.", number);
							msg.Warn(PDBID, msg.buf);
							continue;
						}
						if(validPdbForEnergy)
							protein = parsePdbFile (filename, number);
						else 
							protein = parsePdbFile ( filename, false, 0, number);
						if ( protein ) addProtein (protein, id);
					}
				}
			}
		}
        catch ( std::runtime_error error )
        {
            fl_alert ("Error loading PDB file:\n%s", error.what());
            msg.Error (PDBID, error.what());
       }
    }
//	pred->clean();
//	delete pred;
	delete [] chain;
	delete [] selectedModel;
	delete [] selectedChain;
}

static MD::Protein *loadPdb (const char *filename)
{
    using namespace MD;
    Protein *protein = 0;
    Stride* pred = new Stride;
    if ( filename )
    {
		if(!pred->stridePrediction(filename))
		{
			fl_alert ("Not able to predict the second structure!");
			return 0;
		}
        strncpy (inputfilename, filename, inputfilename_size);
        inputfilename[inputfilename_size-1] = 0;
        try
        {
            /* Check the requirements of energy comuutation for every PDB file: */
            validPdbForEnergy = checkPdbFile (filename);
        }
        catch ( std::runtime_error error )
        {
            /* Silently disable the EnergyDialog */
            validPdbForEnergy = false;
            fl_message ("PDB file not valid for energy computation.");
        }
		if ( energyLibrary && !validPdbForEnergy)
            msg.Warn (PDBID, "This PDB file is not ready for Energy Calculation.");
        try
        {
            /* Load a protein model: */
			if(validPdbForEnergy)
				protein = parsePdbFile (filename);
			else 
				protein = parsePdbFile (filename, false);
        }
        catch ( std::runtime_error error )
        {
            fl_alert ("Error loading PDB file:\n%s", error.what());
            msg.Error (PDBID, error.what());
        }
    }
	//pred->clean();
	delete pred;
    return protein;
}

void loadPdbFileCB(Fl_Menu_* menuBar,void* cbData)
{
    using namespace MD;

    /* Show a file select dialog to choose a PDB file: */
    char* pdbFilename=fl_file_chooser("Load PDB File...","*.pdb",0);
    if(pdbFilename!=0)
    {
		const char* chainList = parseChains(pdbFilename);
		msg.Info(PDBID, "Loading ", pdbFilename);
        // loadPdb for single protein and loadProteins for multiple protein
		if (strlen(chainList) > 1 || parseModels( pdbFilename ) > 0)
		{
			loadProteins (pdbFilename);
		}
		else
		{
			Protein *protein = loadPdb (pdbFilename);
         	if ( protein ) addProtein (protein, selectionDialogBrowser->size());
		}
    }
    updateGui();
}


void savePdbFileCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    /* Show a file select dialog to choose a PDB file: */
    char* pdbFilename=fl_file_chooser("Save PDB File...","*.pdb",0);
    if(pdbFilename!=0)
    {
        MD::writePdbFile(*state->protein,pdbFilename);
        undoBuffer.saveCurrentState();
    }

    updateGui();
}


void d_show_callback(void)
{
  	predfc->show();
  	while (predfc->visible())
    	Fl::wait();
	predFileInput->value(predfc->value());
}

void b_show_callback(void)
{
  	pdbfc->show();
  	while (pdbfc->visible())
    	Fl::wait();
	pdbFileInput->value(pdbfc->value());
}

static void closeWindowCB( Fl_Widget* , void* cbData)
{
    using namespace MD;
	Fl_Window *cb = (Fl_Window *)cbData;

	Protein *protein =0;;
	if( predfc->value() && pdbfc->value() ) // load both pdb and pred files
	{
//	 protein = loadPdb (pdbfc->value());
     loadPrediction (predfc->value(), 0);
	 protein = loadPdb (pdbfc->value());
	
	}
	else if (predfc->value()) // pred file only
	{
		validPdbForEnergy = true;
        protein = loadPrediction (predfc->value(), 1);
	}
	else if (pdbfc->value())
         fl_message ("Load a Pred file first!");
	
	if ( protein )
	{
	   addProtein (protein,selectionDialogBrowser->size());
       initBuild (protein);
       if ( offlineBuildBeta == 1 )
       {
              BuildBeta();
              exit(-1);
       }
    }
    else inputfilename[0] = 0;
    updateGui();

  	delete cb;
}

void copyDialogCB(Fl_Button* button,void* cbData)
{
    int begin, end;
	Fl_Window *cb = (Fl_Window *)cbData;
    ProteinState *state = curProtein();
    if ( !state ) return;
		MD::ProteinFactory *factory = new MD::ProteinFactory;

	if(state->interactor->isValid()) {
		//there is a bug of residue index about 0 base, the fix breaks IK so.... 
		begin = state->interactor->getStructure().getFirstResidueIndex();
		end =	state->interactor->getStructure().getLastResidueIndex();
		MD::Protein *newprotein = factory->copy(begin, end);
		if(newprotein)
			addProtein (newprotein, selectionDialogBrowser->size());
		delete factory;
	}
	else {
		if (copyBeginIndexInput->value()!=copyEndIndexInput->value()) {
			if (copyBeginIndexInput->value() >= copyEndIndexInput->value()) {
				fl_message ("The end index should greater than the begin index.");
				return;
				}
			begin = (int)copyBeginIndexInput->value();
			end =	(int)copyEndIndexInput->value();
			MD::Protein *newprotein = factory->copy(begin, end);
			if(newprotein)
				addProtein (newprotein,selectionDialogBrowser->size());
			delete factory;
		}
	}
	selectProtein ( curProtein ());
    selectionDialogBrowser->deselect();
	selectionDialogBrowser->select(state->proteinId+1);
    renderWindow->centerView();
	renderWindow->redraw();
    updateGui();
	cb->hide();
  	//delete cb;
}

void copyGuiCB(Fl_Menu_*, void*)
{
	Fl_Window* copyWin = new Fl_Window(250, 50, "Copy/Extract ");
    copyBeginIndexInput = new Fl_Value_Input(50, 10, 40, 25,"From: ");
    copyEndIndexInput = new Fl_Value_Input(120, 10, 40, 25, "To: ");
 	
	Fl_Button* copy = new Fl_Button(180, 10, 60, 25, "Copy");
  		copy->callback((Fl_Callback*)copyDialogCB, (void*) copyWin);
 	
	copyWin->end();
 	copyWin->show();
 	Fl::run();
}

void BuildBetaGuiCB(Fl_Menu_*, void*)
{
    //using namespace MD;
	pdbfc= new Fl_File_Chooser(".", "*.pdb", Fl_File_Chooser::SINGLE,	"PDB file");;
	predfc= new Fl_File_Chooser(".", "*.pred", Fl_File_Chooser::SINGLE,	"Pred file");;
	
	Fl_Window* loadWin = new Fl_Window(360, 100, "BuildBeta ");
  	predFileInput = new Fl_Input(50, 10, 270, 25, "Pred: ");
  	pdbFileInput = new Fl_Input(50, 40, 270, 25, "PDB: ");
	//pdbFile->value( fc-> value()); 
  	
	
	Fl_Button *button1 = new Fl_Button(325, 10, 25, 25);
  		button1->labelcolor(FL_BLUE);
  		button1->label(".pred");
  		button1->labelsize(10);
  		button1->callback((Fl_Callback *)d_show_callback);
    
	Fl_Button *button2 = new Fl_Button(325, 40, 25, 25);
  		button2->labelcolor(FL_BLUE);
  		button2->label(".pdb");
  		button2->labelsize(10);
  		button2->callback((Fl_Callback *)b_show_callback);
 	
	Fl_Button* done = new Fl_Button(250, 70, 70, 25, "Done");
  		done->callback((Fl_Callback*)closeWindowCB, (void*) loadWin);
 	
	loadWin->resizable(pdbFileInput);
	loadWin->end();
 	loadWin->show();
 	Fl::run();
}


void saveScreenshotCB(Fl_Menu_* menuBar,void* cbData)
{
    /* Show a file select dialog to choose a PPM file: */
    char* imageFilename=fl_file_chooser("Save Screenshot...","*.ppm",0);
    if(imageFilename!=0)
        renderWindow->saveScreenshot(imageFilename);

    updateGui();
}

void quitCB(Fl_Menu_* menuBar,void* cbData)
{
    if(undoBuffer.isSaved()||fl_ask("Do you really want to exit the program and discard all unsaved changes?"))
    {
        /* Hide all windows (thereby exiting the program quite cleanly): */
        hideAllDialogs();
        mainWindow->hide();
    }
}

void connectCB(Fl_Menu_* menuBar,void* cbData)
{
    /* Show connect dialog: */
    connectDialog->show();
}

void disconnectCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->client!=0)
    {
        /* Disconnect from server: */
        delete state->client;
        state->client=0;
    }

    updateGui();
}

void getOptimizationTreeCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->client!=0)
    {
        try
        {
            state->client->getOptimizationTree("OptimizationTree.txt");
        }
        catch(std::runtime_error error)
        {
            fl_alert(error.what());
        }
    }
}

void getConfigCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->client!=0)
    {
        try
        {
            /* Lock server and get number of configurations in queue: */
            state->client->lockServerQueue();
            int numConfigs=state->client->getNumConfigurations();

            /* Set configuration dialog options: */
            getConfigDialogConfigIndexSlider->range(0,numConfigs-1);
            getConfigDialogConfigIndexSlider->value(0);

            /* Show get configuration dialog: */
            getConfigDialog->show();
        }
        catch(std::runtime_error error)
        {
            fl_alert(error.what());
            try
            {
                /* Unlock server just to be sure: */
                state->client->unlockServerQueue();
            }
            catch(...)
            {
                /* Oops... */
                fl_alert(error.what());
            }
        }
    }
}

void getConfigByIdCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->client!=0)
    {
        /* Show get configuration by ID dialog: */
        getConfigByIdDialog->show();
    }
}

void getBestConfigCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->client!=0)
    {
        try
        {
            state->proteinId=state->client->getBestConfiguration(state->protein);
        }
        catch(std::runtime_error error)
        {
            fl_alert(error.what());
        }
    }

    state->proteinRenderer->updateProtein();
    state->interactor->resetDragBox();
    renderWindow->redraw();

    updateGui();
}

void addConfigCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->client!=0)
    {
        try
        {
            state->proteinId=state->client->addConfiguration(state->proteinId,state->protein);
        }
        catch(std::runtime_error error)
        {
            fl_alert(error.what());
        }
    }
}

void undoCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    undoBuffer.undo();
    state->proteinRenderer->updateProtein();
    state->interactor->resetDragBox();
    renderWindow->redraw();

    updateGui();
}

void redoCB(Fl_Menu_* menuBar,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    undoBuffer.redo();
    state->proteinRenderer->updateProtein();
    state->interactor->resetDragBox();
    renderWindow->redraw();

    updateGui();
}


static void centerView()
{
    renderWindow->centerView();
    renderWindow->redraw();

    updateGui();
}


void centerViewCB(Fl_Menu_* menuBar,void* cbData)
{
    centerView();
}


void pushViewCB(Fl_Menu_* menuBar,void* cbData)
{
    renderWindow->pushViewSpec();
    renderWindow->redraw();

    updateGui();
}

void popViewCB(Fl_Menu_* menuBar,void* cbData)
{
    renderWindow->popViewSpec();
    renderWindow->redraw();

    updateGui();
}

void showResidueDialogCB(Fl_Menu_* menuBar,void* cbData)
{
    if(showResidueDialogToggle->value())
        residueDialog->show();
    else
        residueDialog->hide();

    updateGui();
}

void showStructureDialogCB(Fl_Menu_* menuBar,void* cbData)
{
    if(showStructureDialogToggle->value())
        structureDialog->show();
    else
        structureDialog->hide();

    updateGui();
}

void showSequenceDialogCB(Fl_Menu_* menuBar,void* cbData)
{
    if(showSequenceDialogToggle->value())
        sequenceDialog->show();
    else
        sequenceDialog->hide();

    updateGui();
}

void showRamachandranPlotWindowCB(Fl_Menu_* menuBar,void* cbData)
{
    if(showRamachandranPlotWindowToggle->value())
        ramachandranPlotWindow->show();
    else
        ramachandranPlotWindow->hide();

    updateGui();
}

void showRenderingDialogCB(Fl_Menu_* menuBar,void* cbData)
{
    if(showRenderingDialogToggle->value())
        renderingDialog->show();
    else
        renderingDialog->hide();

    updateGui();
}


void showSelectionDialogCB (Fl_Menu_* menuBar, void *cbData)
{
    if ( showSelectionDialogToggle->value() )
        selectionDialog->show();
    else
        selectionDialog->hide();
    updateGui();
}


void showEnergyDialogCB(Fl_Menu_* menuBar,void* cbData)
{
    if(energyLibrary!=0)
    {
        if(showEnergyDialogToggle->value())
            energyDialog->show();
        else
            hideEnergyDialog();
    }

    updateGui();
}

void aboutCB(Fl_Menu_* menuBar,void* cbData)
{
    aboutDialog->show();
}


#ifdef FLTK_DOUBLE
void connectDialogCB(Fl_Double_Window* window,void* cbData)
#else
void connectDialogCB(Fl_Window* window,void* cbData)
#endif
{
    connectDialog->hide();
    updateGui();
}

void connectDialogOkCB(Fl_Return_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    try
    {
        /* Connect to server: */
        state->client = new ProteinClient (
            connectDialogServerNameText->value(),
            int(connectDialogServerPortInput->value())
        );
    }
    catch(std::runtime_error error)
    {
        fl_alert(error.what());
        delete state->client;
        state->client=0;
    }

    connectDialog->hide();
    updateGui();
}

void connectDialogCancelCB(Fl_Button* button,void* cbData)
{
    connectDialog->hide();
    updateGui();
}

#ifdef FLTK_DOUBLE
void getConfigDialogCB(Fl_Double_Window* window,void* cbData)
#else
void getConfigDialogCB(Fl_Window* window,void* cbData)
#endif
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    try
    {
        /* Unlock server queue: */
        state->client->unlockServerQueue();
    }
    catch(std::runtime_error error)
    {
        fl_alert(error.what());
    }

    getConfigDialog->hide();
    updateGui();
}

void getConfigDialogOkCB(Fl_Return_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    try
    {
        /* Get configuration: */
        state->proteinId = state->client->getConfiguration (
            int(getConfigDialogConfigIndexSlider->value()),
            state->protein
        );
    }
    catch(std::runtime_error error)
    {
        fl_alert(error.what());
    }

    try
    {
        /* Unlock server queue: */
        state->client->unlockServerQueue();
    }
    catch(std::runtime_error error)
    {
        fl_alert(error.what());
    }

    state->proteinRenderer->updateProtein();
    state->interactor->resetDragBox();
    renderWindow->redraw();

    getConfigDialog->hide();
    updateGui();
}

void getConfigDialogCancelCB(Fl_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    try
    {
        /* Unlock server queue: */
        state->client->unlockServerQueue();
    }
    catch(std::runtime_error error)
    {
        fl_alert(error.what());
    }

    getConfigDialog->hide();
    updateGui();
}

#ifdef FLTK_DOUBLE
void getConfigByIdDialogCB(Fl_Double_Window* window,void* cbData)
#else
void getConfigByIdDialogCB(Fl_Window* window,void* cbData)
#endif
{
    getConfigByIdDialog->hide();
    updateGui();
}

void getConfigByIdDialogOkCB(Fl_Return_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    try
    {
        /* Get configuration: */
        unsigned int id=(unsigned int)(getConfigByIdDialogConfigurationIdInput->value());
        state->client->getConfigurationById(id,state->protein);
        state->proteinId=id;
    }
    catch(std::runtime_error error)
    {
        fl_alert(error.what());
    }

    state->proteinRenderer->updateProtein();
    state->interactor->resetDragBox();
    renderWindow->redraw();

    getConfigByIdDialog->hide();
    updateGui();
}

void getConfigByIdDialogCancelCB(Fl_Button* button,void* cbData)
{
    getConfigByIdDialog->hide();
    updateGui();
}

void residueDialogResidueIndexCB(Fl_Value_Slider* slider,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->selectedResidue=state->protein->pickResidue(int(slider->value()));
    renderWindow->redraw();

    updateGui();
}

void residueDialogStructureTypeCB(Fl_Choice* choice,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    MD::Protein::SecondaryStructure::StructureType types[] = {
        MD::Protein::SecondaryStructure::COIL,
        MD::Protein::SecondaryStructure::ALPHA_HELIX,
        MD::Protein::SecondaryStructure::BETA_STRAND
    };

    state->protein->changeResidueStructureType(state->selectedResidue,types[choice->value()]);
    configFile->setCurrentSection("/ProteinRenderer");
    state->proteinRenderer->updateStructureFlags();
    configFile->setCurrentSection("/");
    state->interactor->selectStructure(state->protein->pickStructure(state->selectedResidue));
    renderWindow->redraw();

    updateGui();
}

void residueDialogRandomizeAnglesCB(Fl_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if ( state->selectedResidue != 0 &&
         state->selectedResidue->getType() != MD::Protein::Residue::UNK &&
         state->selectedResidue->getSecondaryStructure()->getStructureType() ==
            MD::Protein::SecondaryStructure::COIL )
    {
        /* Create random displacement values for residue's dihedral angles: */
        MD::Scalar deltaPhi=MD::Scalar(2)*Math::Constants<MD::Scalar>::pi*MD::Scalar(rand())/MD::Scalar(RAND_MAX);
        MD::Scalar deltaPsi=MD::Scalar(2)*Math::Constants<MD::Scalar>::pi*MD::Scalar(rand())/MD::Scalar(RAND_MAX);

        /* Apply the change: */
        undoBuffer.startInteraction (
            state->protein,
            state->selectedResidue,
            state->interactor->getUpdateDirection()
        );
        state->protein->changeDihedralAngles (
            state->protein->getResidueIndex(state->selectedResidue),
            1,
            &deltaPhi,
            &deltaPsi,
            state->interactor->getUpdateDirection()
        );
        undoBuffer.finishInteraction();

        /* Update visualization state: */
        state->interactor->resetDragBox();
        updateVisualization();

   }

    updateGui();
}

void showHydrophobicButtonCB(Fl_Light_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
	state->proteinRenderer->setShowHydrophobic(button->value());
	checkShowToggles(state);
	renderWindow->redraw();
}
void showHydrophilicButtonCB(Fl_Light_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
	state->proteinRenderer->setShowHydrophilic(button->value());
	checkShowToggles(state);
	renderWindow->redraw();
}
void showDisulfideButtonCB(Fl_Light_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
	state->proteinRenderer->setShowDisulfide(button->value());
	checkShowToggles(state);
	renderWindow->redraw();
}
void labelAtomNamesCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
	
	state->proteinRenderer->setDrawAtomNames(labelAtomNamesToggle->value());
    renderWindow->redraw();
    updateGui();
}

void labelResidueNameCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
		
	state->proteinRenderer->setDrawResidueName(labelResidueNameToggle->value());
    renderWindow->redraw();
    updateGui();
}

void showResidueAnglesCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
	
	state->proteinRenderer->setDrawResidueAngles(showResidueAnglesToggle->value());
    renderWindow->redraw();
    updateGui();
}





void structureDialogStructureIndexCB(Fl_Value_Slider* slider,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->interactor->selectStructure(state->protein->pickStructure(int(slider->value())));
    updateVisualization();

    updateGui();
}

void structureDialogCurlCB(MyRoller* roller,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    static double lastCurlValue=0.0;
    MD::Scalar curlDelta=MD::Scalar(Math::rad(roller->value()-lastCurlValue));
    lastCurlValue=roller->value();

    if(state->interactor->isBetaStrand())
    {
        int numResidues=state->interactor->getStructure().getNumResidues();
        MD::Scalar* deltaPhis=new MD::Scalar[numResidues];
        MD::Scalar* deltaPsis=new MD::Scalar[numResidues];
        for(int i=0;i<numResidues;++i)
        {
            deltaPhis[i]=-curlDelta;
            deltaPsis[i]=curlDelta;
            curlDelta=-curlDelta;
        }
        state->interactor->changeDihedralAngles(deltaPhis,deltaPsis);
        delete[] deltaPhis;
        delete[] deltaPsis;

        /* Update visualization state: */
        updateVisualization();

    }

    updateGui();
}

void structureDialogTwistCB(MyRoller* roller,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    static double lastTwistValue=0.0;
    MD::Scalar twistDelta=MD::Scalar(Math::rad(roller->value()-lastTwistValue));
    lastTwistValue=roller->value();

    if(state->interactor->isBetaStrand())
    {
        int numResidues=state->interactor->getStructure().getNumResidues();
        MD::Scalar* deltaPhis=new MD::Scalar[numResidues];
        MD::Scalar* deltaPsis=new MD::Scalar[numResidues];
        for(int i=0;i<numResidues;++i)
        {
            deltaPhis[i]=twistDelta;
            deltaPsis[i]=twistDelta;
        }
        state->interactor->changeDihedralAngles(deltaPhis,deltaPsis);
        delete[] deltaPhis;
        delete[] deltaPsis;

		state->proteinRenderer->updateProtein();
		if(state->energyCalculator!=0)
			state->energyCalculator->updateProtein();

        /* Update visualization state: */
        updateVisualization();
    }

    updateGui();
}

void structureDialogPleatCB(MyRoller* roller,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    static double lastPleatValue=0.0;
    MD::Scalar pleatDelta=MD::Scalar(Math::rad(roller->value()-lastPleatValue));
    lastPleatValue=roller->value();

    if(state->interactor->isBetaStrand())
    {
        int numResidues=state->interactor->getStructure().getNumResidues();
        MD::Scalar* deltaPhis=new MD::Scalar[numResidues];
        MD::Scalar* deltaPsis=new MD::Scalar[numResidues];
        for(int i=0;i<numResidues;++i)
        {
            deltaPhis[i]=-pleatDelta;
            deltaPsis[i]=pleatDelta;
        }
        state->interactor->changeDihedralAngles(deltaPhis,deltaPsis);
        delete[] deltaPhis;
        delete[] deltaPsis;

		state->proteinRenderer->updateProtein();
		if(state->energyCalculator!=0)
			state->energyCalculator->updateProtein();

        /* Update visualization state: */
        updateVisualization();
    }

    updateGui();
}

void structureDialogBraidCB(MyRoller* roller,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    static double lastBraidValue=0.0;
    MD::Scalar braidDelta=MD::Scalar(Math::rad(roller->value()-lastBraidValue));
    lastBraidValue=roller->value();

    if(state->interactor->isBetaStrand())
    {
        int numResidues=state->interactor->getStructure().getNumResidues();
        MD::Scalar* deltaPhis=new MD::Scalar[numResidues];
        MD::Scalar* deltaPsis=new MD::Scalar[numResidues];
        for(int i=0;i<numResidues;++i)
        {
            deltaPhis[i]=braidDelta;
            deltaPsis[i]=braidDelta;
            braidDelta=-braidDelta;
        }
        state->interactor->changeDihedralAngles(deltaPhis,deltaPsis);
        delete[] deltaPhis;
        delete[] deltaPsis;

		state->proteinRenderer->updateProtein();
		if(state->energyCalculator!=0)
			state->energyCalculator->updateProtein();

        /* Update visualization state: */
        updateVisualization();
    }

    updateGui();
}

void structureDialogFlattenBetaStrandCB(Fl_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isBetaStrand())
    {
        int numResidues=state->interactor->getStructure().getNumResidues();
        MD::Scalar* phis=new MD::Scalar[numResidues];
        MD::Scalar* psis=new MD::Scalar[numResidues];
        for(int i=0;i<numResidues;++i)
        {
            phis[i]=Math::rad(-131.85);
            psis[i]=Math::rad(128.15);
        }
        undoBuffer.startInteraction (
            state->interactor->getStructure(),
            state->interactor->getUpdateDirection()
        );
        state->interactor->setDihedralAngles(phis,psis);
        delete[] phis;
        delete[] psis;
        undoBuffer.finishInteraction();

		state->proteinRenderer->updateProtein();
		if(state->energyCalculator!=0)
			state->energyCalculator->updateProtein();


        /* Update visualization state: */
        updateVisualization();
    }

    updateGui();
}

void structureDialogResetBetaStrandCB(Fl_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
    {
        int numResidues=state->interactor->getStructure().getNumResidues();
        MD::Scalar* phis=new MD::Scalar[numResidues];
        MD::Scalar* psis=new MD::Scalar[numResidues];

        /* Reset dihedral angles inside the structure to the default values for that structure type: */
        switch(state->interactor->getStructure().getStructureType())
        {
            case MD::Protein::SecondaryStructure::COIL:
                for(int i=0;i<numResidues;++i)
                {
                    phis[i]=MD::Scalar(Math::rad(-182.0));
                    psis[i]=MD::Scalar(Math::rad(-182.0));
                }
                break;

            case MD::Protein::SecondaryStructure::ALPHA_HELIX:
                for(int i=0;i<numResidues;++i)
                {
                    phis[i]=MD::Scalar(Math::rad(-60.0));
                    psis[i]=MD::Scalar(Math::rad(-40.0));
                }
                break;

            case MD::Protein::SecondaryStructure::BETA_STRAND:
                for(int i=0;i<numResidues;++i)
                {
                    phis[i]=MD::Scalar(Math::rad(-120.0));
                    psis[i]=MD::Scalar(Math::rad(140.0));
                }
                break;

	    default:
	    	break;
        }

        undoBuffer.startInteraction (
            state->interactor->getStructure(),
            state->interactor->getUpdateDirection()
        );
        state->interactor->setDihedralAngles(phis,psis);
        delete[] phis;
        delete[] psis;
        undoBuffer.finishInteraction();

		state->proteinRenderer->updateProtein();
		if(state->energyCalculator!=0)
			state->energyCalculator->updateProtein();


        /* Update visualization state: */
        updateVisualization();
    }

    updateGui();
}

void renderingDialogDrawAtomsCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawAtoms(renderingDialogDrawAtomsToggle->value());
    checkDrawToggles(state);

    renderWindow->redraw();
    updateGui();
}

void renderingDialogDrawBondsCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawBonds(renderingDialogDrawBondsToggle->value());
    checkDrawToggles(state);

    renderWindow->redraw();
    updateGui();
}

void renderingDialogDrawBackboneCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawBackboneRibbon(renderingDialogDrawBackboneToggle->value());

    renderWindow->redraw();
    updateGui();
}

void renderingDialogDrawCartoonCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawCartoon(renderingDialogDrawCartoonToggle->value());
    checkDrawToggles(state);

    renderWindow->redraw();
    updateGui();
}

void renderingDialogDrawHydrogenBondsCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawHydrogenBonds(renderingDialogDrawHydrogenBondsToggle->value());
	checkShowToggles(state);

    renderWindow->redraw();
    updateGui();
}

void renderingDialogDrawHydrogenBondSitesCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawHydrogenBondSites (
        renderingDialogDrawHydrogenBondSitesToggle->value()
    );
	checkShowToggles(state);

    renderWindow->redraw();
    updateGui();
}

void renderingDialogDrawHydrogenCagesCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawHydrogenCages(renderingDialogDrawHydrogenCagesToggle->value());
	checkShowToggles(state);

    renderWindow->redraw();
    updateGui();
}
void renderingDialogDrawCPKCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawCPK(renderingDialogDrawCPKToggle->value());
	checkDrawToggles(state);
	
    renderWindow->redraw();
    updateGui();
}
void renderingDialogDrawTubeCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawTube(renderingDialogDrawTubeToggle->value());
	checkDrawToggles(state);

    renderWindow->redraw();
    updateGui();
}
void renderingDialogDrawLineCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawLine(renderingDialogDrawLineToggle->value());
	checkDrawToggles(state);

    renderWindow->redraw();
    updateGui();
}

void renderingDialogDrawCollisionsCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->proteinRenderer->setDrawCollisions(renderingDialogDrawCollisionsToggle->value());
	checkShowToggles(state);

    renderWindow->redraw();
    updateGui();
}

void renderingDialogStructureDrawAtomsCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawAtoms (
            state->interactor->getStructure(),
            renderingDialogStructureDrawAtomsToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}

void renderingDialogStructureDrawBondsCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawBonds (
            state->interactor->getStructure(),
            renderingDialogStructureDrawBondsToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}

void renderingDialogStructureDrawBackboneCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawBackboneRibbon (
            state->interactor->getStructure(),
            renderingDialogStructureDrawBackboneToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}

void renderingDialogStructureDrawCartoonCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawCartoon (
            state->interactor->getStructure(),
            renderingDialogStructureDrawCartoonToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}
void renderingDialogStructureDrawCPKCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawCPK (
            state->interactor->getStructure(),
            renderingDialogStructureDrawCPKToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}
void renderingDialogStructureDrawTubeCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawTube (
            state->interactor->getStructure(),
            renderingDialogStructureDrawTubeToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}
void renderingDialogStructureDrawLineCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawLine (
            state->interactor->getStructure(),
            renderingDialogStructureDrawLineToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}

void renderingDialogStructureDrawHydrogenBondSitesCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawHydrogenBondSites (
            state->interactor->getStructure(),
            renderingDialogStructureDrawHydrogenBondSitesToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}

void renderingDialogStructureDrawHydrogenCagesCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawHydrogenCages (
            state->interactor->getStructure(),
            renderingDialogStructureDrawHydrogenCagesToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}

void renderingDialogStructureDrawLargeHydrogenCagesCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(state->interactor->isValid())
        state->proteinRenderer->setDrawLargeHydrogenCages (
            state->interactor->getStructure(),
            renderingDialogStructureDrawLargeHydrogenCagesToggle->value()
        );

    renderWindow->redraw();
    updateGui();
}


void revisEnergy (ProteinState *state)
{
    // similar to updateVisualization() without some unnecessary calls
    // used by energy vis to resample and update display of energy cloud
    if ( !state ) return;
    if ( state->visualizeEnergy && state->energyRenderer )
        state->energyRenderer->sample();
    state->proteinRenderer->updateProtein();
    state->interactor->resetDragBox();
    ramachandranPlot->updateProtein();
    renderWindow->redraw();
    ramachandranPlot->redraw();
    updateGui();
}


void recordDialogCounterCB (Fl_Counter *counter, void *cbData)
{
    ProteinState *state = curProtein();
    if ( !state || !proteinRecord ) return;
    recordDialogSlider->value (counter->value());
    proteinRecord->load (uint(counter->value()));
    proteinRecord->apply (*state);
    revisEnergy (state);
}


void recordDialogSliderCB (Fl_Value_Slider *slider, void *cbData)
{
    ProteinState *state = curProtein();
    if ( !state || !proteinRecord ) return;
    recordDialogCounter->value (slider->value());
    proteinRecord->load (uint(slider->value()));
    proteinRecord->apply (*state);
    revisEnergy (state);
}


void recordDialogGenerateCB (Fl_Button *button, void *cbData)
{
    static const char *goLabel = "Generate Video Frames";
    static const char *stopLabel = "Stop Generating Video Frames";
    static const char *convertFmt = "convert frame.ppm frames/frame%05u.png";

    // check for the stop call first
    if ( !strcmp(stopLabel, button->label()) )
        button->label (goLabel);
    else
    {
        ProteinState *state = curProtein();
        if ( !state || !proteinRecord ) return;
        
        // put up the stop signal and initialize loop variables
        button->label (stopLabel);
        recordDialogCounter->deactivate();
        recordDialogSlider->deactivate();
        uint baseRecord = (uint) recordDialogSlider->value();
        uint frameStep = (uint) recordDialogFrameStepInput->value();
        uint frameNumber = 0;
        bool needFinalFrame = false;

        // generate frame 0 from the current state
        printf ("destroying old frames, if any\n");
        system ("rm -rf frames");
        system ("mkdir frames");
        printf ("exporting frame 0, record %u\n", baseRecord);
        renderWindow->saveScreenshot ("frame.ppm");
        snprintf (buf, buflen, convertFmt, frameNumber++);
        system (buf);

        // loop to generate all of the frames
        for ( uint record = baseRecord + 1;
              record < proteinRecord->numRecords();
              ++record)
        {
            // indicate our current progress in the dialog
            recordDialogCounter->value (record);
            recordDialogSlider->value (record);

            // load the next record
            proteinRecord->load (record);
            proteinRecord->apply (*state);
            if ( state->visualizeEnergy && state->energyRenderer )
                state->energyRenderer->updateEnergy();
            if ( (record - baseRecord) % frameStep == 0 )
            {
                // generate the next frame
                revisEnergy (state);
                Fl::flush();
                printf ("exporting frame %u, record %u\n", frameNumber, record);
                renderWindow->saveScreenshot ("frame.ppm");
                snprintf (buf, buflen, convertFmt, frameNumber++);
                system (buf);
                needFinalFrame = false;
            }
            else
                needFinalFrame = true;

            // check to see if it's time to stop
            Fl::check();
            if ( !strcmp(goLabel, button->label()) ||
                 !recordDialog->visible_r() )
            {
                // when it's time to stop, stop -- don't make another frame
                needFinalFrame = false;
                break;
            }
        }
        if ( needFinalFrame )
        {
            // generate that one last frame
            revisEnergy (state);
            Fl::flush();
            printf ("exporting frame %u, final\n", frameNumber);
            renderWindow->saveScreenshot ("frame.ppm");
            snprintf (buf, buflen, convertFmt, frameNumber++);
            system (buf);
        }
        system ("rm frame.ppm");
        recordDialogCounter->activate();
        recordDialogSlider->activate();
        if ( !strcmp(stopLabel, button->label()) )
            button->label (goLabel);
    }
}


void energyDialogRecalculateCB(Fl_Button* button,void* cbData)
{
    recalculateEnergy();
}


void waitForOptimizerCB (void *ptr)
{
    if ( !optState || !optState->energyCalculator ) return;

    Fl_Button *button = (Fl_Button*) ptr;
    if ( button && !strcmp(button->label(), stopLabel) )
    {
        // check for updates and schedule a new callback
        EnergyOptimizer *optimizer = optState->energyCalculator->optimizer();
        if ( optimizer )
        {
            if ( optimizer->isUpdateAvailable() )
            {
                if ( optState->visualizeEnergy && optState->energyRenderer )
                    optState->energyRenderer->updateEnergy();
                if ( (++updateNumber) % updateInterval == 0 )
                {
                    optState->protein->setAtomPositions (
                        optimizer->xCoordinates(),
                        optimizer->yCoordinates(),
                        optimizer->zCoordinates()
                    );
                    revisEnergy (optState);
                }
                optimizer->finishUpdate();
            }
            if ( optimizer->isRunning() )
                Fl::repeat_timeout (0.2, waitForOptimizerCB, button);
            else
            {
                // the optimizer actually finished
                optState = 0;
                button->label (runLabel);
                energyDialogLoadRecordButton->activate();
            }
        }
    }
}


void energyDialogOptimizeCB (Fl_Button *button, void *cbData)
{
    ProteinState *state = curProtein();
    if ( !state || !state->energyCalculator ) return;

    EnergyOptimizer *optimizer = state->energyCalculator->optimizer();
    if ( optimizer && !strcmp(button->label(), runLabel) )
    {
        optState = state;
        button->label (stopLabel);
        updateNumber = 0;
        if ( energyDialogRecordButton->value() )
            optimizer->record (
                configFile->retrieveString(
                    "/EnergyCalculator/recordFileName", "minimization.record"
                ).c_str()
            );
        energyDialogLoadRecordButton->deactivate();
        optimizer->optimize();
        Fl::add_timeout (0.2, waitForOptimizerCB, button);
    }
    else if ( optimizer )
    {
        optimizer->cancel();
        optimizer->record (0);
        optState = 0;
        button->label (runLabel);
        energyDialogLoadRecordButton->activate();
    }
}


void energyDialogVisualizeEnergyCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    if(energyDialogVisualizeEnergyToggle->value())
    {
        state->visualizeEnergy=true;
        state->proteinRenderer->setMapAtomValues(true);
        updateAtomEnergies();
    }
    else
    {
        state->visualizeEnergy=false;
        state->proteinRenderer->setMapAtomValues(false);
    }

    renderWindow->redraw();
    updateGui();
}


void energyUpdateRatInputCB(Fl_Value_Input* input,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
    state->engUpdateRate=energyUpdateRatInput->value();

    updateGui();
}


void energyDialogDisplayEnergyCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
     state->proteinRenderer->setEnergyUpdate(
	 energyDialogDisplayEnergyToggle->value());

    renderWindow->redraw();
    updateGui();
}


void energyDialogRangeMinCB(Fl_Value_Input* input,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->visualizeEnergyMinRange=energyDialogRangeMinInput->value();
    state->proteinRenderer->setMapAtomValueRange (
        state->visualizeEnergyMinRange,
        state->visualizeEnergyMaxRange
    );

    renderWindow->redraw();
    updateGui();
}


void energyDialogRangeMaxCB(Fl_Value_Input* input,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;

    state->visualizeEnergyMaxRange=energyDialogRangeMaxInput->value();
    state->proteinRenderer->setMapAtomValueRange (
        state->visualizeEnergyMinRange,
        state->visualizeEnergyMaxRange
    );

    renderWindow->redraw();
    updateGui();
}


void energyDialogComponentToggleCB(Fl_Check_Button* button,void* cbData)
{
    ProteinState *state = curProtein();
    if ( !state || !state->energyCalculator ) return;

    /* Get component index of this toggle: */
    int componentIndex=*static_cast<int*>(cbData);
    state->energyCalculator->setEnergyComponentState(componentIndex,button->value());

    if(state->visualizeEnergy)
        updateAtomEnergies();

    renderWindow->redraw();
    updateGui();
}


void energyDialogVolumeCB (Fl_Button *button, void *cbData)
{
    volumeDialog->show();
    updateGui();
}


void energyDialogLoadRecordCB (Fl_Button *button, void *cbData)
{
    ProteinState *state = curProtein();
    if ( !state || !energyLibrary ) return;
    char *recordFileName = fl_file_chooser (
        "Load Record File", "*.record", 0, 1
    );
    if ( !recordFileName ) return;
    delete proteinRecord;
    delete backupRecord;
    proteinRecord = new ProteinRecord (recordFileName);
    if ( proteinRecord->isEmpty() ||
         proteinRecord->numAtoms() != state->protein->getNumAtoms() ||
         proteinRecord->numEnergyComponents() !=
            energyLibrary->getNumEnergyComponents() )
    {
        fl_alert (
            "Either the record file '%s' was incompatible with the protein %s, "
            "or there was an I/O error when loading the file.",
            recordFileName, state->name
        );
        delete proteinRecord;
        proteinRecord = 0;
        return;
    }
    backupRecord = new ProteinRecord;
    backupRecord->snapshot (*state);
    proteinRecord->apply (*state);
    recordDialogCounter->maximum (proteinRecord->numRecords() - 1);
    recordDialogCounter->lstep (proteinRecord->numRecords() - 1);
    recordDialogSlider->maximum (proteinRecord->numRecords() - 1);
    snprintf (
        buf, buflen,
        "%s (%u records)",
        recordFileName, proteinRecord->numRecords()
    );
    recordDialogFileNameOutput->value (buf);
    energyDialogLoadRecordButton->deactivate();
    recordDialog->show();
    revisEnergy (state);
}


void volumeDialogSampleCB (Fl_Button *button, void *cbData)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer && state->energyCalculator )
    {
        // sample the energy texture
        state->energyRenderer->sample();
        if ( volumeDialogAutoNormalizeToggle->value() )
        {
            // store the newly generated normalizing range in the dialog
            volumeDialogRangeMinInput->value (
                state->energyRenderer->minNormal()
            );
            volumeDialogRangeMaxInput->value (
                state->energyRenderer->maxNormal()
            );
        }
        updateVisualization (false);
    }
}


void volumeDialogColorChoiceCB (Fl_Choice *choice, void *cbData)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
    {
        state->energyRenderer->setColorFunction (
            volumeDialogClassNumber, choice->value()
        );
        updateVisualization (false);
    }
}


void volumeDialogAutoNormalizeToggleCB (Fl_Check_Button *button, void *cbData)
{
    ProteinState *state = curProtein();
    if ( !state ) return;
    if ( button->value() )
    {
        volumeDialogRangeMinInput->deactivate();
        volumeDialogRangeMaxInput->deactivate();
        if ( state->energyRenderer )
            state->energyRenderer->setAutoNormalizing (true);
    }
    else
    {
        volumeDialogRangeMinInput->activate();
        volumeDialogRangeMaxInput->activate();
        if ( state->energyRenderer )
            state->energyRenderer->setAutoNormalizing (false);
    }
}


void volumeDialogClassifierChoiceCB (Fl_Choice *choice, void *cbData)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
    {
        volumeDialogClassNumber = 0;
        state->energyRenderer->setClassifier (choice->value());
        updateVolumeDialogGui();
    }
}


void volumeDialogNextClassButtonCB (Fl_Button *button, void *cbData)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
    {
        AtomClassifier *classifier = AtomClassifier::get (
            state->energyRenderer->classifier()
        );
        if ( classifier )
        {
            ++volumeDialogClassNumber;
            if ( volumeDialogClassNumber >= classifier->numClasses() )
                volumeDialogClassNumber = 0;
            updateVolumeDialogGui();
        }
    }
}


void volumeDialogPreviousClassButtonCB (Fl_Button *button, void *cbData)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
    {
        AtomClassifier *classifier = AtomClassifier::get (
            state->energyRenderer->classifier()
        );
        if ( classifier )
        {
            if ( volumeDialogClassNumber )
                --volumeDialogClassNumber;
            else
                volumeDialogClassNumber = classifier->numClasses() - 1;
            updateVolumeDialogGui();
        }
    }
}


void volumeDialogRadiusMultiplierCB (Fl_Value_Slider *slider, void *cbData)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
        state->energyRenderer->setRadiusMultiplier (slider->value());
}


void volumeDialogResolutionCB (Fl_Value_Slider *slider, void *cbData)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
        state->energyRenderer->setTexelsPerAngstrom (slider->value());
}


void volumeDialogChannelCB (Fl_Round_Button*, void*)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
        state->energyRenderer->setGradient (
            volumeDialogGradientRadioButton->value()
        );
}


void volumeDialogRadiusTypeCB (Fl_Round_Button*, void*)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
    {
        EnergyRenderer::RadiusType radiusType = EnergyRenderer::UNIFORM_RADIUS;
        if ( volumeDialogRadiusRadioButton->value() )
            radiusType = EnergyRenderer::ATOM_RADIUS;
        else if ( volumeDialogVDWRadioButton->value() )
            radiusType = EnergyRenderer::VAN_DER_WAALS_RADIUS;
        state->energyRenderer->setRadiusType (radiusType);
    }
}


void volumeDialogNormIntervalCB (Fl_Value_Input*, void*)
{
    ProteinState *state = curProtein();
    if ( state && state->energyRenderer )
        state->energyRenderer->setNormalRange (
            volumeDialogRangeMinInput->value(),
            volumeDialogRangeMaxInput->value()
        );
}


void selectionDialogXRollerCB (MyRoller *roller, void *cbData)
{
    typedef Geometry::AffineTransformation<double,3> Transform;

    static double lastValue = 0.0;
    double offset = roller->value() - lastValue;
    lastValue = roller->value();
    ProteinState *state = curProtein();
    if ( !state ) return;
    Transform t (
        Transform::translate (Geometry::Vector<double,3>(offset, 0.0, 0.0))
    );
    state->protein->transform (t);
    if ( state->energyRenderer )
        state->energyRenderer->transform (t);
    updateVisualization (false);
}


void selectionDialogYRollerCB (MyRoller *roller, void *cbData)
{
    typedef Geometry::AffineTransformation<double,3> Transform;

    static double lastValue = 0.0;
    double offset = roller->value() - lastValue;
    lastValue = roller->value();
    ProteinState *state = curProtein();
    if ( !state ) return;
    Transform t (
        Transform::translate (Geometry::Vector<double,3>(0.0, offset, 0.0))
    );
    state->protein->transform (t);
    if ( state->energyRenderer )
        state->energyRenderer->transform (t);
    updateVisualization (false);
}


void selectionDialogZRollerCB (MyRoller *roller, void *cbData)
{
    typedef Geometry::AffineTransformation<double,3> Transform;

    static double lastValue = 0.0;
    double offset = roller->value() - lastValue;
    lastValue = roller->value();
    ProteinState *state = curProtein();
    if ( !state ) return;
    Transform t (
        Transform::translate(Geometry::Vector<double,3>(0.0, 0.0, offset))
    );
    state->protein->transform (t);
    if ( state->energyRenderer )
        state->energyRenderer->transform (t);
    updateVisualization (false);
}


void selectionDialogXYRollerCB (MyRoller *roller, void *cbData)
{
    typedef Geometry::AffineTransformation<double,3> Transform;

    static double lastValue = 0.0;
    double offset = roller->value() - lastValue;
    lastValue = roller->value();
    ProteinState *state = curProtein();
    if ( !state ) return;
    Geometry::Vector<double,3> centroid (state->protein->calcCentroid());
    Transform t (Transform::translate(centroid));
    t *= Transform::rotate (Geometry::Rotation<double,3>::rotateZ(offset));
    t *= Transform::translate (-centroid);
    state->protein->transform (t);
    if ( state->energyRenderer )
        state->energyRenderer->transform (t);
    updateVisualization (false);
}


void selectionDialogYZRollerCB (MyRoller *roller, void *cbData)
{
    typedef Geometry::AffineTransformation<double,3> Transform;

    static double lastValue = 0.0;
    double offset = roller->value() - lastValue;
    lastValue = roller->value();
    ProteinState *state = curProtein();
    if ( !state ) return;
    Geometry::Vector<double,3> centroid (state->protein->calcCentroid());
    Transform t (Transform::translate(centroid));
    t *= Transform::rotate (Geometry::Rotation<double,3>::rotateX(offset));
    t *= Transform::translate (-centroid);
    state->protein->transform (t);
    if ( state->energyRenderer )
        state->energyRenderer->transform (t);
    updateVisualization (false);
}


void selectionDialogXZRollerCB (MyRoller *roller, void *cbData)
{
    typedef Geometry::AffineTransformation<double,3> Transform;

    static double lastValue = 0.0;
    double offset = roller->value() - lastValue;
    lastValue = roller->value();
    ProteinState *state = curProtein();
    if ( !state ) return;
    Geometry::Vector<double,3> centroid (state->protein->calcCentroid());
    Transform t (Transform::translate(centroid));
    t *= Transform::rotate (Geometry::Rotation<double,3>::rotateY(offset));
    t *= Transform::translate (-centroid);
    state->protein->transform (t);
    if ( state->energyRenderer )
        state->energyRenderer->transform (t);
    updateVisualization (false);
}


static void getSelectedProteins (ProteinState *&p1, ProteinState *&p2)
{
    p1 = p2 = 0;
    for ( int i = 1; i <= selectionDialogBrowser->size() && !p2; ++i )
    {
        if ( selectionDialogBrowser->selected(i) )
        {
            if ( !p1 )
                p1 = (ProteinState*) selectionDialogBrowser->data (i);
            else if ( !p2 )
                p2 = (ProteinState*) selectionDialogBrowser->data (i);
        }
    }
}


static void updateSelectionDialog (ProteinState *p1, ProteinState *p2)
{
    // called when the user selects a protein in the selection dialog
    if ( p1 )
    {
        // show the new selection to the rest of the application software
        selectProtein (p1);
        selectionDialogRemoveButton->activate();
		p1->proteinRenderer->grayOutCartoon(false);
    }
    else
        selectionDialogRemoveButton->deactivate();
    if ( p2 )
	{
        selectionDialogAlignButton->activate();
		p2->proteinRenderer->grayOutCartoon(true);
		renderWindow->setRMSDValue(-1);
    }
	else
        selectionDialogAlignButton->deactivate();
    for ( uint i = 1; i <= selectionDialogBrowser->size(); ++i )
    {
        // ensure that p1 and p2 are selected
        ProteinState *state = (ProteinState*) selectionDialogBrowser->data (i);
        if ( (p1 && state == p1) || (p2 && state == p2) )
            selectionDialogBrowser->select (i);
    }
}


void selectionDialogBrowserCB (Fl_Multi_Browser *browser, void *cbData)
{
    ProteinState *p1, *p2, *state = curProtein();
    getSelectedProteins (p1, p2);

    // force at least one protein to always be selected
    if ( state && !p1 && !p2 )
        p1 = state;
    updateSelectionDialog (p1, p2);
}


void selectionDialogRemoveButtonCB (Fl_Button *button, void *cbData)
{
    ProteinState *p1, *p2;
    getSelectedProteins (p1, p2);
    if ( p1 )
    {
        // give user a chance to save changes if necessary
        if ( !undoBuffer.isSaved() &&
             !fl_ask("Do you really want to discard all unsaved changes?") )
            return;

        // remove protein from selector and global data structure
        for ( uint i = 1; i <= selectionDialogBrowser->size(); ++i )
            if ( ((ProteinState*)selectionDialogBrowser->data(i)) == p1 )
                selectionDialogBrowser->remove (i);
        deleteProtein (p1);
    }
    // figure out which protein should be the new selection
    if ( p2 )
    {
        p1 = p2;
        p2 = 0;
    }
    else
        p1 = curProtein();
    updateSelectionDialog (p1, p2);
    if ( !p1 )
    {
        // if no protein remains, ensure that the UI is cleaned up
        undoBuffer.clear();
		hideAllDialogs();
        renderWindow->redraw();
        updateGui();
    }
}


void selectionDialogCenterViewButtonCB (Fl_Button *button, void *cbData)
{
    centerView();
}

void alignBaseAllCB (Fl_Round_Button *button, void *cbData)
{
alignBase = 0;
}
void alignBaseBackboneCB (Fl_Round_Button *button, void *cbData)
{
alignBase = 1;
}
void alignBaseCalphaCB (Fl_Round_Button *button, void *cbData)
{
alignBase = 2;
}
void selectionDialogAlignButtonCB (Fl_Button *button, void *cbData)
{
    using namespace MD;

    ProteinState *p1, *p2, *cur = curProtein();
    getSelectedProteins (p1, p2);
    if ( !p1 || !p2 )
    {
        // button should be disabled
        button->deactivate();
        return;
    }
    if ( cur != p1 )
    {
        // want current protein to re-use code in updateVisualization()
        p2 = p1;
        p1 = cur;
    }
    // decide how many atoms to use in the correspondence (this is primitive)
    integer numAtoms = p1->protein->getNumAtoms();
	if ( numAtoms > p2->protein->getNumAtoms() )
        numAtoms = p2->protein->getNumAtoms();

    // create and initialize variables to use in call to FITSQ routine
    doublereal *cloud1, *cloud2, *rotation, *translation, rootMeanSquare = 0.0;
    doublereal *cloud3, *cloud4;
    cloud1 = new doublereal[numAtoms * 3];
    cloud2 = new doublereal[numAtoms * 3];
    rotation = new doublereal[9];
    translation = new doublereal[3];
    if ( !cloud1 || !cloud2 )
    {
        delete[] cloud1;
        delete[] cloud2;
        fl_alert ("Insufficient memory to complete operation.");
        return;
    }
    memset (rotation, 0, sizeof(doublereal)*9);
    memset (translation, 0, sizeof(doublereal)*3);
    const Protein::ChainAtom* const *a1 = p1->protein->getAtomPointers();
    const Protein::ChainAtom* const *a2 = p2->protein->getAtomPointers();
    uint base1, base2 =0;
	integer actuallyNumber =0;
	
	int caAtomNum, caAtomNum1 =0, caAtomNum2 =0;
	for ( uint i = 0; i < numAtoms*3; ++i )
	{
	//	if(a1[i]->getAtomName())
	//	if(strcmp(a1[i]->getAtomName(),"CA")==0)
	//		caAtomNum1++;
		
	//	if(a2[i]->getAtomName())
	//	if(strcmp(a2[i]->getAtomName(),"CA")==0)
	//		caAtomNum2++;
	}
    
	caAtomNum = caAtomNum1;
	if ( caAtomNum1 > caAtomNum2 )
        caAtomNum = caAtomNum2;
	
	
	
	for ( uint i = 0; i < numAtoms*3; ++i )
	{
	 cloud1[i]=0; cloud2[i]=0;
	}
	switch (alignBase)
	{
		case 0: // all atoms
		{
			for ( uint i = 0, base =0; i < numAtoms; ++i, base+=3 )
    		{
        		assert ( a1[i] && a2[i] );
        		Position pos1 (a1[i]->getPosition());
				Position pos2 (a2[i]->getPosition());
				for ( uint j = 0; j < 3; ++j )
        		{
            		// load the atom coordinates into the arrays
            		cloud1[base+j] = pos1[j];
            		cloud2[base+j] = pos2[j];
        		}
    		}
			actuallyNumber = numAtoms;
			break;
		}
		case 1:	// BackBone atoms only
		{
			const char* BBType[3]={"CA","N","C"};
			int tag =-1;
			int index;
			for ( uint i = 0, base =0; i < numAtoms; ++i, base+=3 )
    		{
				assert ( a1[i] );
				for(index=0;index<3;++index) 
					if(strcmp(BBType[index],a1[i]->getAtomName())==0)
            			break;
        		
				// not a backbone atom
				if (index>2)
					continue;
	
				if(strcmp(BBType[index],a1[i]->getAtomName())==0)
				{
					for ( uint j = ++tag, base =0; j < numAtoms; ++j, base+=3 )
    				{
						assert ( a2[j] );
        				if(strcmp(BBType[index],a2[j]->getAtomName())==0)
        				{
							Position pos1 (a1[i]->getPosition());
							Position pos2 (a2[j]->getPosition());
							for ( uint k = 0; k < 3; ++k )
        					{
            					// load the atom coordinates into the arrays
            					cloud1[base+k] = pos1[k];
            					cloud2[base+k] = pos2[k];
        					}
							actuallyNumber++;
							tag =j;
							break;
						}
					}
				}
	   		}
	   		break;
	   }
	   case 2:	// C-Aplha atoms only
	   {
			cloud3 = new doublereal[caAtomNum * 3];
    		cloud4 = new doublereal[caAtomNum * 3];

			int tag = -1;
			for ( uint i = 0, base =0; i < numAtoms; ++i, base+=3 )
    		{
				assert ( a1[i] );
				if(strcmp(a1[i]->getAtomName(),"CA")==0)
				{
					for ( uint j = ++tag, base =0; j < numAtoms; ++j, base+=3 )
    				{
        				assert ( a2[j] );
						if(strcmp(a2[j]->getAtomName(),"CA")==0)
        				{
							Position pos1 (a1[i]->getPosition());
							Position pos2 (a2[j]->getPosition());
							for ( uint k = 0; k < 3; ++k )
        					{
            					// load the atom coordinates into the arrays
            					//cloud1[base+k] = pos1[k];
            					//cloud2[base+k] = pos2[k];
            					cloud3[base+k] = pos1[k];
            					cloud4[base+k] = pos2[k];
        					}
							actuallyNumber++;
							tag =j;
							break;

						}
					}
				}
		   }
		   break;
	   }
	}
    sprintf (msg.buf, "Aligning proteins based on %d atoms.",actuallyNumber) ;
    msg.Info (RUNID, msg.buf);
	for ( uint i = 0; i < numAtoms*3; ++i )
	{
	 printf("%f %f %f\n", cloud1[i][0], cloud1[i][1], cloud1[i][2];
	}

	for ( uint i = 0; i < numAtoms*3; ++i )
	{
	 printf("%f %f %f\n", cloud2[i][0], cloud2[i][1], cloud2[i][2];
	}

    // call FITSQ and report RMSD
    fitsq_ (&rootMeanSquare, cloud1, cloud2, &actuallyNumber, translation, rotation);
//    fitsq_ (&rootMeanSquare, cloud3, cloud4, &actuallyNumber, translation, rotation);

    delete[] cloud1;
    delete[] cloud2;
	if ( cloud3 || cloud4 )
    {
        delete[] cloud3;
        delete[] cloud4;
	}
    sprintf (msg.buf, "Root mean square distance = %f.",sqrt(fabs(rootMeanSquare)));
    msg.Info (RUNID, msg.buf);
	renderWindow->setRMSDValue(sqrt(fabs(rootMeanSquare)));

    // build a transformation and apply it to the current protein
    Geometry::Matrix<Scalar,4,4> matrix;
    for ( uint row = 0; row < 3; ++row )
        for ( uint col = 0; col < 3; ++col )
            matrix(row,col) = 0.0;
    
	for ( uint row = 0; row < 3; ++row )
    {
        for ( uint col = 0; col < 3; ++col )
            matrix(row,col) = rotation[col*3 + row];
        matrix(row,3) = translation[row];
        matrix(3,row) = 0.0;
    }
    matrix(3,3) = 1.0;
    p1->protein->transform (Geometry::Transformation<Scalar,3>(matrix));
    updateVisualization (false);
	delete [] rotation;
	delete [] translation;
}


void aboutDialogOkCB(Fl_Return_Button* button,void* cbData)
{
    aboutDialog->hide();
}

void applyTransformation(const ProteinInteractor::Transformation& transformation)
{
    ProteinState *state = curProtein();
    if ( !state ) {
       printf("applyTransormation returning early because of no state\n");
       return;
       }

    typedef ProteinInteractor::Transformation Transformation;
    typedef Transformation::Rotation Rotation;
    typedef Transformation::Point Point;
    typedef Transformation::Scalar Scalar;
    typedef Transformation::Vector Vector;

    /* Decompose the given transformation: */
    Vector translation=transformation.getTranslation();
    Vector rotationAxis=transformation.getRotation().getAxis();
    Scalar rotationAngle=transformation.getRotation().getAngle();
    Point rotateCenterStart=state->interactor->getBoxRotateCenter();
    Point rotateCenterEnd=transformation.transform(rotateCenterStart);

    /* Calculate partial transformations to go in several small steps: */
    int numSteps=int(Math::ceil(Geometry::dist(rotateCenterStart,rotateCenterEnd)/2.0)); // At most 5 Angstrom per step
    int numRotateSteps=int(Math::ceil(Math::abs(rotationAngle)/0.2)); // At most 0.5 radians per step
    if(numSteps<numRotateSteps)
        numSteps=numRotateSteps;
    if(state->interactor->startInteraction())
    {
        for(int i=0;i<numSteps;++i)
        {
            Scalar s=Scalar(i+1)/Scalar(numSteps);
            Rotation r=Rotation::rotateAxis(rotationAxis,rotationAngle*s);
            Point rc=Geometry::affineCombination(rotateCenterStart,rotateCenterEnd,s);
            Vector t=rc-r.transform(rotateCenterStart);
            Transformation stepTransformation(t,r);
            state->interactor->drag(stepTransformation);
            state->interactor->applyChanges();
            if(offlineBuildBeta == 0) updateProteinNow();
        }
        state->interactor->finishInteraction();
        if(offlineBuildBeta == 0) updateDihedralAngles();
        updateGui();
    }
	//msg.Warn(PDBID, "Not ready for interaction!");
}

int zipAntiParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float pl[3], float cn[3])
{
    typedef ProteinInteractor::Transformation Transformation;
    float test, dentest;

    /* Zip the two selected residues together: */
    MD::Protein::Dipole amide1=manipulatedResidue->getAmide();
    MD::Protein::Dipole carboxyl1=manipulatedResidue->getCarboxyl();
    MD::Protein::Dipole amide2=anchorResidue->getAmide();
    MD::Protein::Dipole carboxyl2=anchorResidue->getCarboxyl();
    if(amide1.isValid()&&carboxyl1.isValid()&&amide2.isValid()&&carboxyl2.isValid())
    {
        Transformation goalTransformation=Transformation::identity;

        MD::Point a1=amide1.getMajorAtom()->getPosition();
        MD::Point a2=amide1.getBondSite();
        MD::Point c1=carboxyl1.getMajorAtom()->getPosition();
        MD::Point c2=carboxyl1.getBondSite();

        /* Align the two residues' average bonding site directions: */
        MD::Vector d1=Geometry::mid(amide2.getBondSite(),carboxyl2.getBondSite())-Geometry::mid(amide2.getMajorAtom()->getPosition(),carboxyl2.getMajorAtom()->getPosition());
        MD::Vector d2=Geometry::mid(a2,c2)-Geometry::mid(a1,c1);
        MD::Vector axis=Geometry::cross(d1,d2);
	dentest = (d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]) *
	      (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]);
	if (dentest != 0.) {
	   test = (axis[0]*axis[0] + axis[1]*axis[1] +
	      axis[2]*axis[2]) / dentest;
	   if(test > .00000001) {
              MD::Scalar angle1=Math::acos(-(d1*d2)/Math::sqrt(Geometry::sqr(d1)*Geometry::sqr(d2)));
              Transformation trans1=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,angle1));
              goalTransformation.leftMultiply(trans1);
      
              a1=trans1.transform(a1);
              a2=trans1.transform(a2);
              c1=trans1.transform(c1);
              c2=trans1.transform(c2);
              d2=trans1.transform(d2);
	      }
	   }

        /* Align the two planes formed by the residues' bonding sites: */
        MD::Vector p1=Geometry::cross(d1,carboxyl2.getBondSite()-amide2.getBondSite());
	MD::Vector plane = normalize(p1);
	MD::Point center = Geometry::mid(carboxyl2.getBondSite(), amide2.getBondSite());
	for(int i = 0; i < 3; ++i) {
	   pl[i] = - plane[i];
	   cn[i] = - center[i];
	   }
	MD::Vector e2=c2-a2;
	MD::Vector p2=Geometry::cross(d2,e2);
        dentest = (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]) *
           (e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
	if(dentest != 0.) {
	   test = (p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]) / dentest;
	   if(test > .00000001) {
              MD::Scalar angle2=Math::acos((p1*p2)/Math::sqrt(Geometry::sqr(p1)*Geometry::sqr(p2)));
              if(Geometry::cross(p2,p1)*d2<MD::Scalar(0))
                  angle2=-angle2;
              Transformation trans2=Transformation::rotate(Transformation::Rotation::rotateAxis(d2,angle2));
              goalTransformation.leftMultiply(trans2);
      
              a2=trans2.transform(a2);
              c2=trans2.transform(c2);
	      }
	   }

        /* Align the midpoints of the two residues' bonding sites: */
        MD::Vector dist=Geometry::mid(amide2.getBondSite(),carboxyl2.getBondSite())-Geometry::mid(a2,c2);
        Transformation trans3=Transformation::translate(dist);
        goalTransformation.leftMultiply(trans3);

        applyTransformation(goalTransformation);
        // renderWindow->zipIt(goalTransformation);
        return(1);
    }
		else {
		   if(0) printf(": input %d %d %d %d is not valid\n",
		      amide1.isValid(), carboxyl1.isValid(), amide2.isValid(),
		         carboxyl2.isValid());
		   return(0);
		   }
    }

int zipParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float pl[3], float cn[3])
{
    typedef ProteinInteractor::Transformation Transformation;
    float test, dentest;

    /* Zip the two selected residues together: */
    MD::Protein::Dipole amide1=manipulatedResidue->getAmide();
    MD::Protein::Dipole carboxyl1=manipulatedResidue->getCarboxyl();
    MD::Protein::Residue* anchor1=anchorResidue->getPred();
    MD::Protein::Residue* anchor2=anchorResidue->getSucc();
    if(anchor1!=0&&anchor2!=0)
    {
        MD::Protein::Dipole amide2=anchor2->getAmide();
        MD::Protein::Dipole carboxyl2=anchor1->getCarboxyl();
        if(amide1.isValid()&&carboxyl1.isValid()&&amide2.isValid()&&carboxyl2.isValid())
        {
            Transformation goalTransformation=Transformation::identity;

            MD::Point a1=amide1.getBondSite();
            MD::Point c1=carboxyl1.getBondSite();
            MD::Point p1=Geometry::mid(a1,c1);
            MD::Vector d1=c1-a1;
            MD::Vector s1=(a1-amide1.getMajorAtom()->getPosition())+(c1-carboxyl1.getMajorAtom()->getPosition());
            MD::Point a2=amide2.getBondSite();
            MD::Point c2=carboxyl2.getBondSite();
            MD::Point p2=Geometry::mid(a2,c2);
            MD::Vector d2=a2-c2;
            MD::Vector s2=(a2-amide2.getMajorAtom()->getPosition())+(c2-carboxyl2.getMajorAtom()->getPosition());

            /* Align the two residues' average bonding site directions: */
            MD::Vector axis=Geometry::cross(d1,d2);
	    dentest = (d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]) *
	          (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]);
	    if (dentest != 0.) {
	       test = (axis[0]*axis[0] + axis[1]*axis[1] +
	          axis[2]*axis[2]) / dentest;
	       if(test > .00000001) {
                  MD::Scalar angle1=Math::acos((d1*d2)/Math::sqrt(Geometry::sqr(d1)*Geometry::sqr(d2)));
                  Transformation trans1=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,angle1));
                  goalTransformation.leftMultiply(trans1);

                  p1=trans1.transform(p1);
                  d1=trans1.transform(d1);
                  s1=trans1.transform(s1);
		  }
	       }

            /* Align the two planes formed by the residues' bonding sites: */
            MD::Vector pl1=Geometry::cross(d1,s1);
            MD::Vector pl2=Geometry::cross(s2,d2);
            MD::Vector axis2=Geometry::cross(pl1,pl2);
	    MD::Vector plane = normalize(pl2);
	    MD::Point center = p2;
	    for(int i = 0; i < 3; ++i) {
	       pl[i] = -plane[i];
	       cn[i] = -center[i];
	       }
            dentest = (pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]) *
                   (pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
	    if (dentest != 0.) {
	       test = (axis2[0]*axis2[0] + axis2[1]*axis2[1] +
	          axis2[2]*axis2[2]) / dentest;
	       if(test > .00000001) {
                  MD::Scalar angle2=Math::acos((pl1*pl2)/Math::sqrt(Geometry::sqr(pl1)*Geometry::sqr(pl2)));
                  Transformation trans2=Transformation::rotate(Transformation::Rotation::rotateAxis(axis2,angle2));
                  goalTransformation.leftMultiply(trans2);
      
                  p1=trans2.transform(p1);
		  }
	       }

            /* Align the midpoints of the two residues' bonding sites: */
            Transformation trans3=Transformation::translate(p2-p1);
            goalTransformation.leftMultiply(trans3);

            applyTransformation(goalTransformation);
            // renderWindow->zipIt(goalTransformation);
	    return(1);
        }
		else {
		   if(0) printf(": input %d %d %d %d is not valid\n",
		      amide1.isValid(), carboxyl1.isValid(), amide2.isValid(),
		         carboxyl2.isValid());
		   return(0);
		   }
                }
    else return(0);
    }

int altZipAntiParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float pl[3], float cn[3])
	{
	typedef ProteinInteractor::Transformation Transformation;
	float test, dentest;

	/* Zip the two selected residues together: */
	MD::Protein::Residue* anchor1=anchorResidue->getPred();
	MD::Protein::Residue* anchor2=anchorResidue->getSucc();
	MD::Protein::Residue* manip1=manipulatedResidue->getSucc();
	MD::Protein::Residue* manip2=manipulatedResidue->getPred();
	if(anchor1 == 0 || anchor2 == 0 || manip1 == 0 || manip2 == 0) return(0);
	MD::Protein::Dipole amide1=manip1->getAmide();
	MD::Protein::Dipole carboxyl1=manip2->getCarboxyl();
	MD::Protein::Dipole amide2=anchor2->getAmide();
	MD::Protein::Dipole carboxyl2=anchor1->getCarboxyl();
	if(amide1.isValid()&&carboxyl1.isValid()&&amide2.isValid()&&carboxyl2.isValid())
		{
		Transformation goalTransformation=Transformation::identity;

		MD::Point a1=amide1.getMajorAtom()->getPosition();
		MD::Point a2=amide1.getBondSite();
		MD::Point c1=carboxyl1.getMajorAtom()->getPosition();
		MD::Point c2=carboxyl1.getBondSite();

		/* Align the two residues' average bonding site directions: */
		MD::Vector d1=Geometry::mid(amide2.getBondSite(),carboxyl2.getBondSite())-Geometry::mid(amide2.getMajorAtom()->getPosition(),carboxyl2.getMajorAtom()->getPosition());
		MD::Vector d2=Geometry::mid(a2,c2)-Geometry::mid(a1,c1);
		MD::Vector axis=Geometry::cross(d1,d2);
	        dentest = (d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]) *
	              (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]);
		if (dentest != 0.) {
		   test = (axis[0]*axis[0] + axis[1]*axis[1] +
		      axis[2]*axis[2]) / dentest;
		   if(test > .00000001) {
		      MD::Scalar angle1=Math::acos(-(d1*d2)/Math::sqrt(Geometry::sqr(d1)*Geometry::sqr(d2)));
		      Transformation trans1=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,angle1));
		      goalTransformation.leftMultiply(trans1);

		      a1=trans1.transform(a1);
		      a2=trans1.transform(a2);
		      c1=trans1.transform(c1);
		      c2=trans1.transform(c2);
		      d2=trans1.transform(d2);
		      }
		   }

		/* Align the two planes formed by the residues' bonding sites: */
		MD::Vector e1 = carboxyl2.getBondSite()-amide2.getBondSite();
		MD::Vector p1=Geometry::cross(d1,carboxyl2.getBondSite()-amide2.getBondSite());
	        MD::Vector plane = normalize(p1);
	        MD::Point center = Geometry::mid(carboxyl2.getBondSite(), amide2.getBondSite());
		for(int i = 0; i < 3; ++i) {
		   pl[i] = plane[i];
		   cn[i] = center[i];
		   }
		MD::Vector e2=c2-a2;
		MD::Vector p2=Geometry::cross(d2,e2);
	        dentest = (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]) *
	           (e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
		if(dentest != 0.) {
		   test = (p2[0]*p2[0] + p2[1]*p2[1] +
		      p2[2]*p2[2]) / dentest;
		   if(dentest > .00000001) {
		      MD::Scalar angle2=Math::acos((p1*p2)/Math::sqrt(Geometry::sqr(p1)*Geometry::sqr(p2)));
		      if(Geometry::cross(p2,p1)*d2<MD::Scalar(0))
			      angle2=-angle2;
		      Transformation trans2=Transformation::rotate(Transformation::Rotation::rotateAxis(d2,angle2));
		      goalTransformation.leftMultiply(trans2);

		      a2=trans2.transform(a2);
		      c2=trans2.transform(c2);
		      }
		   }

		/* Align the midpoints of the two residues' bonding sites: */
		MD::Vector dist=Geometry::mid(amide2.getBondSite(),carboxyl2.getBondSite())-Geometry::mid(a2,c2);
		Transformation trans3=Transformation::translate(dist);
		goalTransformation.leftMultiply(trans3);

		applyTransformation(goalTransformation);
		// renderWindow->zipIt(goalTransformation);
		return(1);
		}
		else {
		   if(0) printf(": input %d %d %d %d is not valid\n",
		      amide1.isValid(), carboxyl1.isValid(), amide2.isValid(),
		         carboxyl2.isValid());
		   return(0);
		   }
    }

int altZipParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float pl[3], float cn[3])
	{
	typedef ProteinInteractor::Transformation Transformation;
	float test, dentest;

	/* Zip the two selected residues together: */
	MD::Protein::Dipole amide2=anchorResidue->getAmide();
	MD::Protein::Dipole carboxyl2=anchorResidue->getCarboxyl();
	MD::Protein::Residue* manip1=manipulatedResidue->getPred();
	MD::Protein::Residue* manip2=manipulatedResidue->getSucc();
	if(manip1!=0&&manip2!=0)
		{
		MD::Protein::Dipole amide1=manip2->getAmide();
		MD::Protein::Dipole carboxyl1=manip1->getCarboxyl();
		if(amide1.isValid()&&carboxyl1.isValid()&&amide2.isValid()&&carboxyl2.isValid())
			{
			Transformation goalTransformation=Transformation::identity;

			MD::Point a1=amide1.getBondSite();
			MD::Point c1=carboxyl1.getBondSite();
			MD::Point p1=Geometry::mid(a1,c1);
			MD::Vector d1=c1-a1;
			MD::Vector s1=(a1-amide1.getMajorAtom()->getPosition())+(c1-carboxyl1.getMajorAtom()->getPosition());
			MD::Point a2=amide2.getBondSite();
			MD::Point c2=carboxyl2.getBondSite();
			MD::Point p2=Geometry::mid(a2,c2);
			MD::Vector d2=a2-c2;
			MD::Vector s2=(a2-amide2.getMajorAtom()->getPosition())+(c2-carboxyl2.getMajorAtom()->getPosition());

			/* Align the two residues' average bonding site directions: */
			MD::Vector axis=Geometry::cross(d1,d2);

	                dentest = (d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]) *
	                   (d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2]);
			if (dentest != 0.) {
			   test = (axis[0]*axis[0] + axis[1]*axis[1] +
			      axis[2]*axis[2]) / dentest;
			   if(test > .00000001) {
			      MD::Scalar angle1=Math::acos((d1*d2)/Math::sqrt(Geometry::sqr(d1)*Geometry::sqr(d2)));
			      Transformation trans1=Transformation::rotate(Transformation::Rotation::rotateAxis(axis,angle1));
			      goalTransformation.leftMultiply(trans1);
			      p1=trans1.transform(p1);
			      d1=trans1.transform(d1);
			      s1=trans1.transform(s1);
			      }
			   }

			/* Align the two planes formed by the residues' bonding sites: */
			MD::Vector pl1=Geometry::cross(d1,s1);
			MD::Vector pl2=Geometry::cross(s2,d2);
			MD::Vector axis2=Geometry::cross(pl1,pl2);
	                MD::Vector plane = normalize(pl2);
	                MD::Point center = p2;
			for(int i = 0; i < 3; ++i) {
			   pl[i] = plane[i];
			   cn[i] = center[i];
			   }
	                dentest = (pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]) *
	                   (pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
			if (dentest != 0.) {
			   test = (axis2[0]*axis2[0] + axis2[1]*axis2[1] +
			      axis2[2]*axis2[2]) / dentest;
			   if(test > .00000001) {
			      MD::Scalar angle2=Math::acos((pl1*pl2)/Math::sqrt(Geometry::sqr(pl1)*Geometry::sqr(pl2)));
			      Transformation trans2=Transformation::rotate(Transformation::Rotation::rotateAxis(axis2,angle2));
			      goalTransformation.leftMultiply(trans2);

			      p1=trans2.transform(p1);
			      }
			   }

			/* Align the midpoints of the two residues' bonding sites: */
			Transformation trans3=Transformation::translate(p2-p1);
			goalTransformation.leftMultiply(trans3);

			applyTransformation(goalTransformation);
			// renderWindow->zipIt(goalTransformation);
                        return(1);
                }
                else {
                   if(0) printf(": input %d %d %d %d is not valid\n",
                      amide1.isValid(), carboxyl1.isValid(), amide2.isValid(),
                         carboxyl2.isValid());
                   return(0);
                   }
                }
        else return(0);
        }

void* ikUpdateThread(void*)
{
    ikUpdateDone=true;
    pthread_mutex_lock(&proteinMutex);
    while(true)
    {
        /* Wait for the next update request: */
        if(!ikUpdateRequested)
        {
            ikUpdateDone=true;
            pthread_cond_signal(&ikUpdateDoneCond);
            pthread_cond_wait(&ikUpdateRequestedCond,&proteinMutex);
        }
        ProteinState *state = curProtein();
        if ( state )
        {
            /* Grab the current goal transformation from the drag box: */
            ProteinInteractor::Transformation goalTransformation =
                state->interactor->getDragTransformation();

            // ikUpdateRequested=false; // Keep going even when user keeps mouse still
            ikUpdateDone=false;
            pthread_mutex_unlock(&proteinMutex);

            /* Perform IK steps: */
            state->interactor->drag(goalTransformation);

            /* Apply changes to protein: */
            pthread_mutex_lock(&proteinMutex);
            state->interactor->applyChanges();
            state->proteinRenderer->updateProtein();
            if(state->energyCalculator!=0)
                state->energyCalculator->updateProtein();
            ramachandranPlot->updateProtein();
            updateDihedralAngles();

            /* Notify the main thread to redraw its window: */
            ikUpdatePosted=true;
        }
    }
}


int main (int argc, char **argv)
{
    using namespace MD;
    try
    {
        /* Open configuration file: */
        configFile = new ConfigurationFile ("ProteinShop.cfg");
    }
    catch ( std::runtime_error error )
    {
        msg.Error (CONFIGID,"",error.what());
        fl_alert ("Could not load configuration file ProteinShop.cfg");
        return 1;
    }

    // Introduction Message
	printf("\n");
	msg.Info (INITID,"=============================================================");
	msg.Info (INITID, "ProteinShop version 3.0 Copyright (c) 2005");
	msg.Info (INITID, "The Regents of the University of California.");
	msg.Info (INITID, "All Rights Reserved.	http://proteinshop.lbl.gov");
	msg.Info (INITID, "");
	msg.Info (INITID, "	ProteinShop is open-source software.");
	msg.Info (INITID, "	ProteinShop comes with ABSOLUTELY NO WARRANTY;");
	msg.Info (INITID, "	for details check the LICENSE file or hit H.");
	msg.Info (INITID, "");
	msg.Info (INITID, "	Please cite ProteinShop in published work:");
	msg.Info (INITID, "	S. Crivelli, O. Kreylos, B. Hamann, N. Max, and W. Bethel.");
	msg.Info (INITID, "	ProteinShop: A Tool for Interactive Protein Manipulation.");
	msg.Info (INITID, "	J of Computer-Aided Molecular Design, 2004, 18, 271-285.");
	msg.Info (INITID, "");
	msg.Info (INITID,"=============================================================");
	msg.Info (INITID, "");
	msg.Info (INITID, "");



    /* Initialize IK update thread: */
    pthread_mutex_init(&proteinMutex,0);
    pthread_cond_init(&ikUpdateRequestedCond,0);
    pthread_cond_init(&ikUpdateDoneCond,0);
    pthread_create(&ikUpdateThreadId,0,ikUpdateThread,0);

    /* Create the GUI widget tree: */
    makeWindow();
   	
    /* Initialize energy calculation module: */
    try
    {
        /* Try loading the energy calculator DSO given in the configuration file: */
        loadEnergyLibrary(configFile->retrieveString("/EnergyCalculator/dsoName").c_str());
        updateInterval = retrieveValue (
            *configFile, "/EnergyCalculator/updateInterval", 20
        );
        snprintf (buf, buflen, "loaded update interval = %d", updateInterval);
        msg.Info (CONFIGID, buf);
    }
    catch(TagNotFoundError error)
    {
        /* Silently disable the energy calculation module: */
        energyLibrary=0;
    }
    Protein *firstProtein = 0;
    if ( argc > 1 )
    {
        // load protein files from the command line
	int extra = 0;
        for ( int i = 1; i < argc; ++i )
        {
	  printf("%s   previous extra %d\n", argv[i], extra);
	  if(!strcasecmp(argv[i], "extra") )
	    extra = 1;
	  else
	  {
            const char *ext = strrchr (argv[i], '.');
            Protein *protein = 0;
            if ( ext && !strcasecmp(ext, ".pdb") )
			{
				const char* chainList = parseChains(argv[i]);
				if (strlen(chainList) > 1 || parseModels(argv[i]) > 0)
					loadProteins (argv[i]);
				else
                	protein = loadPdb (argv[i]);
            }
			else if ( ext && !strcasecmp(ext, ".pred") )
            {
			    protein = loadPrediction (argv[i],1 - extra);
				extra = 0;
			}
            else
                fl_message ("Unrecognized file name extension:  %s\n"
                            "while attempting to load %s.", ext, argv[i]);
            if ( protein )
            {
				addProtein (protein, selectionDialogBrowser->size());
                if ( !firstProtein )
                    firstProtein = protein;
            }
          }
        }
    }
    if ( !firstProtein ) inputfilename[0] = 0;
    else
    {
        initBuild (firstProtein);
	printf("returned from command line initiated initBuild\n");
        if ( offlineBuildBeta == 1 )
        {
           BuildBeta();
           exit(-1);
        }
    }

    /* Open the main window and run: */
    mainWindow->show();
    updateGui();
    Fl::run();

    /* Clean up and exit: */
    deleteAllProteins();
    delete configFile;
    delete energyLibrary;
    return 0;
}
