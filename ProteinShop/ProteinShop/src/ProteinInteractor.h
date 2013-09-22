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
ProteinInteractor - Class encapsulating the manipulation state of a
protein.
***********************************************************************/


class ProteinInteractor;


#ifndef PROTEININTERACTOR_INCLUDED
#define PROTEININTERACTOR_INCLUDED

#ifdef __GNUC__
/***************************************************************
This must be some bug in the STL implementation coming with g++:
operator== for lists is declared as a friend function, but not
using the new(?) C++-syntax for template friend functions.
**************************************************************/

namespace std {
template <class _Tp, class _Alloc>
class list;
template <class _Tp, class _Alloc>
bool operator==(const list<_Tp,_Alloc>& __x,const list<_Tp,_Alloc>& __y);
}
#endif
#include <list>
#include <utility>

#include "DragBox.h"
#include "Protein.h"
#include "ProteinRenderer.h"
#include "UndoBuffer.h"
#include "IK.h"
#include "EnergyAPI.h"

class ProteinInteractor
{
    /* Embedded classes: */
    public:
    typedef IK::Transformation Transformation; // Type for transformations
    typedef std::list<MD::Protein::StructureSelector> StructureList; // Type for list of secondary structures
    typedef std::pair<int,MD::Protein::BackboneIterator::BondType> JointMapItem;
    
    enum UpdateDirection // Enumerated type for update directions
    {
        LEFT,RIGHT
    };
    
    private:
    class CoilRegionAngles // Class to hold dihedral angles for a coil region
    {
        /* Elements: */
        public:
        MD::Protein::StructureSelector structure; // Selector for the coil region
        MD::Scalar* phis; // Pointer to C-style array of phi angles
        MD::Scalar* psis; // Pointer to C-style array of psi angles
        
        /* Constructors and destructors: */
        CoilRegionAngles(const MD::Protein::StructureSelector& sStructure) // Creates angle arrays for secondary structure
            :structure(sStructure),
             phis(new MD::Scalar[structure.getNumResidues()]),
             psis(new MD::Scalar[structure.getNumResidues()])
        {
        };
        ~CoilRegionAngles(void)
        {
            delete[] phis;
            delete[] psis;
        };
        
        /* Methods: */
        int getFirstResidueIndex(void) const
        {
            return structure.getFirstResidueIndex();
        };
        int getNumResidues(void) const
        {
            return structure.getNumResidues();
        };
        void getCurrentAngles(void) // Retrieves current dihedral angles from structure
        {
            structure.getDihedralAngles(phis,psis);
        };
        void changeDihedralAngles(UpdateDirection updateDirection)
        {
            structure.changeDihedralAngles(phis,psis,updateDirection==LEFT?-1:1);
        };
    };
    
    /* Elements: */
    private:
    MD::Protein* protein; // Pointer to the protein this interactor works on
    MD::ProteinRenderer* renderer; // Renderer associated with the protein worked on
    EnergyCalculator* energyCalculator; // Pointer to energy calculator associated with protein
    UndoBuffer* undoBuffer; // Pointer to the global undo buffer
    
    UpdateDirection updateDirection; // Default update direction when selecting a new structure
    
    MD::Protein::StructureSelector selectedStructure; // The selected secondary structure
    MD::Protein::Residue* selectedResidue; // A selected residue inside the protein
    
    DragBox box; // 3D interaction widget for the selected secondary structure
    StructureList activeCoils; // List of active coil regions in the protein
    
    int numJoints;
    IK* ik; // Pointer to an IK object calculating dihedral angle updates for the protein
    JointMapItem* jointMap; // Array mapping from joint indices to backbone dihedral angles
    int firstResidueIndex;
    int numResidues;
    MD::Scalar* currentPhis;
    MD::Scalar* currentPsis;
    MD::Scalar* deltaPhis;
    MD::Scalar* deltaPsis;
    int ikUpdateDirection;
    int effectorJointIndex;
    Transformation initialTransformation;
    Transformation correctionTransformation;
    IK::Scalar effectorStepSize;
    IK::Scalar correctionStepSize;
    
    /* Constructors and destructors: */
    public:
    ProteinInteractor(MD::Protein* sProtein,MD::ProteinRenderer* sRenderer,EnergyCalculator* sEnergyCalculator,UndoBuffer* sUndoBuffer);
    ~ProteinInteractor(void);
    
    /* Methods: */
    void resetDragBox(void); // Resets the dragbox to the currently selected structure
    void clearCoilRegions(void); // Deactivates all coil regions
    void resetCoilRegions(void); // Resets the active coil regions to the default
    void toggleUpdateDirection(void); // Changes the default update direction
    int getUpdateDirection(void) const // Returns the current update direction for set/changeDihedralAngles
    {
        if(updateDirection==LEFT)
            return -1;
        else
            return 1;
    };
    void selectStructure(const MD::Protein::StructureSelector& newSelectedStructure); // Selects a new secondary structure
    MD::Protein::StructureSelector getStructure(void) const // Returns the currently selected structure
    {
        return selectedStructure;
    };
    bool isValid(void) const // Returns true if a secondary structure is selected
    {
        return selectedStructure.isValid();
    };
    bool isAlphaHelix(void) const // Returns true if the currently selected structure is an alpha helix
    {
        return selectedStructure.isValid()&&selectedStructure.getStructureType()==MD::Protein::SecondaryStructure::ALPHA_HELIX;
    };
    bool isBetaStrand(void) const // Returns true if the currently selected structure is a beta strand
    {
        return selectedStructure.isValid()&&selectedStructure.getStructureType()==MD::Protein::SecondaryStructure::BETA_STRAND;
    };
    bool isCoil(void) const // Returns true if the currently selected structure is a coil region
    {
        return selectedStructure.isValid()&&selectedStructure.getStructureType()==MD::Protein::SecondaryStructure::COIL;
    };
    DragBox& getBox(void) // Returns the drag box
    {
        return box;
    };
    const MD::Point& getBoxRotateCenter(void) const // Returns center of rotation of drag box
    {
        return box.getRotateCenter();
    };
    void setDihedralAngles(const MD::Scalar phis[],const MD::Scalar psis[]); // Sets dihedral angles inside the selected structure
    void changeDihedralAngles(const MD::Scalar deltaPhis[],const MD::Scalar deltaPsis[]); // Changes dihedral angles inside the selected structure
    void toggleCoil(const MD::Protein::StructureSelector& coil); // Toggles the "active" state of a coil region
    const std::list<MD::Protein::StructureSelector>& getLeftCoils(void) const // Returns left active coil regions
    {
        return activeCoils;
    };
    const std::list<MD::Protein::StructureSelector>& getRightCoils(void) const // Returns right active coil regions
    {
        return activeCoils;
    };
    bool startInteraction(void); // Prepares for performing IK on the active coil regions
    bool pickDragBox(const DragBox::HTransformation& modelView,const DragBox::HTransformation& projection,const DragBox::Point& mouseClip);
    bool pickDragBox(const DragBox::Transformation& initialDragTransformation);
    void dragBox(const DragBox::HTransformation& modelView,const DragBox::HTransformation& projection,const DragBox::Point& mouseClip)
    {
        box.drag(modelView,projection,mouseClip);
    };
    void dragBox(const DragBox::Transformation& dragTransformation)
    {
        box.drag(dragTransformation);
    };
    const Transformation& getDragTransformation(void) const
    {
        return box.getDragTransformation();
    };
    void drag(const Transformation& goalTransformation); // Performs a step of IK on the active coil regions
    void applyChanges(void); // Applies IK changes to protein
    void releaseDragBox(void)
    {
        box.release();
    };
    void finishInteraction(void); // Finishes performing IK
    void glRenderAction(GLContextData& contextData) const; // Renders the interaction object
    void glTransformProtein(GLContextData& contextData) const; // Loads a transformation matrix for visually offsetting the protein
    bool readCoilRegionAngles(const char* angleFilename); // Reads dihedral angles from a file
    bool writeCoilRegionAngles(const char* angleFilename) const; // Writes dihedral angles of active coil regions to a file
};

#endif
