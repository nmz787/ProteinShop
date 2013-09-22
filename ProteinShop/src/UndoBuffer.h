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
UndoBuffer - Class to encapsulate undo and redo functionality for
dihedral angle updates on proteins.
***********************************************************************/

#ifndef UNDOBUFFER_INCLUDED
#define UNDOBUFFER_INCLUDED

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
#include <Geometry/OrthonormalTransformation.h>

#include "Protein.h"

class UndoBuffer
{
    /* Embedded classes: */
    private:
    class DihedralAngleSequence // Class to store dihedral angles for a sequence of residues
    {
        /* Elements: */
        private:
        int firstResidue,numResidues; // Index of first residue and number of residues in sequence
        public:
        MD::Scalar* phis; // C-style array of phi angles
        MD::Scalar* psis; // C-style array of psi angles
        
        /* Constructors and destructors: */
        public:
        DihedralAngleSequence(void) // Creates empty angle sequence
            :firstResidue(-1),numResidues(0),
             phis(0),psis(0)
        {
        };
        ~DihedralAngleSequence(void)
        {
            delete[] phis;
            delete[] psis;
        };
        
        /* Methods: */
        int getFirstResidue(void) const // Returns index of first residue in sequence
        {
            return firstResidue;
        };
        int getNumResidues(void) const // Returns number of residues in sequence
        {
            return numResidues;
        };
        void resize(int newFirstResidue,int newNumResidues) // Sets new sequence start and length
        {
            if(newNumResidues!=numResidues)
            {
                /* Delete old arrays: */
                delete[] phis;
                delete[] psis;
                
                /* Allocate new arrays: */
                numResidues=newNumResidues;
                phis=new MD::Scalar[numResidues];
                psis=new MD::Scalar[numResidues];
            }
            
            /* Change start index: */
            firstResidue=newFirstResidue;
        };
    };
    
    class UndoQueueEntry // Entry in the undo queue
    {
        /* Elements: */
        private:
        MD::Protein* protein; // Pointer to manipulated protein
        Geometry::OrthonormalTransformation<double,3> proteinTransformation; // Transformation to apply to entire protein before dihedral angle undo
        std::list<DihedralAngleSequence> sequences[2]; // Lists of sequences of angle deltas (updated to left/right, respectively)
        
        /* Constructors and destructors: */
        public:
        UndoQueueEntry(void) // Creates invalid undo queue entry
            :protein(0)
        {
        };
        
        /* Methods: */
        void set(MD::Protein* sProtein,const Geometry::OrthonormalTransformation<double,3>& sProteinTransformation,const std::list<DihedralAngleSequence> previousSequences[2]); // Calculates angle deltas from given previous values and current values
        MD::Protein* getProtein(void) const
        {
            return protein;
        };
        void undo(void);
        void redo(void);
    };
    
    /* Elements: */
    std::list<UndoQueueEntry> queue; // Undo queue
    std::list<UndoQueueEntry>::iterator queueTail; // Iterator to tail of queue (can be in the middle due to redos)
    MD::Protein* currentProtein; // Pointer to currently manipulated protein
    std::list<DihedralAngleSequence> previousSequences[2]; // Lists of sequences of dihedral angles before update (updated to left/right, respectively)
    int lastSavedDepth; // Distance of last saved state from currently active state
    
    /* Constructors and destructors: */
    public:
    UndoBuffer(void); // Creates empty undo buffer
    ~UndoBuffer(void);
    
    /* Methods: */
    bool isSaved(void) const // Returns true if the current undo state is the one last written to file
    {
        return lastSavedDepth==0;
    };
    bool canUndo(void) const // Returns true if at least one interaction can be undone
    {
        return std::list<UndoQueueEntry>::const_iterator(queueTail)!=queue.begin();
    };
    bool canRedo(void) const // Returns true if at least one undone interaction can be redone
    {
        return std::list<UndoQueueEntry>::const_iterator(queueTail)!=queue.end();
    };
    void saveCurrentState(void) // Marks the current state as "saved to file"
    {
        lastSavedDepth=0;
    };
    void startInteraction(MD::Protein* protein); // Starts interaction on given protein, with empty saved angle list
    void addResidueSequence(int startResidue,int numResidues,int direction); // Saves the given sequence of dihedral angles for subsequent interaction
    void startInteraction(MD::Protein* protein,MD::Protein::Residue* selectedResidue,int direction); // Prepares for interaction on the selected residue
    void startInteraction(MD::Protein* protein,int startResidue,int numResidues,int direction); // Prepares for interaction on the selected residue range
    void startInteraction(const MD::Protein::StructureSelector& selectedStructure,int direction); // Prepares for interaction on the selected structure
    void startInteraction(const std::list<MD::Protein::StructureSelector>& leftStructures,const std::list<MD::Protein::StructureSelector>& rightStructures); // Prepares for interaction on two lists of selected structures (to left and right, respectively)
    void finishInteraction(void); // Finishes last interaction
    void finishInteraction(const Geometry::OrthonormalTransformation<double,3>& proteinTransformation); // Ditto, with overall transformation
    void undo(void); // Undoes the last interaction
    void redo(void); // Redoes the last undone interaction
    void clear(void); // clean up all
    void deleteProtein(const MD::Protein* protein); // Removes all changes in the given protein from the undo queue
};

#endif
