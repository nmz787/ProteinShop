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
Protein - Class to represent a protein as a single chain of amino acid
residues.
***********************************************************************/

#ifndef PROTEIN_INCLUDED
#define PROTEIN_INCLUDED

#include <utility>
#include <vector>
#include <Math/Constants.h>

#include "DragBox.h"
#include "DistanceRange.h"
#include "MDGeometry.h"
#include "Atom.h"

namespace MD {

class Protein
{
    /* Embedded classes: */
    public:
    class ChainAtom;
    class Residue;
    class SecondaryStructure;
    class Dipole;
    class ResidueCreator;
    class ConstAtomIterator;
    class ConstResidueIterator;
    class ConstStructureIterator;
    class BackboneIterator;
    class StructureSelector;
    class BondRotator;
    
    class ChainAtom:public Atom // Class for atoms inside proteins
    {
        friend class Protein;
        friend class Residue;
        friend class Dipole;
        friend class ResidueCreator;
        friend class ConstAtomIterator;
        friend class BackboneIterator;
        friend class StructureSelector;
        friend class BondRotator;
        friend class ProteinRenderer;
        friend class ProteinFactory;
        
        /* Elements: */
        private:
        int atomIndex; // Index of atom in protein chain
        char placement[4]; // Atom placement tag taken from the PDB file
        char name[5]; // Atom name taken from the PDB file
        Residue* residue; // Pointer to amino acid residue this atom belongs to
        int backboneIndex; // Index of atom in residue backbone triple; -1 for non-backbone atoms
        ChainAtom* pred; // Pointer to previous chain atom in list
        ChainAtom* succ; // Pointer to next chain atom in list
        unsigned int rotationMarker; // Marker to stop molecule traversal during rotations
        
        /* Constructors and destructors: */
        public:
        ChainAtom(Element sType,const Position& sPosition,int sAtomIndex,const char* sPlacement,Residue* sResidue,int sBackboneIndex =-1);
        
        /* Methods: */
        int getAtomIndex(void) const
        {
            return atomIndex;
        };
        const char* getPlacement(void) const
        {
            return placement;
        };
        const char* getAtomName(void) const
        {
            return name;
        };
        Residue* getResidue(void) const
        {
            return residue;
        };
        int getBackboneIndex(void) const
        {
            return backboneIndex;
        };
    };
    
    class Residue // Class for amino acid residues
    {
        friend class Protein;
        friend class SecondaryStructure;
        friend class Dipole;
        friend class ResidueCreator;
        friend class ConstResidueIterator;
        friend class BackboneIterator;
        friend class StructureSelector;
        friend class BondRotator;
        friend class ProteinRenderer;
        friend class ProteinFactory;
        
        /* Embedded classes: */
        public:
        enum AminoAcid // Enumerated type for standard amino acids (PDB abbreviations)
        {
            ACE=0,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,
            HIS,ILE,LEU,LYS,MET,NME,PHE,PRO,SER,THR,TRP,TYR,UNK,VAL
        };
        enum ResidueType // Enumerated residue type for standard amino acids 
        {
			ALIPHATIC=0,AROMATIC,CHARGED,NEGATIVE,PHILIC,PHOBIC,POLAR,
			POSITIVE,BASIC,DISULFIDE
        };
        
        /* Static elements: */
        static const char* const commonNames[23]; // Common names of amino acids
        static const char abbreviatedNames[23][4]; // PDB abbreviations of common names
        static const char abbreviatedLcNames[23][4]; // Lower case PDB abbreviations of common names
        static const char singleLetterNames[23]; // PDB short abbreviations of common names
		static const char* resTypePhilic[12];  // Hydropfilic
		static const char* resTypePhobic[12];  // Hydrophobic
        
        /* Elements: */
        private:
        AminoAcid type; // Type of this residue
        ResidueType resType; // Type of this residue
        char residueName[4]; // Residue name taken from the PDB file
        int residueIndex; // Index of this residue
        ChainAtom* beginAtom; // Pointer to first atom in residue
        ChainAtom* endAtom; // Pointer behind last atom in residue
        std::vector<ChainAtom*> backbone; // List of backbone atoms (usually N-C-C)
        SecondaryStructure* secondaryStructure; // Pointer to secondary structure this residue belongs to
        Residue* pred; // Pointer to previous residue in list
        Residue* succ; // Pointer to next residue in list
		char chainID[1];
		int modelID;
        
        /* Constructors and destructors: */
        Residue(const char* sResidueName,int sResidueIndex,SecondaryStructure* sSecondaryStructure,char* chainId,int modelId);
        
        /* Methods: */
        public:
        static AminoAcid parseType(const char* abbreviatedName); // Converts an abbreviated name into an amino acid type
        static ResidueType parseResType(const char* abbreviatedName); // Converts an abbreviated name into an res type
        AminoAcid getType(void) const // Returns the residue's type
        {
            return type;
        };
        ResidueType getResType(void) const // Returns the residue's type
        {
            return resType;
        };
        const char* getPdbResidueName(void) const // Returns the residue's PDB name
        {
            return residueName;
        };
        int getPdbResidueIndex(void) const // Returns the residue's PDB index
        {
            return residueIndex;
        };
        SecondaryStructure* getSecondaryStructure(void) const
        {
            return secondaryStructure;
        };
        Residue* getPred(void) // Returns previous residue
        {
            return pred;
        };
        Residue* getSucc(void) // Returns next residue
        {
            return succ;
        };
		int getModelId(void) const { return modelID; }
		const char* getChainId(void) const { return chainID; }
        int getNumAtoms(void) const; // Returns number of atoms in residue
        Dipole getAmide(void) const; // Returns a dipole object for the residue's amide group
        Dipole getCarboxyl(void) const; // Returns a dipole object for the residue's carboxyl group
    };
    
    class SecondaryStructure // Class to separate proteins into secondary structures
    {
        friend class Protein;
        friend class ResidueCreator;
        friend class ConstStructureIterator;
        friend class StructureSelector;
        friend class ProteinRenderer;

        /* Embedded classes: */
        public:
        enum StructureType
        {
            NONE,COIL,ALPHA_HELIX,BETA_STRAND
        };
        
        /* Elements: */
        private:
        StructureType structureType; // Type of this structure
        Residue* residueBegin; // Pointer to first residue in structure
        Residue* residueEnd; // Pointer to one past last residue in structure
        std::vector<Dipole>::iterator amidesBegin; // Iterator to the first amide in structure
        std::vector<Dipole>::iterator amidesEnd; // Iterator one past the last amide in structure
        std::vector<Dipole>::iterator carboxylsBegin; // Iterator to the first carboxyl in structure
        std::vector<Dipole>::iterator carboxylsEnd; // Iterator to one past the last carboxyl in structure
        SecondaryStructure* pred; // Pointer to previous structure in list
        SecondaryStructure* succ; // Pointer to next structure in list
        
        /* Constructors and destructors: */
        SecondaryStructure(StructureType sStructureType) // Starts building a secondary structure
            :structureType(sStructureType),residueBegin(0),residueEnd(0),pred(0),succ(0)
        {
        };
        
        /* Methods: */
        public:
        StructureType getStructureType(void) const
        {
            return structureType;
        };
        int getFirstResidueIndex(void) const; // Returns index of first residue in structure
        int getNumResidues(void) const; // Returns number of residues in structure
    };
    
    class Dipole // Class for atom pairs potentially part of hydrogen bonds
    {
        friend class Protein;
        friend class ProteinRenderer;
        
        /* Elements: */
        private:
        static const DistanceRange hydrogenBondRangeMajor; // Distance range between major1/minor2 hydrogen bond atoms
        static const DistanceRange hydrogenBondRangeMinor; // Distance range between minor1/minor2 hydrogen bond atoms
        static const Scalar NHdist; // Distance between N and H atoms in amide group?
        static const Scalar cosAngleMin; // Cosine of maximum angle between N-H and H--O bonds
        
        ChainAtom* majorAtom; // Pointer to major atom (N in N-H, C in C=O)
        ChainAtom* minorAtom; // Pointer to minor atom (H in N-H, O in C=O)
        int polarity; // Polarity indicator (+1 for N-H, -1 for C=O)
        
        /* Constructors and destructors: */
        public:
        Dipole(void) // Constructs invalid dipole
            :majorAtom(0),minorAtom(0)
        {
        };
        Dipole(ChainAtom* sMajorAtom,ChainAtom* sMinorAtom,int sPolarity) // Elementwise
            :majorAtom(sMajorAtom),minorAtom(sMinorAtom),polarity(sPolarity)
        {
        };
        
        /* Methods: */
        bool isValid(void) const // Checks if the dipole is correctly formed
        {
            return majorAtom!=0&&minorAtom!=0;
        };
        const ChainAtom* getMajorAtom(void) const
        {
            return majorAtom;
        };
        const ChainAtom* getMinorAtom(void) const
        {
            return minorAtom;
        };
        const ChainAtom* getAlphaCarbon(void) const // Returns pointer to the alpha carbon closest to the major atom
        {
            /*******************************************************
            If the atom is not part of the backbone, we still return
            the correct alpha carbon atom. It just might not be very
            useful for the caller. So there.
            *******************************************************/
            
            return majorAtom->residue->backbone[1];
        };
        Point getBondSite(void) const; // Returns a point halfway between this dipole and a potential bonding partner
        friend bool formHydrogenBond(const Dipole& amide,const Dipole& carboxyl); // Checks if two dipoles can form a hydrogen bond
    };
    
    class ResidueCreator // Class to create residues (and proteins) one atom at a time
    {
        /* Static elements: */
        private:
        static const double covalentBondThreshold; // Maximum distance between two atoms to form a covalent bond
        
        /* Elements: */
        Protein* protein; // Pointer to protein being constructed
        ChainAtom* lastAtom; // Pointer to last atom currently in protein
        ChainAtom* lastBackboneAtom; // Last backbone atom added to the protein
        SecondaryStructure* currentSecondaryStructure; // The secondary structure being constructed
        Residue* currentResidue; // Residue being constructed
        bool firstAtom; // True before first atom has been added to current residue
        bool firstResidue; // True before first residue has been added to current secondary structure
        bool addHydrogens; // true if missing amide protons should be added (David Bernick)
        Scalar fabricatedAtomDistance; // distance of fabricated hydrogen atom from the backbone nitrogen
        int fabricatedAtomIndex; // used when missing amide protons are being inserted
        
        // take any steps necessary to finalize a residue when we are done creating it
        // do nothing if currentResidue is null, otherwise, operate on currentResidue
        void finalizeResidue();
        
        /* Constructors and destructors: */
        public:
        ResidueCreator(Protein* sProtein); // Constructs a creator attached to an empty protein
        ~ResidueCreator(void);
        
        /* Methods: */
        void newSecondaryStructure(SecondaryStructure::StructureType sStructureType); // Starts creating a new secondary structure
        void newResidue(const char* newAbbreviatedName,int newResidueIndex,char*chainId=0, int modelId=0); // Starts adding a new residue
        void addAtom(const char* elementName,int atomIndex,const Position& atomPosition,const char* placement); // Adds a new atom to the protein
        void finishProtein(void); // Finishes creating a protein
		int padHydrogenAtoms(int residueIndex,int numAtoms, const Position&	atomPos, int type);
     };
    
    class ConstAtomIterator // Class to iterate through a protein's atoms
    {
        friend class Protein;
        friend class ConstResidueIterator;
        
        /* Element: */
        private:
        const ChainAtom* atom; // Currently referenced atom
        
        /* Constructors and destructors: */
        public:
        ConstAtomIterator(void) // Creates invalid iterator
            :atom(0)
        {
        };
        private:
        ConstAtomIterator(const ChainAtom* sAtom) // Constructs iterator from pointer to atom
            :atom(sAtom)
        {
        };
        
        /* Methods: */
        public:
        friend bool operator==(const ConstAtomIterator& it1,const ConstAtomIterator& it2) // Equality operator
        {
            return it1.atom==it2.atom;
        };
        friend bool operator!=(const ConstAtomIterator& it1,const ConstAtomIterator& it2) // Inequality operator
        {
            return it1.atom!=it2.atom;
        };
        const ChainAtom& operator*(void) const // Indirection operator
        {
            return *atom;
        };
        const ChainAtom* operator->(void) const // Element access operator
        {
            return atom;
        };
        ConstAtomIterator& operator++(void) // Pre-increment operator
        {
            atom=atom->succ;
            return *this;
        };
        ConstAtomIterator operator++(int) // Post-increment operator
        {
            ConstAtomIterator result=*this;
            atom=atom->succ;
            return result;
        };
        ConstAtomIterator& operator--(void) // Pre-decrement operator
        {
            atom=atom->pred;
            return *this;
        };
        ConstAtomIterator operator--(int) // Post-decrement operator
        {
            ConstAtomIterator result=*this;
            atom=atom->pred;
            return result;
        };
    };
    
    class ConstResidueIterator // Class to iterate through a protein's residues
    {
        friend class Protein;
        friend class ConstStructureIterator;
        
        /* Element: */
        private:
        const Residue* residue; // Currently referenced residue
        
        /* Constructors and destructors: */
        public:
        ConstResidueIterator(void) // Creates invalid iterator
            :residue(0)
        {
        };
        private:
        ConstResidueIterator(const Residue* sResidue) // Constructs iterator from pointer to residue
            :residue(sResidue)
        {
        };
        
        /* Methods: */
        public:
        friend bool operator==(const ConstResidueIterator& it1,const ConstResidueIterator& it2) // Equality operator
        {
            return it1.residue==it2.residue;
        };
        friend bool operator!=(const ConstResidueIterator& it1,const ConstResidueIterator& it2) // Inequality operator
        {
            return it1.residue!=it2.residue;
        };
        const Residue& operator*(void) const // Indirection operator
        {
            return *residue;
        };
        const Residue* operator->(void) const // Element access operator
        {
            return residue;
        };
        ConstResidueIterator& operator++(void) // Pre-increment operator
        {
            residue=residue->succ;
            return *this;
        };
        ConstResidueIterator operator++(int) // Post-increment operator
        {
            ConstResidueIterator result=*this;
            residue=residue->succ;
            return result;
        };
        ConstResidueIterator& operator--(void) // Pre-decrement operator
        {
            residue=residue->pred;
            return *this;
        };
        ConstResidueIterator operator--(int) // Post-decrement operator
        {
            ConstResidueIterator result=*this;
            residue=residue->pred;
            return result;
        };
        int getNumAtoms(void) const // Returns number of atoms in residue
        {
            return residue->getNumAtoms();
        };
        ConstAtomIterator atomsBegin(void) const // Returns iterator to first atom in residue
        {
            return ConstAtomIterator(residue->beginAtom);
        };
        ConstAtomIterator atomsEnd(void) const // Returns iterator behind last atom in residue
        {
            return ConstAtomIterator(residue->endAtom);
        };
    };
    
    class ConstStructureIterator // Class to iterate through a protein's secondary structures
    {
        friend class Protein;
        
        /* Element: */
        private:
        const SecondaryStructure* structure; // Currently referenced secondary structure
        
        /* Constructors and destructors: */
        public:
        ConstStructureIterator(void) // Creates invalid iterator
            :structure(0)
        {
        };
        private:
        ConstStructureIterator(const SecondaryStructure* sStructure) // Constructs iterator from pointer to secondary structure
            :structure(sStructure)
        {
        };
        
        /* Methods: */
        public:
        friend bool operator==(const ConstStructureIterator& it1,const ConstStructureIterator& it2) // Equality operator
        {
            return it1.structure==it2.structure;
        };
        friend bool operator!=(const ConstStructureIterator& it1,const ConstStructureIterator& it2) // Inequality operator
        {
            return it1.structure!=it2.structure;
        };
        const SecondaryStructure& operator*(void) const // Indirection operator
        {
            return *structure;
        };
        const SecondaryStructure* operator->(void) const // Element access operator
        {
            return structure;
        };
        ConstStructureIterator& operator++(void) // Pre-increment operator
        {
            structure=structure->succ;
            return *this;
        };
        ConstStructureIterator operator++(int) // Post-increment operator
        {
            ConstStructureIterator result=*this;
            structure=structure->succ;
            return result;
        };
        ConstStructureIterator& operator--(void) // Pre-decrement operator
        {
            structure=structure->pred;
            return *this;
        };
        ConstStructureIterator operator--(int) // Post-decrement operator
        {
            ConstStructureIterator result=*this;
            structure=structure->pred;
            return result;
        };
        int getNumResidues(void) const // Returns number of residues in secondary structure
        {
            return structure->getNumResidues();
        };
        ConstResidueIterator residuesBegin(void) const // Returns iterator to first residue in secondary structure
        {
            return ConstResidueIterator(structure->residueBegin);
        };
        ConstResidueIterator residuesEnd(void) const // Returns iterator behind last residue in secondary structure
        {
            return ConstResidueIterator(structure->residueEnd);
        };
    };
    
    class BackboneIterator // Class to iterate through a protein's backbone bonds
    {
        friend class BondRotator;
        
        /* Embedded classes: */
        public:
        enum BondType // Enumerated type for types of backbone bonds
        {
            UNKNOWN,PEPTIDE,PHI,PSI
        };
        
        /* Elements: */
        private:
        std::vector<ChainAtom*>::iterator atom1,atom2; // The two (backbone) atoms participating in the bond
        
        /* Constructors and destructors: */
        public:
        BackboneIterator(void) // Creates invalid backbone iterator
        {
        };
        BackboneIterator(const std::vector<ChainAtom*>::iterator& sAtom1,const std::vector<ChainAtom*>::iterator& sAtom2) // Creates backbone iterator pointing to the bond between the two given (backbone) atoms
            :atom1(sAtom1),atom2(sAtom2)
        {
        };
        
        /* Methods: */
        bool isPeptideBond(void) const // Returns true if the bond straddles two residues
        {
            return (*atom1)->residue!=(*atom2)->residue;
        };
        BondType getBondType(void) const // Returns type of this bond
        {
            if((*atom1)->residue!=(*atom2)->residue)
                return PEPTIDE;
            else if((*atom1)->residue->type==Residue::UNK)
                return UNKNOWN;
            else if((*atom1)->backboneIndex==0&&(*atom2)->backboneIndex==1)
                return PHI;
            else if((*atom1)->backboneIndex==1&&(*atom2)->backboneIndex==2)
                return PSI;
            else
                return UNKNOWN;
        };
        bool isRotatable(void) const // Returns true if the bond can be freely rotated
        {
            if(isPeptideBond()) // Peptide bonds are always rigid
                return false;
            
            if((*atom1)->residue->type==Residue::UNK) // Flag all backbone bonds of unknown residues as rigid
                return false;
            
            /* Check if the bond is part of a ring structure (as in Proline, for example): */
            if((*atom1)->residue->type==Residue::PRO&&(*atom1)->backboneIndex==0&&(*atom2)->backboneIndex==1)
                return false;
            
            return true;
        };
        friend bool operator==(const BackboneIterator& it1,const BackboneIterator& it2)
        {
            return it1.atom1==it2.atom1&&it1.atom2==it2.atom2;
        };
        friend bool operator!=(const BackboneIterator& it1,const BackboneIterator& it2)
        {
            return it1.atom1!=it2.atom1&&it1.atom2!=it2.atom2;
        };
        ChainAtom* getAtom1(void) const // Returns first atom participating in bond
        {
            return (*atom1);
        };
        ChainAtom* getAtom2(void) const // Returns second atom participating in bond
        {
            return (*atom2);
        };
        int getResidueIndex(void) const; // Returns index of residue containing first atom
        Residue* getResidue(void) const // Returns pointer to residue containing first atom
        {
            return (*atom1)->residue;
        };
        BackboneIterator& operator++(void) // Pre-increment operator
        {
            atom1=atom2;
            Residue* rPtr=(*atom2)->residue;
            ++atom2;
            if(atom2==rPtr->backbone.end())
            {
                rPtr=rPtr->succ;
                if(rPtr!=0)
                    atom2=rPtr->backbone.begin();
            }
            
            return *this;
        };
        BackboneIterator operator++(int) // Post-increment operator
        {
            BackboneIterator result=*this;
            operator++();
            return result;
        };
        BackboneIterator& operator--(void) // Pre-decrement operator
        {
            atom2=atom1;
            Residue* rPtr=(*atom1)->residue;
            if(atom1==rPtr->backbone.begin())
            {
                rPtr=rPtr->pred;
                if(rPtr!=0)
                    atom1=rPtr->backbone.end();
            }
            --atom1;
            
            return *this;
        };
        BackboneIterator operator--(int) // Post-decrement operator
        {
            BackboneIterator result=*this;
            operator--();
            return result;
        };
    };
    
    class StructureSelector // Class to encapsulate selection/modification of secondary structures
    {
        friend class Protein;
        friend class ProteinRenderer;

        /* Elements: */
        private:
        Protein* protein; // Pointer to the protein
        SecondaryStructure* structure; // Pointer to the selected structure
        int structureIndex; // Index of the selected structure inside the protein
        int firstResidueIndex,numResidues; // Index of first residue and number of residues in structure
        BackboneIterator begin,end; // Backbone iterators for the selected structure's residues
        
        /* Constructors and destructors: */
        public:
        StructureSelector(void) // Creates an invalid structure selector
            :protein(0),structure(0),structureIndex(-1),firstResidueIndex(-1),numResidues(0)
        {
        };
        private:
        StructureSelector(Protein* sProtein,SecondaryStructure* sStructure,int sStructureIndex);
        
        /* Methods: */
        public:
        bool isValid(void) const // Checks if a selector is valid
        {
            return protein!=0&&structure!=0;
        };
        Protein* getProtein(void) const // Returns the protein containing the referenced structure
        {
            return protein;
        };
        bool containsResidue(const Residue* rPtr) const // Tests if the given residue is contained in the selected structure
        {
            /* Validity checks: */
            if(!isValid()||rPtr==0)
                return false;

            return rPtr->secondaryStructure==structure;
        };
        SecondaryStructure::StructureType getStructureType(void) const // Returns type of the indicated secondary structure
        {
            return structure->getStructureType();
        };
        StructureSelector getPred(void) const // Returns previous structure
        {
            return StructureSelector(protein,structure->pred,structureIndex-1);
        };
        StructureSelector getSucc(void) const // Returns next structure
        {
            return StructureSelector(protein,structure->succ,structureIndex+1);
        };
        int getStructureIndex(void) const // Returns index of this structure
        {
            return structureIndex;
        };
        int getStructureTypeIndex(void) const; // Returns index of this structure amongst its own type
        int getFirstResidueIndex(void) const // Returns index of the first residue in this structure
        {
            return firstResidueIndex;
        };
        int getNumResidues(void) const // Returns the number of residues in this structure
        {
            return numResidues;
        };
        int getLastResidueIndex(void) const // Returns index of the last residue in this structure
        {
            return firstResidueIndex + numResidues -1;
        };
        BackboneIterator beginBackbone(void) const // Returns backbone iterator for the beginning of this structure
        {
            return begin;
        };
        BackboneIterator endBackbone(void) const // Returns backbone iterator for the end of this structure
        {
            return end;
        };
        void getDihedralAngles(Scalar phis[],Scalar psis[]) const // Returns dihedral angles of all residues in structure
        {
            protein->getDihedralAngles(firstResidueIndex,numResidues,phis,psis);
        };
        void changeDihedralAngles(const Scalar deltaPhis[],const Scalar deltaPsis[],int direction) const // Changes dihedral angles of all residues in structure
        {
            protein->changeDihedralAngles(firstResidueIndex,numResidues,deltaPhis,deltaPsis,direction);
        };
        void setDihedralAngles(const Scalar phis[],const Scalar psis[],int direction) const // Sets dihedral angles of all residues in structure
        {
            protein->setDihedralAngles(firstResidueIndex,numResidues,phis,psis,direction);
        };
        DragBox createDragBox(void) const; // Returns a drag box enclosing the selected structure
    };
    
    class BondRotator // Class to encapsulate rotating of backbone bonds for protein manipulation
    {
        friend class Protein;
        
        /* Embedded classes: */
        public:
        enum IntersectionType // Enumerated type for different intersections
        {
            NONE,NON_INTERACTIVE,TOP,BOTTOM,SIDE,POINT_INSIDE
        };
        
        /* Elements: */
        private:
        Protein* protein; // Pointer to the protein
        std::vector<ChainAtom*>::iterator atom1,atom2; // The two (backbone) atoms participating in the bond
        IntersectionType intersectionType; // The intersection type (NONE if bond is not valid)
        Point lastIntersection; // The last intersection point; used to calculate rotations
        Vector axis; // Direction of rotation axis
        
        /* Constructors and destructors: */
        public:
        BondRotator(void) // Constructs invalid rotator
            :intersectionType(NONE)
        {
        };
        BondRotator(Protein* sProtein,std::vector<ChainAtom*>::iterator sAtom1,std::vector<ChainAtom*>::iterator sAtom2) // Creates a non-interactive rotator for the bond between the two specified atoms
            :protein(sProtein),atom1(sAtom1),atom2(sAtom2),intersectionType(NON_INTERACTIVE)
        {
            if(isValid())
            {
                /* Calculate the rotation axis: */
                axis=Vector((*atom2)->getPosition()-(*atom1)->getPosition());
                axis.normalize();
            }
        };
        BondRotator(Protein* sProtein,std::vector<ChainAtom*>::iterator sAtom1,std::vector<ChainAtom*>::iterator sAtom2,IntersectionType sIntersectionType,const Point& sLastIntersection) // Creates a rotator for the bond between the two specified atoms
            :protein(sProtein),atom1(sAtom1),atom2(sAtom2),intersectionType(sIntersectionType),
             lastIntersection(sLastIntersection)
        {
            if(isValid())
            {
                /* Calculate the rotation axis: */
                axis=Vector((*atom2)->getPosition()-(*atom1)->getPosition());
                axis.normalize();
            }
        };
        BondRotator(Protein* sProtein,const BackboneIterator& backboneIt) // Constructs a bond rotator from a backbone iterator
            :protein(sProtein),atom1(backboneIt.atom1),atom2(backboneIt.atom2),intersectionType(NON_INTERACTIVE)
        {
            if(isValid())
            {
                /* Calculate the rotation axis: */
                axis=Vector((*atom2)->getPosition()-(*atom1)->getPosition());
                axis.normalize();
            }
        };
        
        /* Methods: */
        bool isValid(void) const // Checks if the rotator is valid
        {
            return protein!=0&&intersectionType!=NONE;
        };
        Point getStart(void) const // Returns the bond's start point
        {
            return Point((*atom1)->getPosition());
        };
        Point getEnd(void) const // Returns the bond's end point
        {
            return Point((*atom2)->getPosition());
        };
        Vector getAxis(void) const // Returns the bond's rotation axis
        {
            return axis;
        };
        void rotate(Scalar angle); // Rotates the protein around the given bond's axis
        Scalar getDihedralAngle(void) const; // Returns the dihedral angle of a backbone bond
        void setDihedralAngle(Scalar angle) // Sets the dihedral angle of a backbone bond to a specific value
        {
            /* Rotate the bond by the difference between the current and the new angle: */
            Scalar angleDist=angle-getDihedralAngle();
            if(angleDist<-Math::Constants<Scalar>::pi)
                angleDist+=Scalar(2)*Math::Constants<Scalar>::pi;
            else if(angleDist>Math::Constants<Scalar>::pi)
                angleDist-=Scalar(2)*Math::Constants<Scalar>::pi;
            rotate(angleDist);
        };
    };
    
    private:
    struct MoleculeTraversal // Structure to traverse the molecule during rotations
    {
        /* Elements: */
        public:
        ChainAtom* atom; // The atom that is to be visited next
        int rotationIndex; // Index into an array of rotation matrices
        
        /* Constructors and destructors: */
        MoleculeTraversal(ChainAtom* sAtom,int sRotationIndex)
            :atom(sAtom),rotationIndex(sRotationIndex)
        {
        };
    };
    
    friend class ResidueCreator;
    friend class StructureSelector;
    friend class BondRotator;
    friend class ProteinRenderer;
    
    /* Elements: */
    private:
    static Scalar bondRadius; // Radius of cylinders representing bonds and other related stuff
    int numAtoms; // Number of atoms in protein
    ChainAtom* atoms; // List of all atoms in protein
    Residue* residues; // List of all residues in protein in order along the chain
    SecondaryStructure* secondaryStructures; // List of all secondary structures in protein in order along the chain
    BackboneIterator beginIterator; // Iterator to the first backbone bond
    BackboneIterator endIterator; // Iterator one past the last backbone bond
    std::vector<Dipole> amides; // Array holding backbone amide groups
    std::vector<Dipole> carboxyls; // Array holding backbone carboxyl groups
    unsigned int rotationMarker; // Marker to stop molecule traversal during rotations
    ChainAtom** atomPointers; // A standard C array holding pointers to all atoms, in the order they were read from file
    
    /* Constructors and destructors: */
    public:

    typedef std::vector<Dipole>::const_iterator ConstDipoleIterator;
    
    Protein(void);
    ~Protein(void);
    
    /* Methods: */
    static void setBondRadius(Scalar newBondRadius) // Changes the bond radius
    {
        bondRadius=newBondRadius;
    };
    int getNumAtoms(void) const // Returns the number of atoms in the protein
    {
        return numAtoms;
    };
    const ChainAtom* const* getAtomPointers(void) const // Returns an array of pointers to atoms
    {
        return atomPointers;
    };
    ChainAtom** getAtomPointers(void) // Ditto
    {
        return atomPointers;
    };
    Point calcCentroid(void) const; // Calculates the centroid of the protein
    double calcRadius(void) const; // Calculates the radius of the protein around its centroid
    ConstDipoleIterator amidesBegin() const // Returns an iterator to the first amide group
    {
        return amides.begin();
    };
    ConstDipoleIterator amidesEnd() const // Returns an iterator to one past the last amide group
    {
        return amides.end();
    };
    ConstDipoleIterator carboxylsBegin() const // Returns an iterator to the first carboxyl group
    {
        return carboxyls.begin();
    };
    ConstDipoleIterator carboxylsEnd() const // Returns an iterator to one past the last carboxyl group
    {
        return carboxyls.end();
    };
    ConstAtomIterator atomsBegin(void) const // Returns an iterator to the first atom
    {
        return ConstAtomIterator(atoms);
    };
    ConstAtomIterator atomsEnd(void) const // Returns an iterator one past the last atom
    {
        return ConstAtomIterator(0);
    };
    int getNumResidues(void) const; // Returns the number of residues in the protein
    std::pair<int,int> getResidueIndexRange(void) const; // Returns range of PDB residue indices
    ConstResidueIterator residuesBegin(void) const // Returns an iterator to the first residue
    {
        return ConstResidueIterator(residues);
    };
    ConstResidueIterator residuesEnd(void) const // Returns an iterator one past the last residue
    {
        return ConstResidueIterator(0);
    };
    int getNumStructures(void) const; // Returns the number of secondary structures in the protein
    ConstStructureIterator structuresBegin(void) const // Returns an iterator to the first secondary structure
    {
        return ConstStructureIterator(secondaryStructures);
    };
    ConstStructureIterator structuresEnd(void) const // Returns an iterator one past the last secondary structure
    {
        return ConstStructureIterator(0);
    };
    BackboneIterator backboneBegin(void) // Returns an iterator to the first backbone bond
    {
        return beginIterator;
    };
    BackboneIterator backboneEnd(void) // Returns an iterator one past the last backbone bond
    {
        return endIterator;
    };
    ConstAtomIterator pickAtom(const Ray& ray) const; // Selects an atom by shooting a ray
    const Residue* pickResidue(int pdbResidueIndex) const; // Returns pointer to a residue by PDB index
    Residue* pickResidue(int pdbResidueIndex); // Ditto
    Residue* pickResidue(const Ray& ray); // Selects a residue by shooting a ray
    int getResidueIndex(const Residue* rPtr) const; // Returns protein index of a residue
    void changeResidueStructureType(Residue* rPtr,SecondaryStructure::StructureType newType); // Changes the structure type of a residue, therefore changing the secondary structure sequence
    StructureSelector pickStructure(int structureIndex); // Selects a secondary structure by number
    StructureSelector pickStructure(SecondaryStructure::StructureType structureType,int structureIndex); // Selects a secondary structure of given type by number
    StructureSelector pickStructure(const Residue* residue); // Selects the secondary structure containing a given residue
    StructureSelector pickStructure(const Point& p); // Selects a secondary structure by specifying a point
    StructureSelector pickStructure(const Ray& ray); // Selects a secondary structure by shooting a ray
    BondRotator pickBond(const Point& p); // Selects a backbone bond by specifying a point
    BondRotator pickBond(const Ray& ray); // Selects a backbone bond by specifying a ray
    void getDihedralAngles(int startResidue,int numResidues,Scalar phis[],Scalar psis[]) const; // Retrieves dihedral angles of multiple residues
    void changeDihedralAngles(int startResidue,int numResidues,const Scalar deltaPhis[],const Scalar deltaPsis[],int direction,bool stopAtLastResidue =false); // Changes dihedral angles of multiple residues
    void setDihedralAngles(int startResidue,int numResidues,const Scalar phis[],const Scalar psis[],int direction); // Sets dihedral angles of multiple residues
    template <class TransformationParam>
    void transform(const TransformationParam& t) // Transforms complete protein
    {
        for(ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
            aPtr->setPosition(t.transform(aPtr->getPosition()));
    };
    // directly set the locations of the atoms in the protein
    void setAtomPositions (const Scalar *x, const Scalar *y, const Scalar *z);
};

// Nelson's Stuff
	bool initBuild(const MD::Protein* protein);
	void BuildBeta();
    
}   // end namespace MD

#endif
