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

#define PROTEIN_IMPLEMENTATION

#include <assert.h>
#include <deque>
#include <string.h>
#include <ctype.h>
#include <Math/Math.h>
#include <Math/Constants.h>
#include <Geometry/AffineCombiner.h>
#include <Geometry/Matrix.h>
#include <Geometry/AffineTransformation.h>
#include <Geometry/Cylinder.h>

#include "Globals.h"
#include "Protein.h"
//#define _VERBOSE_
namespace MD {

/***********************************
Methods of class Protein::ChainAtom:
***********************************/

Protein::ChainAtom::ChainAtom(Atom::Element sType,const Position& sPosition,int sAtomIndex,const char* sPlacement,Protein::Residue* sResidue,int sBackboneIndex)
    :Atom(sType,sPosition),atomIndex(sAtomIndex),residue(sResidue),backboneIndex(sBackboneIndex),
     pred(0),succ(0),rotationMarker(0)
{
    strncpy(placement,sPlacement,4);
    placement[3]='\0';
    name[0]='\0';
    strcat(name,getElementName(sType));
    strcat(name,sPlacement);
    name[4]='\0';
}

/*****************************************
Static elements of class Protein::Residue:
*****************************************/

const char* const Protein::Residue::commonNames[23]={
    "ACE","Alanine","Arginine","Asparagine","Aspartic acid","Cysteine",
    "Glutamine","Glutamic acid","Glycine","Histidine","Isoleucine",
    "Leucine","Lysine","Methionine","NME","Phenylalanine","Proline","Serine","Threonine",
    "Tryptophan","Tyrosine","Unknown","Valine"};

const char Protein::Residue::abbreviatedNames[23][4]={
    "ACE","ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","NME","PHE","PRO","SER","THR","TRP","TYR","UNK","VAL"};

const char Protein::Residue::abbreviatedLcNames[23][4]={
    "Ace","Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly",
    "His","Ile","Leu","Lys","Met","Nme","Phe","Pro","Ser","Thr","Trp","Tyr","Unk","Val"};

const char Protein::Residue::singleLetterNames[23]={
    'B','A','R','N','D','C','Q','E','G',
    'H','I','L','K','M','J','F','P','S','T','W','Y','U','V'};

const char* Protein::Residue::resTypePhilic[12]={"ARG","ASN","ASP","GLN",
"GLU","HIS","HID","HIE","HIP","LYS","SER","THR"};
const char* Protein::Residue::resTypePhobic[12]={"ALA","ILE","LEU","MET",
"PHE","PRO","TRP","VAL","TYR","GLY","CYC","CYX"};

/*********************************
Methods of class Protein::Residue:
*********************************/

Protein::Residue::Residue(const char* sResidueName,int sResidueIndex,Protein::SecondaryStructure* sSecondaryStructure, char* chainId, int modelId)
    :residueIndex(sResidueIndex),beginAtom(0),endAtom(0),secondaryStructure(sSecondaryStructure),
     pred(0),succ(0),modelID(modelId)
{
    /* Keep the residue name around (it might be one we don't know): */
    strncpy(residueName,sResidueName,4);
    residueName[3]='\0';
    
	if(chainId)
	{
		strncpy(chainID,chainId,1);
		chainID[1] = '\0';
	}
	else
		chainID[0] = '\0';

    /* Try to parse the residue name: */
    type=parseType(residueName);
	resType = parseResType(residueName);
	sprintf(msg.buf, "%s resType %d \n", residueName, ResidueType(resType));
	msg.Debug(RESIDUEID, msg.buf);
}

Protein::Residue::AminoAcid Protein::Residue::parseType(const char* abbreviatedName)
{
    /* Search through the list of abbreviated names to find the given one: */
    int index;
    for(index=0;index<25;++index)
        if(strcmp(abbreviatedNames[index],abbreviatedName)==0)
            break;
	msg.Debug(RESIDUEID, abbreviatedName, " -- ", commonNames[index]);

    /* If nothing matched, punt: */
    if(index>=25)
        return UNK;
    else return AminoAcid(index);
}

int Protein::Residue::getNumAtoms(void) const
{
    int result=0;
    for(const ChainAtom* aPtr=beginAtom;aPtr!=endAtom;aPtr=aPtr->succ)
        ++result;
    return result;
}

Protein::Residue::ResidueType Protein::Residue::parseResType(const char* abbreviatedName)
{
    /* Search through the list of abbreviated names to match the residue type: */
    int index;
    for(index=0;index<12;++index)
	{
        if(strcmp(resTypePhilic[index],abbreviatedName)==0)
		return ResidueType(4);
        
		if(strcmp(resTypePhobic[index],abbreviatedName)==0)
		return ResidueType(5);
		
		if(strcmp(abbreviatedNames[5],abbreviatedName)==0)
		return ResidueType(9);
	}
	
	return ResidueType(8);
}

Protein::Dipole Protein::Residue::getAmide(void) const
{
    /* Assume that the residue type is neither unknown nor proline: */
    if(type!=UNK&&type!=PRO)
    {
        ChainAtom* major=backbone[0];
        ChainAtom* minor=0;
     
        /* Find hydrogen atom: */
        for(std::vector<Atom*>::const_iterator bIt=major->getBonds().begin();bIt!=major->getBonds().end();++bIt)
            if((*bIt)->getType()==Atom::H)
            {
                minor=static_cast<ChainAtom*>(*bIt);
                break;
            }

        return Dipole(major,minor,1);
    }
    else
        return Dipole(); // Return invalid dipole
}

Protein::Dipole Protein::Residue::getCarboxyl(void) const
{
    /* Assume that the residue type is not unknown: */
    if(type!=UNK)
    {
        ChainAtom* major=backbone[2];
        ChainAtom* minor=0;

        /* Find oxygen atom: */
        for(std::vector<Atom*>::const_iterator bIt=major->getBonds().begin();bIt!=major->getBonds().end();++bIt)
            if((*bIt)->getType()==Atom::O)
            {
                minor=static_cast<ChainAtom*>(*bIt);
                break;
            }

        return Dipole(major,minor,-1);
    }
    else
        return Dipole(); // Return invalid dipole
}

/********************************************
Methods of class Protein::SecondaryStructure:
********************************************/

int Protein::SecondaryStructure::getFirstResidueIndex(void) const
{
    int result=-1;
    for(const Residue* rPtr=residueBegin;rPtr!=0;rPtr=rPtr->pred)
        ++result;
    return result;
}

int Protein::SecondaryStructure::getNumResidues(void) const
{
    int result=0;
    for(const Residue* rPtr=residueBegin;rPtr!=residueEnd;rPtr=rPtr->succ)
        ++result;
    return result;
}

/****************************************
Static elements of class Protein::Dipole:
****************************************/

const DistanceRange Protein::Dipole::hydrogenBondRangeMajor(2.7,3.2);
const DistanceRange Protein::Dipole::hydrogenBondRangeMinor(1.7,2.35);
const Scalar Protein::Dipole::NHdist=Scalar(1.006);
const Scalar Protein::Dipole::cosAngleMin=Scalar(0.707);

/********************************
Methods of class Protein::Dipole:
********************************/

Point Protein::Dipole::getBondSite(void) const
{
    /* Calculate a direction vector from the major to the minor atom: */
    Point minorPos=minorAtom->getPosition();
    Vector dir=minorPos-majorAtom->getPosition();

    /* Offset the minor atom position by half the hydrogen bond distance: */
    double averageBondDistance=(hydrogenBondRangeMinor.getMin()+hydrogenBondRangeMinor.getMax())*0.5;
    return minorPos+dir*Scalar(averageBondDistance*0.5/Geometry::mag(dir));
}

bool formHydrogenBond(const Protein::Dipole& amide,const Protein::Dipole& carboxyl)
{
    /**************************************************************
    This function silently assumes that the first given dipole is
    an amide group and the second given dipole is a carboxyl group.
    **************************************************************/

    /* Check for distance requirements: */
    if(Protein::Dipole::hydrogenBondRangeMajor.areInRange(amide.majorAtom->getPosition(),carboxyl.minorAtom->getPosition())&&
         Protein::Dipole::hydrogenBondRangeMinor.areInRange(amide.minorAtom->getPosition(),carboxyl.minorAtom->getPosition()))
    {
        /* Check for alignment requirements: */
        Vector A=amide.minorAtom->getPosition()-amide.majorAtom->getPosition();
        Vector B=carboxyl.minorAtom->getPosition()-amide.minorAtom->getPosition();
        Scalar cosAngle=(A*B)/Math::sqrt(Geometry::sqr(A)*Geometry::sqr(B));
        return cosAngle>=Protein::Dipole::cosAngleMin;
    }
    else
        return false;
}

/************************************************
Static elements of class Protein::ResidueCreator:
************************************************/

const double Protein::ResidueCreator::covalentBondThreshold=1.85;

/****************************************
Methods of class Protein::ResidueCreator:
****************************************/

Protein::ResidueCreator::ResidueCreator(Protein* sProtein)
    :protein(sProtein),lastAtom(0),lastBackboneAtom(0),currentSecondaryStructure(0),currentResidue(0),
     firstAtom(false),firstResidue(false),
     addHydrogens( retrieveValue(*configFile, "/ProteinCreator/addHydrogens", false) ),
     fabricatedAtomDistance( retrieveValue(*configFile, "/ProteinCreator/fabricatedAtomDistance", 1.1) ),
     fabricatedAtomIndex(1)
{
}

Protein::ResidueCreator::~ResidueCreator(void)
{
}

void Protein::ResidueCreator::finalizeResidue()
{
    // do nothing on the first call
    if ( !currentResidue ) return;

    // check second and subsequent residues for a missing hydrogen atom (amide proton)
    if ( addHydrogens &&
         currentResidue->type != Residue::PRO &&
         currentResidue->backbone.size() > 1 )
    {
        ChainAtom *tN = currentResidue->backbone[0];
        ChainAtom *tCA = currentResidue->backbone[1];
        if ( tN->getType() == Atom::N && tCA->getType() == Atom::C )
        {
            // see if the nitrogen is bonded to a hydrogen
            for ( uint i = 0; i < tN->getBonds().size(); ++i )
                if ( tN->getBonds()[i]->getType() == Atom::H )
                    return;
            Position fabricatedAtomPosition;
            if ( currentResidue->pred && currentResidue->pred->backbone.size() > 0 )
            {
                // fudge the position of the fabricated hydrogen atom
                ChainAtom *pred = currentResidue->pred->backbone.back();
                Point midpoint = mid (pred->getPosition(), tCA->getPosition());
                Vector radial = tN->getPosition() - midpoint;
                fabricatedAtomPosition = midpoint + (
                    radial * (1.0 + fabricatedAtomDistance / radial.mag())
                );
            }
            else
            {
                // special case handling for the first residue
                Vector radial = tN->getPosition() - tCA->getPosition();
                radial.normalize();
                radial *= fabricatedAtomDistance;
                fabricatedAtomPosition = tN->getPosition() + radial;
            }

            // fabricate a new atom for the missing hydrogen
            ++protein->numAtoms;
            ChainAtom *fabH = new ChainAtom (Atom::H, fabricatedAtomPosition,
                                             fabricatedAtomIndex++, "N",
                                             currentResidue);
            fabH->pred = lastAtom;
            if ( lastAtom ) lastAtom->succ = fabH;
            else protein->atoms = fabH;
            lastAtom = fabH;

            // bond the hydrogen to its nitrogen and fabricate a new amide group
            protein->amides.push_back (Dipole(tN, fabH, 1));
            bond (*fabH, *tN);
        }
    }
}

void Protein::ResidueCreator::newSecondaryStructure(Protein::SecondaryStructure::StructureType sStructureType)
{
    /* Create a new secondary structure: */
    SecondaryStructure* newSecondaryStructure=new SecondaryStructure(sStructureType);
    newSecondaryStructure->pred=currentSecondaryStructure;
    if(currentSecondaryStructure!=0)
        currentSecondaryStructure->succ=newSecondaryStructure;
    else
        protein->secondaryStructures=newSecondaryStructure;
    currentSecondaryStructure=newSecondaryStructure;
    
    /* Reset addition flags: */
    firstResidue=true;
}

void Protein::ResidueCreator::newResidue(const char* newAbbreviatedName,int newResidueIndex, char* chainId, int modelId)
{
    // finalize the existing residue, if any
    finalizeResidue();
    #ifdef _VERBOS_
	printf("New residue %s %d\n", newAbbreviatedName,  newResidueIndex);
    #endif
    
    /* Add a new residue to the chain of residues: */
    Residue* newResidue=new	Residue(newAbbreviatedName,newResidueIndex,currentSecondaryStructure, chainId, modelId);
    newResidue->pred=currentResidue;
    if(currentResidue!=0)
        currentResidue->succ=newResidue;
    else
        protein->residues=newResidue;
    currentResidue=newResidue;
    
    /* Reset addition flags: */
    firstAtom=true;
    
    /* Start adding residues to a new secondary structure if necessary: */
    if(firstResidue)
    {
        if(currentSecondaryStructure->pred!=0)
            currentSecondaryStructure->pred->residueEnd=currentResidue;
        currentSecondaryStructure->residueBegin=currentResidue;
        firstResidue=false;
    }
}

void Protein::ResidueCreator::addAtom(const char* elementName,int atomIndex,const Position& atomPosition,const char* placement)
{
    // fabricate atom indices when inserting hydrogen atoms for missing amide protons
    if ( addHydrogens )
        atomIndex = fabricatedAtomIndex++;
    
	// swap the atom name to our system (1HB => HB1)
	if (Atom::parseType(elementName)== 118)
	{
		char* swap = new char[strlen(placement)+1];
		sprintf(swap,"%s%s",placement, elementName);
		strncpy((char*)elementName, swap, 1); 
		placement = ++swap;
	}
    /* Add the new atom to the atom list: */
    ++protein->numAtoms;
	ChainAtom* newAtom=new ChainAtom(Atom::parseType(elementName),atomPosition,atomIndex,placement,currentResidue);
    newAtom->pred=lastAtom;
    if(lastAtom!=0)
        lastAtom->succ=newAtom;
    else
        protein->atoms=newAtom;
    lastAtom=newAtom;
    
    /* Add the new atom to the current residue: */
    if(firstAtom)
    {
        if(currentResidue->pred!=0)
            currentResidue->pred->endAtom=newAtom;
        currentResidue->beginAtom=newAtom;
        firstAtom=false;
    }
    
    /* Add the new atom to the backbone if appropriate: */
    if(newAtom->getType()==Atom::N||newAtom->getType()==Atom::C)
    {
		if((newAtom->getType()==Atom::C&&(placement[0]=='\0'||placement[0]=='A'||placement[0]=='Y'))
			||(newAtom->getType()==Atom::N&&(placement[0]=='\0'|| placement[0]=='T')))
		{
			newAtom->backboneIndex=currentResidue->backbone.size();
            if(lastBackboneAtom!=0)
                bond(*newAtom,*lastBackboneAtom);
            currentResidue->backbone.push_back(newAtom);
            lastBackboneAtom=newAtom;
	 		#ifdef _VERBOSE_
			printf("%s%s  add %d th Atom backbone size %d\n", elementName, placement, atomIndex, currentResidue->backbone.size());
			#endif
		}
    }
    
    /* Bond the new atom with all appropriate existing ones: */
    for(ChainAtom* bondSearch=newAtom;bondSearch!=currentResidue->beginAtom;bondSearch=bondSearch->pred)
    {
        ChainAtom* testAtom=bondSearch->pred;
        
        /* Skip all existing hydrogen atoms - they must already be bonded: */
        if(testAtom->getType()==Atom::H)
            continue;
        
        /* Bond two atoms if their distance is less than the covalent bond threshold: */
        if(dist(*newAtom,*testAtom)<covalentBondThreshold)
        {
            /* Create backbone amide/carboxyl groups on-the-fly: */
            if(newAtom->getType()==Atom::H&&testAtom->getType()==Atom::N&&testAtom->getBackboneIndex()==0)
            {
                /* Build an amide group: */
                protein->amides.push_back(Dipole(testAtom,newAtom,1));
            }
            else if(newAtom->getType()==Atom::O&&testAtom->getType()==Atom::C&&testAtom->getBackboneIndex()==2)
            {
                /* Build a carboxyl group: */
                protein->carboxyls.push_back(Dipole(testAtom,newAtom,-1));
            }
            
            /* Bond the atoms: */
            bond(*newAtom,*testAtom);
            if(newAtom->getType()==Atom::H)
                break; // Only one bond per H atom, and to the closest "heavy" atom back in the chain
        }
    }
}

void Protein::ResidueCreator::finishProtein(void)
{
    // finalize the final residue
    finalizeResidue();
    
    /* Finish the last secondary structure: */
    if(currentSecondaryStructure!=0)
        currentSecondaryStructure->residueEnd=0;
    
    /* Associate dipoles with secondary structures: */
    SecondaryStructure* sPtr=0;
    for(std::vector<Dipole>::iterator nh=protein->amides.begin();nh!=protein->amides.end();++nh)
    {
        if(nh->getMajorAtom()->residue->secondaryStructure!=sPtr)
        {
            if(sPtr!=0)
                sPtr->amidesEnd=nh;
            sPtr=nh->getMajorAtom()->residue->secondaryStructure;
            sPtr->amidesBegin=nh;
        }
    }
    if(sPtr!=0)
        sPtr->amidesEnd=protein->amides.end();
    sPtr=0;
    for(std::vector<Dipole>::iterator co=protein->carboxyls.begin();co!=protein->carboxyls.end();++co)
    {
        if(co->getMajorAtom()->residue->secondaryStructure!=sPtr)
        {
            if(sPtr!=0)
                sPtr->carboxylsEnd=co;
            sPtr=co->getMajorAtom()->residue->secondaryStructure;
            sPtr->carboxylsBegin=co;
        }
    }
    if(sPtr!=0)
        sPtr->carboxylsEnd=protein->carboxyls.end();
    
    /* Create backbone bond iterators: */
    Residue* rPtr=protein->residues;
    std::vector<ChainAtom*>::iterator atomIt1,atomIt2;
    atomIt1=atomIt2=rPtr->backbone.begin();
    ++atomIt2;
    protein->beginIterator=BackboneIterator(atomIt1,atomIt2);
    while(rPtr->succ!=0)
        rPtr=rPtr->succ;
    atomIt1=atomIt2=rPtr->backbone.end();
    --atomIt1;
    protein->endIterator=BackboneIterator(atomIt1,atomIt2);
    
    /* Create array of pointers to atoms: */
    protein->atomPointers=new ChainAtom*[protein->numAtoms];
    int index=0;
    for(ChainAtom* aPtr=protein->atoms;aPtr!=0;aPtr=aPtr->succ,++index)
        protein->atomPointers[index]=aPtr;
}

/************************************************
Static elements of class Protein::ResidueCreator:
************************************************/

const char atomNamesInACE[6][4]={"CAY","HY1","HY2","HY3"," CY"," OY"};
const char atomNamesInNME[6][4]={" NT","HNT","CAT","HT1","HT2","HT3"};

int Protein::ResidueCreator::padHydrogenAtoms(int residueIndex, int numAtoms, const Position& atomPos, int type)
{
	 double atomPosition[3] = { atomPos[0], atomPos[1], atomPos[2] };
 	 #ifdef _VERBOSE_
	 printf("End Point %f %f %f \n", atomPosition[0],atomPosition[1], atomPosition[2]);
	 #endif

	 if (type==0) 
     	newResidue("ACE",--residueIndex);
     else 
     	newResidue("NME",++residueIndex);

     /* atom index starts with 1 */
     numAtoms++;
	 int index, i;
	 double offsetPos;
     for( index=0; index<6; index++) {
     	char elementName[2];
     	char* atomNamePtr = new char[4];
     
     	/* Parse the element name and offset structure position: */
    	if (type==0)
     		strcpy(atomNamePtr,atomNamesInACE[index]);
     	else
     		strcpy(atomNamePtr,atomNamesInNME[index]);
     
     	while(isspace(*atomNamePtr))
     		++atomNamePtr;
     	elementName[0]=*atomNamePtr;
     	++atomNamePtr;
     	elementName[1]='\0';
	 	
		/* Offset the first/last atom position as the new coord. of ACE/NME */
		i = (int)fmod(index,3);
		offsetPos = (double)(i+index+1)*0.1;
		atomPosition[i]+= offsetPos;
 	 	#ifdef _VERBOSE_
 		printf("Final Point %f %f %f \n", atomPosition[0],atomPosition[1], atomPosition[2]);
	 	#endif

     	addAtom(elementName,numAtoms+index,Position(atomPosition),atomNamePtr);
	}
	return index;
}

/******************************************
Methods of class Protein::BackboneIterator:
******************************************/

int Protein::BackboneIterator::getResidueIndex(void) const
{
    int result=-1;
    for(const Residue* rPtr=(*atom1)->residue;rPtr!=0;rPtr=rPtr->pred)
        ++result;
    return result;
}

/*******************************************
Methods of class Protein::StructureSelector:
*******************************************/

Protein::StructureSelector::StructureSelector(Protein* sProtein,SecondaryStructure* sStructure,int sStructureIndex)
    :protein(sProtein),structure(sStructure),structureIndex(sStructureIndex),
     firstResidueIndex(-1),numResidues(0)
{
    if(isValid())
    {
        /* Cache first residue index and number of residues: */
        firstResidueIndex=structure->getFirstResidueIndex();
        numResidues=structure->getNumResidues();
        
        /* Create iterators for beginning and end of backbone: */
        std::vector<ChainAtom*>::iterator atomIt1,atomIt2;
        
        /* Create the begin iterator out of the first residue's first two backbone atoms: */
        atomIt1=atomIt2=structure->residueBegin->backbone.begin();
        ++atomIt2;
        begin=BackboneIterator(atomIt1,atomIt2);
        
        /* Find the last residue in the structure: */
        Residue* rPtr;
        for(rPtr=structure->residueBegin;rPtr->succ!=structure->residueEnd;rPtr=rPtr->succ)
            ;
        
        /* Create an iterator out of the last residue's first two backbone atoms: */
        atomIt1=atomIt2=rPtr->backbone.begin();
        ++atomIt2;
        end=BackboneIterator(atomIt1,atomIt2);
        
        /* Increment the end iterator twice to point past the last residue's end: */
        ++end;
        ++end;
    }
}

int Protein::StructureSelector::getStructureTypeIndex(void) const
{
    int result=0;
    for(const SecondaryStructure* sPtr=protein->secondaryStructures;sPtr!=structure;sPtr=sPtr->succ)
        if(sPtr->getStructureType()==structure->getStructureType())
            ++result;
    return result;
}

DragBox Protein::StructureSelector::createDragBox(void) const
{
    /* Calculate centroid of selected structure: */
    Geometry::AffineCombiner<double,3> centroidCombiner;
    int numAtoms=0;
    for(const Residue* rPtr=structure->residueBegin;rPtr!=structure->residueEnd;rPtr=rPtr->succ)
    {
        #if 0
        for(const ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ,++numAtoms)
            centroidCombiner.addPoint(aPtr->getPosition());
        #else
        for(std::vector<ChainAtom*>::const_iterator aIt=rPtr->backbone.begin();aIt!=rPtr->backbone.end();++aIt,++numAtoms)
            centroidCombiner.addPoint((*aIt)->getPosition());
        #endif
    }
    Geometry::Point<double,3> centroid=centroidCombiner.getPoint();
    
    /* Create variance/covariance matrix for the selected structure: */
    Geometry::Matrix<double,3,3> cov=Geometry::Matrix<double,3,3>::zero;
    for(const Residue* rPtr=structure->residueBegin;rPtr!=structure->residueEnd;rPtr=rPtr->succ)
    {
        #if 0
        for(const ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
        {
            Geometry::Vector<double,3> d=aPtr->getPosition()-centroid;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    cov(i,j)+=d[i]*d[j];
        }
        #else
        for(std::vector<ChainAtom*>::const_iterator aIt=rPtr->backbone.begin();aIt!=rPtr->backbone.end();++aIt)
        {
            Geometry::Vector<double,3> d=(*aIt)->getPosition()-centroid;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    cov(i,j)+=d[i]*d[j];
        }
        #endif
    }
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            cov(i,j)/=double(numAtoms-1);
    
    /* Calculate the variance/covariance matrix's major eigenvector: */
    Geometry::Vector<double,3> a[3];
    a[0]=Geometry::Vector<double,3>(1.0,0.0,0.0);
    for(int i=0;i<10000;++i)
    {
        a[0]=Geometry::Vector<double,3>(cov*a[0]);
        a[0].normalize();
    }
    
    const Residue* rPtr=structure->residueBegin;
    while(rPtr!=structure->residueEnd&&rPtr->type==Residue::UNK) // Stay away from unknown residue types
        rPtr=rPtr->succ;
    if(rPtr!=structure->residueEnd) // There are known residues in the structure
    {
        /* Second axis is cross product of first and last backbone bond directions in structure: */
        Vector d1=rPtr->backbone[1]->getPosition()-rPtr->backbone[0]->getPosition();
        d1.normalize();
        while(rPtr->succ!=structure->residueEnd)
            rPtr=rPtr->succ;
        while(rPtr->type==Residue::UNK) // Move away from any unknown residues
            rPtr=rPtr->pred;
        Vector d2=rPtr->backbone[2]->getPosition()-rPtr->backbone[1]->getPosition();
        d2.normalize();
        a[1]=Geometry::cross(d2,d1);
        a[1]-=(a[0]*a[1])*a[0];
        double a1Len=Geometry::mag(a[1]);
        if(a1Len>1.0e-4)
            a[1]/=a1Len;
        else
            a[1]=Geometry::normalize(Geometry::normal(a[0]));
    }
    else
    {
        /* Second axis is just any axis orthogonal to the first one: */
        a[1]=Geometry::normalize(Geometry::normal(a[0]));
    }
    a[2]=Geometry::normalize(Geometry::cross(a[0],a[1]));
    
    /* Calculate the structure's extents along the three axes: */
    double s[3];
    for(int i=0;i<3;++i)
    {
        double min=Math::Constants<double>::max;
        double max=Math::Constants<double>::min;
        for(const Residue* rPtr=structure->residueBegin;rPtr!=structure->residueEnd;rPtr=rPtr->succ)
        {
            #if 0
            for(const ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
            {
                double d=a[i]*aPtr->getPosition();
                if(min>d)
                    min=d;
                if(max<d)
                    max=d;
            }
            #else
            for(std::vector<ChainAtom*>::const_iterator aIt=rPtr->backbone.begin();aIt!=rPtr->backbone.end();++aIt)
            {
                double d=a[i]*(*aIt)->getPosition();
                if(min>d)
                    min=d;
                if(max<d)
                    max=d;
            }
            #endif
        }
        s[i]=max-min;
        
        /* Adjust the centroid: */
        centroid+=(0.5*(min+max)-centroid*a[i])*a[i];
    }
    
    /* Return a drag box: */
    DragBox result;
    result.setCenter(centroid);
    for(int i=0;i<3;++i)
    {
        result.setAxis(i,a[i]);
        result.setSize(i,s[i]+2.0);
    }
    result.setEdgeRadius(double(bondRadius));
    return result;
}

/*************************************
Methods of class Protein::BondRotator:
*************************************/

void Protein::BondRotator::rotate(Scalar angle)
{
    /* Calculate rotation transformations for both parts of the molecule: */
    Geometry::AffineTransformation<Scalar,3> rotation[2];
    for(int i=0;i<2;++i)
    {
        rotation[i]=Geometry::AffineTransformation<Scalar,3>::identity;
        rotation[i]*=Geometry::AffineTransformation<Scalar,3>::translate(getStart()-Point::origin);
        rotation[i]*=Geometry::AffineTransformation<Scalar,3>::rotate(Geometry::Rotation<Scalar,3>::rotateAxis(getAxis(),i==1?Scalar(0.5)*angle:-Scalar(0.5)*angle));
        rotation[i]*=Geometry::AffineTransformation<Scalar,3>::translate(Point::origin-getStart());
    }
    
    /* Perform a breadth-first traversal on all atoms: */
    ++protein->rotationMarker;
    std::deque<MoleculeTraversal> traversalQueue;
    (*atom1)->rotationMarker=protein->rotationMarker;
    traversalQueue.push_back(MoleculeTraversal(*atom1,0));
    (*atom2)->rotationMarker=protein->rotationMarker;
    traversalQueue.push_back(MoleculeTraversal(*atom2,1));
    while(!traversalQueue.empty())
    {
        /* Pick the first element from the queue: */
        MoleculeTraversal mt=traversalQueue.front();
        traversalQueue.pop_front();
        
        /* Rotate it: */
        mt.atom->setPosition(rotation[mt.rotationIndex].transform(mt.atom->getPosition()));
        
        /* Put all its neighbours into the queue: */
        for(std::vector<Atom*>::iterator aIt=mt.atom->getBonds().begin();aIt!=mt.atom->getBonds().end();++aIt)
        {
            ChainAtom* aPtr=static_cast<ChainAtom*>(*aIt);
            if(aPtr->rotationMarker!=protein->rotationMarker)
            {
                aPtr->rotationMarker=protein->rotationMarker;
                traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex));
            }
        }
    }
}

Scalar Protein::BondRotator::getDihedralAngle(void) const
{
    /* Find the atoms before atom1 and after atom2 in the backbone chain: */
    assert((*atom1)->residue==(*atom2)->residue); // Both atoms must be in the same residue (peptide bonds don't rotate)

	const Residue* res=(*atom1)->residue;
    std::vector<ChainAtom*>::const_iterator atomIt;
    
    /* Find the atom before atom1 by going backwards: */
    const ChainAtom* atom0;
    atomIt=atom1;
    if(atomIt!=res->backbone.begin())
    {
        --atomIt;
        atom0=*atomIt;
    }
    else if(res->pred!=0)
        atom0=res->pred->backbone.back();
    else
        return Scalar(0); // Ugly hack: The first phi angle is set to zero
    
    /* Find the atom behind atom2 by going forward: */
    const ChainAtom* atom3;
    atomIt=atom2;
    ++atomIt;
    if(atomIt!=res->backbone.end())
        atom3=*atomIt;
    else if(res->succ!=0)
        atom3=res->succ->backbone.front();
    else
        return Scalar(0); // Ugly hack: The last psi angle is set to zero

	/* Now calculate the normal vectors for the planes defined by the two atom triples: */
    Vector d1=(*atom1)->getPosition()-atom0->getPosition();
    Vector d2=(*atom2)->getPosition()-(*atom1)->getPosition();
    Vector d3=atom3->getPosition()-(*atom2)->getPosition();
    Vector n1=Geometry::cross(d1,d2);
    Vector n2=Geometry::cross(d2,d3);
    
    /* The dihedral angle is the angle between the two normals: */
    Scalar angle=Math::acos((n1*n2)/(Geometry::mag(n1)*Geometry::mag(n2)));
    
    /* Determine the angle's sign by comparing the cross product of n1 and n2 with the bond axis: */
    if(Geometry::cross(n1,n2)*d2<Scalar(0))
        angle=-angle;
    
    return angle;
}

/********************************
Static elements of class Protein:
********************************/

Scalar Protein::bondRadius=Scalar(0.25);

/************************
Methods of class Protein:
************************/

Protein::Protein(void)
    :numAtoms(0),atoms(0),residues(0),secondaryStructures(0),rotationMarker(0),atomPointers(0)
{
}

Protein::~Protein(void)
{
    /* Delete all secondary structures: */
    while(secondaryStructures!=0)
    {
        SecondaryStructure* succ=secondaryStructures->succ;
        delete secondaryStructures;
        secondaryStructures=succ;
    }
    
    /* Delete all residues: */
    while(residues!=0)
    {
        Residue* succ=residues->succ;
        delete residues;
        residues=succ;
    }
    
    /* Delete all atoms: */
    while(atoms!=0)
    {
        ChainAtom* succ=atoms->succ;
        delete atoms;
        atoms=succ;
    }
    
    /* Delete the atom pointers: */
    delete[] atomPointers;
}

Position Protein::calcCentroid(void) const
{
    /* Iterate through all atoms to calculate the protein's centroid: */
    Point::AffineCombiner centroidCombiner;
    for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
        centroidCombiner.addPoint(aPtr->getPosition());
    
    return centroidCombiner.getPoint();
}

double Protein::calcRadius(void) const
{
    /* Iterate through all atoms to calculate the protein's centroid: */
    Point::AffineCombiner centroidCombiner;
    for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
        centroidCombiner.addPoint(aPtr->getPosition());
    Point centroid=centroidCombiner.getPoint();
    
    /* Calculate the maximum distance of any atom from the centroid: */
    double radius2=0.0;
    for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
    {
        double d2=double(Geometry::sqrDist(Point(aPtr->getPosition()),centroid));
        if(radius2<d2)
            radius2=d2;
    }
    
    return Math::sqrt(radius2);
}

int Protein::getNumResidues(void) const
{
    int result=0;
    for(const Residue* rPtr=residues;rPtr!=0;rPtr=rPtr->succ)
        ++result;
    return result;
}

std::pair<int,int> Protein::getResidueIndexRange(void) const
{
    const Residue* rPtr=residues;
    int rangeMin,rangeMax;
    rangeMin=rangeMax=rPtr->residueIndex;
    for(rPtr=rPtr->succ;rPtr!=0;rPtr=rPtr->succ)
    {
        if(rangeMin>rPtr->residueIndex)
            rangeMin=rPtr->residueIndex;
        else if(rangeMax<rPtr->residueIndex)
            rangeMax=rPtr->residueIndex;
    }
    
    return std::pair<int,int>(rangeMin,rangeMax);
}

int Protein::getNumStructures(void) const
{
    int result=0;
    for(const SecondaryStructure* sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
        ++result;
    return result;
}

Protein::ConstAtomIterator Protein::pickAtom(const Ray& ray) const
{
    Scalar lambdaMin=Math::Constants<Scalar>::max;
    const ChainAtom* selected=0;
    
    /* Iterate through all atoms in the protein and intersect the query ray with each: */
    for(const ChainAtom* aPtr=atoms;aPtr!=0;aPtr=aPtr->succ)
    {
        /* Intersect the query ray with the atom: */
        Scalar lambda=aPtr->intersectRay(ray);
        if(lambdaMin>lambda) // Intersection is the closest so far
        {
            selected=aPtr;
            lambdaMin=lambda;
        }
    }
    
    return ConstAtomIterator(selected);
}

const Protein::Residue* Protein::pickResidue(int pdbResidueIndex) const
{
    const Residue* rPtr;
    for(rPtr=residues;rPtr!=0&&rPtr->residueIndex!=pdbResidueIndex;rPtr=rPtr->succ)
        ;
    return rPtr;
}

Protein::Residue* Protein::pickResidue(int pdbResidueIndex)
{
    Residue* rPtr;
    for(rPtr=residues;rPtr!=0&&rPtr->residueIndex!=pdbResidueIndex;rPtr=rPtr->succ)
        ;
    return rPtr;
}

Protein::Residue* Protein::pickResidue(const Ray& ray)
{
    Scalar lambdaMin=Math::Constants<Scalar>::max;
    Residue* selected=0;
    
    /* Iterate through all backbone bonds in the protein and intersect the query ray with each: */
    for(BackboneIterator bbIt=beginIterator;bbIt!=endIterator;++bbIt)
    {
        /* Construct a cylinder around the bond and intersect the query ray with it: */
        Geometry::Cylinder<Scalar,3> bond(bbIt.getAtom1()->getPosition(),bbIt.getAtom2()->getPosition(),bondRadius);
        Scalar lambda=bond.intersectRay(ray).getParameter();
        if(lambdaMin>lambda)
        {
            lambdaMin=lambda;
            if(bbIt.isPeptideBond())
            {
                /* Check which side of the peptide bond was selected: */
                Point intersection=ray(lambda);
                if(Geometry::sqrDist(bbIt.getAtom1()->getPosition(),intersection)<=Geometry::sqrDist(bbIt.getAtom2()->getPosition(),intersection))
                    selected=bbIt.getAtom1()->getResidue();
                else
                    selected=bbIt.getAtom2()->getResidue();
            }
            else
                selected=bbIt.getResidue();
        }
    }
    
    return selected;
}

int Protein::getResidueIndex(const Protein::Residue* rPtr) const
{
    int result=0;
    const Residue* rPtr2;
    for(rPtr2=residues;rPtr2!=0&&rPtr2!=rPtr;rPtr2=rPtr2->succ,++result)
        ;
    if(rPtr2==0)
        result=-1;
    return result;
}

void Protein::changeResidueStructureType(Protein::Residue* rPtr,Protein::SecondaryStructure::StructureType newType)
{
    SecondaryStructure* sPtr=rPtr->secondaryStructure;
    if(newType!=sPtr->structureType)
    {
        if(rPtr!=sPtr->residueBegin&&rPtr->succ!=sPtr->residueEnd) // Residue is in the middle of structure
        {
            /* Split current structure into three: */
            SecondaryStructure* leftStructure=new SecondaryStructure(sPtr->structureType);
            leftStructure->residueBegin=sPtr->residueBegin;
            leftStructure->residueEnd=sPtr->residueBegin=rPtr;
            SecondaryStructure* rightStructure=new SecondaryStructure(sPtr->structureType);
            rightStructure->residueEnd=sPtr->residueEnd;
            rightStructure->residueBegin=sPtr->residueEnd=rPtr->succ;
            sPtr->structureType=newType;
            leftStructure->pred=sPtr->pred;
            if(leftStructure->pred!=0)
                leftStructure->pred->succ=leftStructure;
            else
                secondaryStructures=leftStructure;
            leftStructure->succ=sPtr;
            sPtr->pred=leftStructure;
            rightStructure->succ=sPtr->succ;
            if(rightStructure->succ!=0)
                rightStructure->succ->pred=rightStructure;
            rightStructure->pred=sPtr;
            sPtr->succ=rightStructure;
        }
        else if(rPtr==sPtr->residueBegin&&rPtr->succ!=sPtr->residueEnd)
        {
            /* Change first residue in structure: */
            if(sPtr->pred==0||sPtr->pred->structureType!=newType)
            {
                /* Insert new structure before sPtr: */
                SecondaryStructure* leftStructure=new SecondaryStructure(newType);
                leftStructure->residueBegin=sPtr->residueBegin;
                leftStructure->residueEnd=sPtr->residueBegin=rPtr->succ;
                leftStructure->pred=sPtr->pred;
                if(sPtr->pred!=0)
                    sPtr->pred->succ=leftStructure;
                else
                    secondaryStructures=leftStructure;
                leftStructure->succ=sPtr;
                sPtr->pred=leftStructure;
            }
            else
            {
                /* Move separation between sPtr->pred and sPtr: */
                sPtr->pred->residueEnd=sPtr->residueBegin=rPtr->succ;
            }
        }
        else if(rPtr!=sPtr->residueBegin&&rPtr->succ==sPtr->residueEnd)
        {
            /* Change last residue in structure: */
            if(sPtr->succ==0||sPtr->succ->structureType!=newType)
            {
                /* Insert new structure after sPtr: */
                SecondaryStructure* rightStructure=new SecondaryStructure(newType);
                rightStructure->residueEnd=sPtr->residueEnd;
                rightStructure->residueBegin=sPtr->residueEnd=rPtr;
                rightStructure->succ=sPtr->succ;
                if(sPtr->succ!=0)
                    sPtr->succ->pred=rightStructure;
                rightStructure->pred=sPtr;
                sPtr->succ=rightStructure;
            }
            else
            {
                /* Move separation between sPtr and sPtr->succ: */
                sPtr->residueEnd=sPtr->succ->residueBegin=rPtr;
            }
        }
        else
        {
            /* Change entire structure: */
            sPtr->structureType=newType;
            SecondaryStructure* left=sPtr->pred;
            if(left!=0&&left->structureType==newType)
            {
                /* Merge with structure to the left: */
                sPtr->residueBegin=left->residueBegin;
                sPtr->pred=left->pred;
                if(left->pred!=0)
                    left->pred->succ=sPtr;
                else
                    secondaryStructures=sPtr;
                delete left;
            }
            SecondaryStructure* right=sPtr->succ;
            if(right!=0&&right->structureType==newType)
            {
                /* Merge with structure to the right: */
                sPtr->residueEnd=right->residueEnd;
                sPtr->succ=right->succ;
                if(right->succ!=0)
                    right->succ->pred=sPtr;
                delete right;
            }
        }
        
        /* Re-associate residues with secondary structures: */
        for(sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
            for(Residue* rPtr=sPtr->residueBegin;rPtr!=sPtr->residueEnd;rPtr=rPtr->succ)
                rPtr->secondaryStructure=sPtr;
        
        /* Re-associate dipoles with secondary structures: */
        std::vector<Dipole>::iterator nh=amides.begin();
        for(sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
        {
            sPtr->amidesBegin=nh;
            while(nh!=amides.end()&&nh->getMajorAtom()->residue->secondaryStructure==sPtr)
                ++nh;
            sPtr->amidesEnd=nh;
        }
        std::vector<Dipole>::iterator co=carboxyls.begin();
        for(sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
        {
            sPtr->carboxylsBegin=co;
            while(co!=carboxyls.end()&&co->getMajorAtom()->residue->secondaryStructure==sPtr)
                ++co;
            sPtr->carboxylsEnd=co;
        }
    }
}

Protein::StructureSelector Protein::pickStructure(int index)
{
    SecondaryStructure* sPtr;
    int i;
    for(sPtr=secondaryStructures,i=0;sPtr!=0&&i!=index;sPtr=sPtr->succ,++i)
        ;
    return StructureSelector(this,sPtr,i);
}

Protein::StructureSelector Protein::pickStructure(Protein::SecondaryStructure::StructureType structureType,int index)
{
    SecondaryStructure* sPtr=secondaryStructures;
    int i=0;
    while(true)
    {
        while(sPtr!=0&&sPtr->getStructureType()!=structureType)
        {
            sPtr=sPtr->succ;
            ++i;
        }
        --index;
        if(sPtr==0||index<0)
            break;
        sPtr=sPtr->succ;
        ++i;
    }
    return StructureSelector(this,sPtr,i);
}

Protein::StructureSelector Protein::pickStructure(const Protein::Residue* residue)
{
    int i=0;
    for(SecondaryStructure* sPtr=secondaryStructures;sPtr!=residue->getSecondaryStructure();++i,sPtr=sPtr->succ)
        ;
    return StructureSelector(this,residue->getSecondaryStructure(),i);
}

Protein::StructureSelector Protein::pickStructure(const Point& p)
{
    SecondaryStructure* selected=0;
    int selectedIndex=-1;
    Scalar distance2Min=Math::Constants<Scalar>::max;
    
    /* Iterate through all atoms in the protein and check each against the query point: */
    int i=0;
    for(SecondaryStructure* sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ,++i)
    {
        for(Residue* rPtr=sPtr->residueBegin;rPtr!=sPtr->residueEnd;rPtr=rPtr->succ)
        {
            for(ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
            {
                /* Intersect the query ray with the covalent bond sphere around the atom: */
                Scalar distance2=Geometry::sqrDist(aPtr->getPosition(),p);
                if(distance2<=Math::sqr(aPtr->getCovalentRadius()*Scalar(1.5))&&distance2<distance2Min)
                {
                    distance2Min=distance2;
                    selected=sPtr;
                    selectedIndex=i;
                }
            }
        }
    }
    
    return StructureSelector(this,selected,selectedIndex);
}

Protein::StructureSelector Protein::pickStructure(const Ray& ray)
{
    SecondaryStructure* selected=0;
    int selectedIndex=-1;
    Scalar lambdaMin=Math::Constants<Scalar>::max;
    
    /* Iterate through all atoms in the protein and intersect the query ray with each: */
    int i=0;
    for(SecondaryStructure* sPtr=secondaryStructures;sPtr!=0;sPtr=sPtr->succ,++i)
    {
        for(Residue* rPtr=sPtr->residueBegin;rPtr!=sPtr->residueEnd;rPtr=rPtr->succ)
        {
            for(ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
            {
                /* Intersect the query ray with the atom: */
                Scalar lambda=aPtr->intersectRay(ray);
                if(lambdaMin>lambda) // Intersection is the closest so far: */
                {
                    lambdaMin=lambda;
                    selected=sPtr;
                    selectedIndex=i;
                }
            }
        }
    }
    
    return StructureSelector(this,selected,selectedIndex);
}

Protein::BondRotator Protein::pickBond(const Point& p)
{
    Scalar radius2=Math::sqr(bondRadius);
    
    /* Iterate through all backbone bonds and find the one closest to the given point: */
    std::vector<ChainAtom*>::iterator atom1,atom2;
    BondRotator::IntersectionType intersectionType=BondRotator::NONE;
    Scalar minDist2=Math::Constants<Scalar>::max;
    
    /***************************************************************************
    Due to the two-level hierarchy of storing atoms, iterating through all
    backbone bonds is slightly more complicated than one might initially expect.
    However, the fact that every third bond is a rigid peptide bond falls right
    out of the data structure, so that is a benefit at least.
    ***************************************************************************/
    
    Residue* rPtr=residues;
    std::vector<ChainAtom*>::iterator atomIt=rPtr->backbone.begin();
    while(true)
    {
        /* Pick two neighbouring atoms: */
        bool peptideBond=false;
        std::vector<ChainAtom*>::iterator a1=atomIt;
        ++atomIt;
        if(atomIt==rPtr->backbone.end()) // The second atom is in another residue
        {
            rPtr=rPtr->succ;
            if(rPtr==0) // Reached the end of the backbone
                break;
            atomIt=rPtr->backbone.begin();
            peptideBond=true; // This must therefore be a peptide bond
        }
        std::vector<ChainAtom*>::iterator a2=atomIt;
        if(!peptideBond) // Peptide bonds are rigid and cannot be picked!
        {
            /* Construct an upright cylinder connecting the two atoms, and test it against the point: */
            Point p1=Point((*a1)->getPosition());
            Point p2=Point((*a2)->getPosition());
            Vector d=p2-p1;
            Scalar dLen=Scalar(Geometry::mag(d));
            d/=dLen;
            Vector pp1=p-p1;
            Scalar beta=pp1*d;
            Scalar dist2=Geometry::sqr(pp1-beta*d);
            if(dist2<minDist2&&dist2<=radius2&&beta>=Scalar(0)&&beta<=dLen)
            {
                atom1=a1;
                atom2=a2;
                intersectionType=BondRotator::POINT_INSIDE;
                minDist2=dist2;
            }
        }
    }
    
    return BondRotator(this,atom1,atom2,intersectionType,p);
}

Protein::BondRotator Protein::pickBond(const Ray& ray)
{
    /* Iterate through all backbone bonds and find the closest intersection with the given ray: */
    std::vector<ChainAtom*>::iterator atom1,atom2;
    BondRotator::IntersectionType intersectionType=BondRotator::NONE;
    Scalar lambdaMin=Math::Constants<Scalar>::max;
    Point intersection;
    
    /***************************************************************************
    Due to the two-level hierarchy of storing atoms, iterating through all
    backbone bonds is slightly more complicated than one might initially expect.
    However, the fact that every third bond is a rigid peptide bond falls right
    out of the data structure, so that is a benefit at least.
    ***************************************************************************/
    
    Residue* rPtr=residues;
    std::vector<ChainAtom*>::iterator atomIt=rPtr->backbone.begin();
    while(true)
    {
        /* Pick two neighbouring atoms: */
        bool peptideBond=false;
        std::vector<ChainAtom*>::iterator a1=atomIt;
        ++atomIt;
        if(atomIt==rPtr->backbone.end()) // The second atom is in another residue
        {
            rPtr=rPtr->succ;
            if(rPtr==0) // Reached the end of the backbone
                break;
            atomIt=rPtr->backbone.begin();
            peptideBond=true; // This must therefore be a peptide bond
        }
        std::vector<ChainAtom*>::iterator a2=atomIt;
        if(!peptideBond) // Peptide bonds are rigid and cannot be picked!
        {
            /* Construct an upright cylinder connecting the two atoms, and intersect it with the ray: */
            Geometry::Cylinder<Scalar,3> bond((*a1)->getPosition(),(*a2)->getPosition(),bondRadius);
            Geometry::Cylinder<Scalar,3>::HitResult result=bond.intersectRay(ray);
            if(result.isValid()&&result.getParameter()<lambdaMin)
            {
                switch(result.getPart())
                {
                    case Geometry::Cylinder<Scalar,3>::HitResult::MANTEL:
                        intersectionType=BondRotator::SIDE;
                        break;
                    
                    case Geometry::Cylinder<Scalar,3>::HitResult::BOTTOMCAP:
                        intersectionType=BondRotator::BOTTOM;
                        break;
                    
                    case Geometry::Cylinder<Scalar,3>::HitResult::TOPCAP:
                        intersectionType=BondRotator::TOP;
                        break;

                    case Geometry::Cylinder<Scalar,3>::HitResult::INVALID_PART:
			break;

                    default:
                        /* Do nothing */;
                        break;
                }
                lambdaMin=result.getParameter();
                intersection=ray(lambdaMin);
            }
        }
    }
    
    return BondRotator(this,atom1,atom2,intersectionType,intersection);
}

void Protein::getDihedralAngles(int startResidue,int numResidues,Scalar phis[],Scalar psis[]) const
{
    /* Iterate to the first queried residue: */
    const Residue* rPtr=residues;
    for(int i=0;i<startResidue;++i)
        rPtr=rPtr->succ;
    
    /* Retrieve all residues' phi and psi angles: */
    for(int i=0;i<numResidues;++i,rPtr=rPtr->succ)
    {
        if(rPtr->type!=Residue::UNK && rPtr->type!=Residue::ACE && rPtr->type!=Residue::NME ) // Don't assume anything about unknown residue types
        {
			/* Get pointers to the five atoms defining the residue's three planes: */
            const ChainAtom* atom0=rPtr->pred!=0?rPtr->pred->backbone[rPtr->pred->backbone.size()-1]:rPtr->backbone[0];
            const ChainAtom* atom1=rPtr->backbone[0];
            const ChainAtom* atom2=rPtr->backbone[1];
            const ChainAtom* atom3=rPtr->backbone[2];
            const ChainAtom* atom4=rPtr->succ!=0?rPtr->succ->backbone[0]:rPtr->backbone[2];

            /* Calculate the planes' normal vectors: */
            Vector d01=atom1->getPosition()-atom0->getPosition();
            Vector d12=atom2->getPosition()-atom1->getPosition();
            Vector d23=atom3->getPosition()-atom2->getPosition();
            Vector d34=atom4->getPosition()-atom3->getPosition();
            Vector n1=Geometry::cross(d01,d12);
            double n1Len=Geometry::mag(n1);
            Vector n2=Geometry::cross(d12,d23);
            double n2Len=Geometry::mag(n2);
            Vector n3=Geometry::cross(d23,d34);
            double n3Len=Geometry::mag(n3);

            /* Phi is the angle between n1 and n2: */
            if(rPtr->pred!=0)
            {
                phis[i]=Math::acos((n1*n2)/Scalar(n1Len*n2Len));
                if(Geometry::cross(n1,n2)*d12<Scalar(0)) // Angle is positive if clockwise along main bond
                    phis[i]=-phis[i];
            }
            else
                phis[i]=0.0; // Force first phi angle in protein to zero

            /* Psi is the angle between n2 and n3: */
            if(rPtr->succ!=0)
            {
                psis[i]=Math::acos((n2*n3)/Scalar(n2Len*n3Len));
                if(Geometry::cross(n2,n3)*d23<Scalar(0)) // Angle is positive if clockwise along main bond
                    psis[i]=-psis[i];
            }
            else
                psis[i]=0.0; // Force last psi angle in protein to zero
        }
        else
            phis[i]=psis[i]=0.0; // Force dihedral angles of unknown residues to zero
    }
}

void Protein::changeDihedralAngles(int startResidue,int numResidues,const Scalar deltaPhis[],const Scalar deltaPsis[],int direction,bool stopAtLastResidue)
{
    /* Calculate cumulative transformation matrices along the backbone: */
    std::vector<Geometry::AffineTransformation<Scalar,3> > rotations;
    rotations.reserve(numResidues*2);
    if(direction>0)
    {
        /* Iterate to the first affected residue: */
        Residue* startResiduePtr=residues;
        for(int i=0;i<startResidue;++i)
            startResiduePtr=startResiduePtr->succ;
        
        /* Build rotation matrices: */
        Geometry::AffineTransformation<Scalar,3> currentRotation=Geometry::AffineTransformation<Scalar,3>::identity;
        Residue* rPtr=startResiduePtr;
        for(int i=0;i<numResidues;++i,rPtr=rPtr->succ)
        {
            if(rPtr->type==Residue::UNK) // Again, punt on unknown residue types
            {
                /* Create as many dummy rotations as there are backbone bonds inside this residue: */
                for(int j=1;j<rPtr->backbone.size();++j)
                    rotations.push_back(currentRotation);
            }
            else
            {
                /* Get atom positions: */
                Point p[3];
                for(int j=0;j<3;++j)
                    p[j]=rPtr->backbone[j]->getPosition();

                /* Concatenate phi rotation: */
                if(rPtr->type!=Residue::PRO) // Never rotate proline around its phi axis - it forms a rigid ring structure there
                {
                    currentRotation*=Geometry::AffineTransformation<Scalar,3>::translate(p[0]-Point::origin);
                    currentRotation*=Geometry::AffineTransformation<Scalar,3>::rotate(Geometry::Rotation<Scalar,3>::rotateAxis(p[1]-p[0],deltaPhis[i]));
                    currentRotation*=Geometry::AffineTransformation<Scalar,3>::translate(Point::origin-p[0]);
                }
                rotations.push_back(currentRotation);

                /* Concatenate psi rotation: */
                currentRotation*=Geometry::AffineTransformation<Scalar,3>::translate(p[1]-Point::origin);
                currentRotation*=Geometry::AffineTransformation<Scalar,3>::rotate(Geometry::Rotation<Scalar,3>::rotateAxis(p[2]-p[1],deltaPsis[i]));
                currentRotation*=Geometry::AffineTransformation<Scalar,3>::translate(Point::origin-p[1]);
                rotations.push_back(currentRotation);
            }
        }
        
        /* Perform a breadth-first traversal on all atoms to the right of the start residue: */
        int lastRotation=rotations.size()-1;
        ++rotationMarker;
        std::deque<MoleculeTraversal> traversalQueue;
        startResiduePtr->backbone[0]->rotationMarker=rotationMarker;
        startResiduePtr->backbone[1]->rotationMarker=rotationMarker;
        if(stopAtLastResidue&&rPtr!=0)
            rPtr->backbone[0]->rotationMarker=rotationMarker;
        traversalQueue.push_back(MoleculeTraversal(startResiduePtr->backbone[1],0));
        while(!traversalQueue.empty())
        {
            /* Pick the first element from the queue: */
            MoleculeTraversal mt=traversalQueue.front();
            traversalQueue.pop_front();

            /* Rotate it: */
            mt.atom->setPosition(rotations[mt.rotationIndex].transform(mt.atom->getPosition()));

            /* Put all its neighbours into the queue: */
            for(std::vector<Atom*>::iterator aIt=mt.atom->getBonds().begin();aIt!=mt.atom->getBonds().end();++aIt)
            {
                ChainAtom* aPtr=static_cast<ChainAtom*>(*aIt);
                if(aPtr->rotationMarker!=rotationMarker)
                {
                    aPtr->rotationMarker=rotationMarker;
                    if(mt.rotationIndex<lastRotation&&aPtr->backboneIndex>=1&&mt.atom->backboneIndex==aPtr->backboneIndex-1) // Is traversed bond a phi/psi backbone bond inside the residue sequence?
                        traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex+1));
                    else
                        traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex));
                }
            }
        }
    }
    else
    {
        /* Iterate to the last affected residue: */
        Residue* finishResiduePtr=residues;
        for(int i=0;i<startResidue+numResidues-1;++i)
            finishResiduePtr=finishResiduePtr->succ;
        
        /* Build rotation matrices: */
        Geometry::AffineTransformation<Scalar,3> currentRotation=Geometry::AffineTransformation<Scalar,3>::identity;
        Residue* rPtr=finishResiduePtr;
        for(int i=numResidues-1;i>=0;--i,rPtr=rPtr->pred)
        {
            if(rPtr->type==Residue::UNK) // Again, punt on unknown residue types
            {
                /* Create as many dummy rotations as there are backbone bonds inside this residue: */
                for(int j=1;j<rPtr->backbone.size();++j)
                    rotations.push_back(currentRotation);
            }
            else
            {
                /* Get atom positions: */
                Point p[3];
                for(int j=0;j<3;++j)
                    p[j]=rPtr->backbone[j]->getPosition();
                
                /* Concatenate psi rotation: */
                currentRotation*=Geometry::AffineTransformation<Scalar,3>::translate(p[2]-Point::origin);
                currentRotation*=Geometry::AffineTransformation<Scalar,3>::rotate(Geometry::Rotation<Scalar,3>::rotateAxis(p[1]-p[2],deltaPsis[i]));
                currentRotation*=Geometry::AffineTransformation<Scalar,3>::translate(Point::origin-p[2]);
                rotations.push_back(currentRotation);
                
                /* Concatenate phi rotation: */
                if(rPtr->type!=Residue::PRO) // Never rotate proline around its phi axis - it forms a rigid ring structure there
                {
                    currentRotation*=Geometry::AffineTransformation<Scalar,3>::translate(p[1]-Point::origin);
                    currentRotation*=Geometry::AffineTransformation<Scalar,3>::rotate(Geometry::Rotation<Scalar,3>::rotateAxis(p[0]-p[1],deltaPhis[i]));
                    currentRotation*=Geometry::AffineTransformation<Scalar,3>::translate(Point::origin-p[1]);
                }
                rotations.push_back(currentRotation);
            }
        }
        
        /* Perform a breadth-first traversal on all atoms to the left of the finish residue: */
        int lastRotation=rotations.size()-1;
        ++rotationMarker;
        std::deque<MoleculeTraversal> traversalQueue;
        finishResiduePtr->backbone[finishResiduePtr->backbone.size()-1]->rotationMarker=rotationMarker;
        finishResiduePtr->backbone[finishResiduePtr->backbone.size()-2]->rotationMarker=rotationMarker;
        if(stopAtLastResidue&&rPtr!=0)
            rPtr->backbone[0]->rotationMarker=rotationMarker;
        traversalQueue.push_back(MoleculeTraversal(finishResiduePtr->backbone[finishResiduePtr->backbone.size()-2],0));
        while(!traversalQueue.empty())
        {
            /* Pick the first element from the queue: */
            MoleculeTraversal mt=traversalQueue.front();
            traversalQueue.pop_front();

            /* Rotate it: */
            mt.atom->setPosition(rotations[mt.rotationIndex].transform(mt.atom->getPosition()));

            /* Put all its neighbours into the queue: */
            for(std::vector<Atom*>::iterator aIt=mt.atom->getBonds().begin();aIt!=mt.atom->getBonds().end();++aIt)
            {
                ChainAtom* aPtr=static_cast<ChainAtom*>(*aIt);
                if(aPtr->rotationMarker!=rotationMarker)
                {
                    aPtr->rotationMarker=rotationMarker;
                    if(mt.rotationIndex<lastRotation&&mt.atom->backboneIndex>=1&&aPtr->backboneIndex==mt.atom->backboneIndex-1) // Is traversed bond a phi/psi backbone bond inside the residue sequence?
                        traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex+1));
                    else
                        traversalQueue.push_back(MoleculeTraversal(aPtr,mt.rotationIndex));
                }
            }
        }
    }
}

void Protein::setDihedralAngles(int startResidue,int numResidues,const Scalar phis[],const Scalar psis[],int direction)
{
    /* Calculate the current dihedral angles: */
    Scalar* deltaPhis=new Scalar[numResidues];
    Scalar* deltaPsis=new Scalar[numResidues];
    getDihedralAngles(startResidue,numResidues,deltaPhis,deltaPsis);
    
    /* Calculate angle increments: */
    for(int i=0;i<numResidues;++i)
    {
        deltaPhis[i]=phis[i]-deltaPhis[i];
        deltaPsis[i]=psis[i]-deltaPsis[i];
    }
    
    /* Update the dihedral angles: */
    changeDihedralAngles(startResidue,numResidues,deltaPhis,deltaPsis,direction);
    
    /* Clean up: */
    delete[] deltaPhis;
    delete[] deltaPsis;
}

void Protein::setAtomPositions (const Scalar *x, const Scalar *y, const Scalar *z)
{
    ChainAtom *aPtr = atoms;
    for ( int i = 0; i < numAtoms && aPtr; ++i, aPtr = aPtr->succ )
    {
        Position p (x[i], y[i], z[i]);
        aPtr->setPosition (p);
    }
}


}   // end namespace MD
