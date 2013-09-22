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

#include <Geometry/AffineTransformation.h>

#include "UndoBuffer.h"

/*******************************************
Methods of class UndoBuffer::UndoQueueEntry:
*******************************************/

void UndoBuffer::UndoQueueEntry::set(MD::Protein* sProtein,const Geometry::OrthonormalTransformation<double,3>& sProteinTransformation,const std::list<DihedralAngleSequence> previousSequences[2])
{
    /* Copy parameters: */
    protein=sProtein;
    proteinTransformation=sProteinTransformation;
    sequences[0].clear();
    sequences[1].clear();
    
    /* Process both lists of sequences: */
    for(int i=0;i<2;++i)
    {
        /* Calculate and store dihedral angle changes for all left/right sequences: */
        for(std::list<DihedralAngleSequence>::const_iterator pIt=previousSequences[i].begin();pIt!=previousSequences[i].end();++pIt)
        {
            /* Add sequence to sequence list: */
            DihedralAngleSequence& deltas=*sequences[i].insert(sequences[i].end(),DihedralAngleSequence());

            /* Calculate (negative) angle changes: */
            deltas.resize(pIt->getFirstResidue(),pIt->getNumResidues());
            protein->getDihedralAngles(deltas.getFirstResidue(),deltas.getNumResidues(),deltas.phis,deltas.psis);
            for(int j=0;j<deltas.getNumResidues();++j)
            {
                deltas.phis[j]=pIt->phis[j]-deltas.phis[j];
                deltas.psis[j]=pIt->psis[j]-deltas.psis[j];
            }
        }
    }
}

void UndoBuffer::UndoQueueEntry::undo(void)
{
    /* Apply inverse protein transformation: */
    protein->transform(Geometry::AffineTransformation<MD::Scalar,3>(Geometry::invert(proteinTransformation)));
    
    /* Process both lists of sequences: */
    for(int i=0;i<2;++i)
    {
        /* Determine update direction: */
        int direction=i==0?-1:1;
        
        /* Process all left/right sequences: */
        for(std::list<DihedralAngleSequence>::const_iterator pIt=sequences[i].begin();pIt!=sequences[i].end();++pIt)
        {
            /* Apply negative angle changes: */
            protein->changeDihedralAngles(pIt->getFirstResidue(),pIt->getNumResidues(),pIt->phis,pIt->psis,direction);
        }
    }
}

void UndoBuffer::UndoQueueEntry::redo(void)
{
    /* Process both lists of sequences: */
    for(int i=0;i<2;++i)
    {
        /* Determine update direction: */
        int direction=i==0?-1:1;
        
        /* Process all left/right sequences: */
        for(std::list<DihedralAngleSequence>::const_iterator pIt=sequences[i].begin();pIt!=sequences[i].end();++pIt)
        {
            /* Negate angle deltas: */
            for(int j=0;j<pIt->getNumResidues();++j)
            {
                pIt->phis[j]=-pIt->phis[j];
                pIt->psis[j]=-pIt->psis[j];
            }
            
            /* Apply positive angle changes: */
            protein->changeDihedralAngles(pIt->getFirstResidue(),pIt->getNumResidues(),pIt->phis,pIt->psis,direction);
            
            /* Negate angle deltas again: */
            for(int j=0;j<pIt->getNumResidues();++j)
            {
                pIt->phis[j]=-pIt->phis[j];
                pIt->psis[j]=-pIt->psis[j];
            }
        }
    }
    
    /* Apply protein transformation: */
    protein->transform(Geometry::AffineTransformation<MD::Scalar,3>(proteinTransformation));
}

/***************************
Methods of class UndoBuffer:
***************************/

UndoBuffer::UndoBuffer(void)
    :queueTail(queue.end()),
     currentProtein(0),lastSavedDepth(0)
{
}

UndoBuffer::~UndoBuffer(void)
{
}

void UndoBuffer::startInteraction(MD::Protein* protein)
{
    /* Retrieve parameters from given structure: */
    currentProtein=protein;
}

void UndoBuffer::addResidueSequence(int startResidue,int numResidues,int direction)
{
    /* Create a single residue sequence: */
    int i=direction==-1?0:1;
    DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
    prev.resize(startResidue,numResidues);
    
    /* Retrieve current dihedral angles: */
    currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.phis,prev.psis);
}

void UndoBuffer::startInteraction(MD::Protein* protein,MD::Protein::Residue* selectedResidue,int direction)
{
    /* Retrieve parameters from given structure: */
    currentProtein=protein;
    
    /* Create a single residue sequence: */
    int i=direction==-1?0:1;
    DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
    prev.resize(protein->getResidueIndex(selectedResidue),1);
    
    /* Retrieve current dihedral angles: */
    currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.phis,prev.psis);
}

void UndoBuffer::startInteraction(MD::Protein* protein,int startResidue,int numResidues,int direction)
{
    /* Retrieve parameters from given structure: */
    currentProtein=protein;
    
    /* Create a single residue sequence: */
    int i=direction==-1?0:1;
    DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
    prev.resize(startResidue,numResidues);
    
    /* Retrieve current dihedral angles: */
    currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.phis,prev.psis);
}

void UndoBuffer::startInteraction(const MD::Protein::StructureSelector& selectedStructure,int direction)
{
    /* Retrieve parameters from given structure: */
    currentProtein=selectedStructure.getProtein();
    
    /* Create a single residue sequence: */
    int i=direction==-1?0:1;
    DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
    prev.resize(selectedStructure.getFirstResidueIndex(),selectedStructure.getNumResidues());
    
    /* Retrieve current dihedral angles: */
    currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.phis,prev.psis);
}

void UndoBuffer::startInteraction(const std::list<MD::Protein::StructureSelector>& leftStructures,const std::list<MD::Protein::StructureSelector>& rightStructures)
{
    /* Retrieve parameters from given structure: */
    currentProtein=0;
    
    /* Create left/right residue sequences: */
    for(int i=0;i<2;++i)
    {
        const std::list<MD::Protein::StructureSelector>& str=i==0?leftStructures:rightStructures;
        for(std::list<MD::Protein::StructureSelector>::const_iterator strIt=str.begin();strIt!=str.end();++strIt)
        {
            /* Get protein from structure (silently assume all structures belong to the same protein): */
            currentProtein=strIt->getProtein();
            
            /* Create a single residue sequence: */
            DihedralAngleSequence& prev=*previousSequences[i].insert(previousSequences[i].end(),DihedralAngleSequence());
            prev.resize(strIt->getFirstResidueIndex(),strIt->getNumResidues());

            /* Retrieve current dihedral angles: */
            currentProtein->getDihedralAngles(prev.getFirstResidue(),prev.getNumResidues(),prev.phis,prev.psis);
        }
    }
}
void UndoBuffer::clear(void)
{
	if(queue.size()>0)
		queue.clear();
    lastSavedDepth=0;
	queueTail=queue.end();
    /* Clean up: */
    currentProtein=0;
    if(previousSequences[0].size()>0)
		previousSequences[0].clear();
    if(previousSequences[1].size()>0)
    	previousSequences[1].clear();
}

void UndoBuffer::finishInteraction(void)
{
    if(currentProtein!=0)
    {
        /* Insert a new undo queue entry before the current tail: */
        std::list<UndoQueueEntry>::iterator newIt=queue.insert(queueTail,UndoQueueEntry());

		/* Remove any queue entries behind the current tail (therefore flushing the redo buffer): */
        queue.erase(queueTail,queue.end());
		queueTail=queue.end();

        if(lastSavedDepth<0) // The last saved state has been flushed from the undo buffer
        {
            /* Invalidate the marker: */
            lastSavedDepth=queue.size(); // One-behind-last will never reach zero
        }
        else
            ++lastSavedDepth;

        /* Set the undo queue entry's data: */
        newIt->set(currentProtein,Geometry::OrthonormalTransformation<double,3>::identity,previousSequences);
    }
//	else
//		queue.clear();
    
    /* Clean up: */
    currentProtein=0;
    previousSequences[0].clear();
    previousSequences[1].clear();
}

void UndoBuffer::finishInteraction(const Geometry::OrthonormalTransformation<double,3>& proteinTransformation)
{
    if(currentProtein!=0)
    {
        /* Insert a new undo queue entry before the current tail: */
        std::list<UndoQueueEntry>::iterator newIt=queue.insert(queueTail,UndoQueueEntry());

        /* Remove any queue entries behind the current tail (therefore flushing the redo buffer): */
        queue.erase(queueTail,queue.end());
        queueTail=queue.end();

        if(lastSavedDepth<0) // The last saved state has been flushed from the undo buffer
        {
            /* Invalidate the marker: */
            lastSavedDepth=queue.size(); // One-behind-last will never reach zero
        }
        else
            ++lastSavedDepth;

        /* Set the undo queue entry's data: */
        newIt->set(currentProtein,proteinTransformation,previousSequences);
    }
//	else
//		queue.clear();
    
    /* Clean up: */
    currentProtein=0;
    previousSequences[0].clear();
    previousSequences[1].clear();
}

void UndoBuffer::undo(void)
{
    /* Rewind the queue tail by one and un-apply the skipped entry's dihedral angle change: */
    if(queueTail!=queue.begin())
    {
        --queueTail;
        queueTail->undo();
        --lastSavedDepth;
    }
}

void UndoBuffer::redo(void)
{
    /* Advance the queue tail by one and re-apply the skipped entry's dihedral angle change: */
    if(queueTail!=queue.end())
    {
        queueTail->redo();
        ++queueTail;
        ++lastSavedDepth;
    }
}

void UndoBuffer::deleteProtein(const MD::Protein* protein)
{
    std::list<UndoQueueEntry>::iterator qIt=queue.begin();
    while(qIt!=queue.end())
    {
        if(qIt->getProtein()==protein)
        {
            /* Delete this entry: */
            std::list<UndoQueueEntry>::iterator del=qIt;
            ++qIt;
            queue.erase(del);
        }
        else
            ++qIt;
    }
}
