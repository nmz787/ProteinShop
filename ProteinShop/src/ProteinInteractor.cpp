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

#include <vector>
#include <errno.h>
#include <stdio.h>
#include <Math/Math.h>
#include <Geometry/AffineTransformation.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLTransformations.h>
#include <GLContextData.h>
#include <ProteinInteractor.h>

// the following semi-global variables are interdependency hacks
extern char inputfilename[];
extern double lastIKResidual;

/* External global GUI functions: */
extern void updateProteinNow(void);

/**********************************
Methods of class ProteinInteractor:
**********************************/

ProteinInteractor::ProteinInteractor(MD::Protein* sProtein,MD::ProteinRenderer* sRenderer,EnergyCalculator* sEnergyCalculator,UndoBuffer* sUndoBuffer)
    :protein(sProtein),renderer(sRenderer),energyCalculator(sEnergyCalculator),undoBuffer(sUndoBuffer),
     updateDirection(LEFT),
     selectedResidue(0),
     ik(0),jointMap(0),
     currentPhis(0),currentPsis(0),deltaPhis(0),deltaPsis(0)
{
}

ProteinInteractor::~ProteinInteractor(void)
{
    delete ik;
    delete[] jointMap;
    delete[] currentPhis;
    delete[] currentPsis;
    delete[] deltaPhis;
    delete[] deltaPsis;
}

void ProteinInteractor::resetDragBox(void)
{
    if(selectedStructure.isValid())
    {
        /* Create a drag box: */
        box=selectedStructure.createDragBox();
        
        /* Set the drag box's center of rotation: */
        if(selectedStructure.containsResidue(selectedResidue)&&selectedResidue->getType()!=MD::Protein::Residue::UNK)
        {
            if(selectedResidue->getType()==MD::Protein::Residue::PRO)
            {
                /* Set the drag box's center of rotation to the selected residue's hydrogen bond site: */
                MD::Protein::Dipole co=selectedResidue->getCarboxyl();
                box.setRotateCenter(co.getBondSite());
            }
            else
            {
                /* Set the drag box's center of rotation to the midpoint of the selected residue's hydrogen bond sites: */
                MD::Protein::Dipole nh=selectedResidue->getAmide();
                MD::Protein::Dipole co=selectedResidue->getCarboxyl();
                box.setRotateCenter(Geometry::mid(nh.getBondSite(),co.getBondSite()));
            }
        }
        else if(updateDirection==LEFT)
        {
            /* Set the drag box's center of rotation to the first structure atom: */
            MD::Protein::BackboneIterator bbIt=selectedStructure.beginBackbone();
            box.setRotateCenter(bbIt.getAtom1()->getPosition());
        }
        else
        {
            /* Set the drag box's center of rotation to the last structure atom: */
            MD::Protein::BackboneIterator bbIt=selectedStructure.endBackbone();
            box.setRotateCenter(bbIt.getAtom2()->getPosition());
        }
        
        /* Set drag box' picking mode: */
        box.setPickMode(DragBox::OPAQUE);
    }
}

void ProteinInteractor::clearCoilRegions(void)
{
    /* Reset the active coil regions: */
    activeCoils.clear();
    renderer->resetAllBackboneColor();
}

void ProteinInteractor::resetCoilRegions(void)
{
    /* Reset the active coil regions: */
    activeCoils.clear();
    renderer->resetAllBackboneColor();
    if(selectedStructure.isValid())
    {
        /* Set the default active coil region: */
        if(updateDirection==LEFT)
        {
            /* Find the next coil region to the left: */
            MD::Protein::StructureSelector coil=selectedStructure;
            while(coil.isValid()&&coil.getStructureType()!=MD::Protein::SecondaryStructure::COIL)
                coil=coil.getPred();
            if(coil.isValid())
            {
                activeCoils.push_back(coil);
                renderer->selectStructure(coil);
            }
        }
        else
        {
            /* Find the next coil region to the right: */
            MD::Protein::StructureSelector coil=selectedStructure;
            while(coil.isValid()&&coil.getStructureType()!=MD::Protein::SecondaryStructure::COIL)
                coil=coil.getSucc();
            if(coil.isValid())
            {
                activeCoils.push_back(coil);
                renderer->selectStructure(coil);
            }
        }
    }
}

void ProteinInteractor::toggleUpdateDirection(void)
{
    /* Change the update direction: */
    if(updateDirection==LEFT)
        updateDirection=RIGHT;
    else
        updateDirection=LEFT;
    
    /* Reset drag box and active coil regions: */
    resetDragBox();
    resetCoilRegions();
}

void ProteinInteractor::selectStructure(const MD::Protein::StructureSelector& newSelectedStructure)
{
    if(newSelectedStructure.getProtein()==protein)
    {
        selectedStructure=newSelectedStructure;

        /* Create a new drag box: */
        resetDragBox();
        
        /* Reset the active coil regions: */
        resetCoilRegions();
    }
}

void ProteinInteractor::setDihedralAngles(const MD::Scalar phis[],const MD::Scalar psis[])
{
    selectedStructure.setDihedralAngles(phis,psis,updateDirection==LEFT?1:-1);
    resetDragBox();
}

void ProteinInteractor::changeDihedralAngles(const MD::Scalar deltaPhis[],const MD::Scalar deltaPsis[])
{
    selectedStructure.changeDihedralAngles(deltaPhis,deltaPsis,updateDirection==LEFT?1:-1);
    resetDragBox();
}

void ProteinInteractor::toggleCoil(const MD::Protein::StructureSelector& coil)
{
    if(selectedStructure.isValid()&&coil.isValid()&&coil.getProtein()==protein&&coil.getStructureType()==MD::Protein::SecondaryStructure::COIL)
    {
        /* Find coil region in active coil list: */
        std::list<MD::Protein::StructureSelector>::iterator it;
        for(it=activeCoils.begin();it!=activeCoils.end()&&it->getStructureIndex()<coil.getStructureIndex();++it)
            ;
        if(it!=activeCoils.end()&&it->getStructureIndex()==coil.getStructureIndex())
        {
            /* Remove coil from list: */
            activeCoils.erase(it);
            renderer->deselectStructure(coil);
        }
        else
        {
            /* Insert coil into list: */
            activeCoils.insert(it,coil);
            renderer->selectStructure(coil);
        }
    }
}

bool ProteinInteractor::startInteraction(void)
{
    bool ok=true;

    /* Count number of rotatable backbone bonds in the update region: */
    numJoints=0;
    for(StructureList::iterator acIt=activeCoils.begin();acIt!=activeCoils.end();++acIt)
        for(MD::Protein::BackboneIterator bbIt=acIt->beginBackbone();bbIt!=acIt->endBackbone();++bbIt)
            if(bbIt.isRotatable())
                ++numJoints;
	// no structure, no joint
	if(!numJoints)
		return false;
    /* Create IK object: */
    ik=new IK(numJoints);

    /* Create translation structures from IK object to protein: */
    jointMap=new JointMapItem[numJoints];
    firstResidueIndex=activeCoils.front().getFirstResidueIndex();
    numResidues=activeCoils.back().getFirstResidueIndex()+activeCoils.back().getNumResidues()-firstResidueIndex;
    currentPhis=new MD::Scalar[numResidues];
    currentPsis=new MD::Scalar[numResidues];
    deltaPhis=new MD::Scalar[numResidues];
    deltaPsis=new MD::Scalar[numResidues];
    
    /* Determine the IK update direction: */
    if(activeCoils.back().getStructureIndex()<selectedStructure.getStructureIndex())
    {
        /* Unidirectional IK towards the C=O-terminus (to the right): */
        ikUpdateDirection=1;
    }
    else if(activeCoils.front().getStructureIndex()>selectedStructure.getStructureIndex())
    {
        /* Unidirectional IK towards the N-H-terminus (to the left): */
        ikUpdateDirection=-1;
    }
    else
    {
        /* Bidirectional IK: */
        ikUpdateDirection=updateDirection==LEFT?1:-1;
    }
    
    /* Initialize the joint chain: */
    if(ikUpdateDirection==1)
    {
        effectorJointIndex=numJoints;
        StructureList::iterator coilIt=activeCoils.begin();
        MD::Protein::BackboneIterator bbIt;
        int jointIndex=0;
        MD::Vector currentOffset=MD::Vector::zero;
        while(coilIt!=activeCoils.end())
        {
            if(effectorJointIndex==numJoints&&coilIt->getStructureIndex()>selectedStructure.getStructureIndex())
                effectorJointIndex=jointIndex;
            
            /* Add joints for all residues in this coil region: */
            bbIt=coilIt->beginBackbone();
            while(bbIt!=coilIt->endBackbone())
            {
                /* Add one rotation joint for each rotatable bond in the coil: */
                if(bbIt.isRotatable())
                {
                    MD::Protein::BondRotator rotator(protein,bbIt);
                    MD::Point rotPoint=rotator.getStart();
                    MD::Point offset=rotPoint-currentOffset;
                    currentOffset=rotPoint-MD::Point::origin;
                    ik->setJoint(jointIndex,offset,rotator.getAxis());
                    ik->setJointAngle(jointIndex,IK::Scalar(0));
                    jointMap[jointIndex]=std::make_pair(bbIt.getResidueIndex()-firstResidueIndex,bbIt.getBondType());
                    ++jointIndex;
                }
                
                ++bbIt;
            }
            
            ++coilIt;
        }
        
        /* Set leaf offset: */
        MD::Point offset=bbIt.getAtom1()->getPosition()-currentOffset;
        ik->setLeafOffset(offset);
    }
    else
    {
        effectorJointIndex=numJoints;
        StructureList::iterator coilIt=activeCoils.end();
        MD::Protein::BackboneIterator bbIt;
        int jointIndex=0;
        MD::Vector currentOffset=MD::Vector::zero;
        while(coilIt!=activeCoils.begin())
        {
            --coilIt;
            
            if(effectorJointIndex==numJoints&&coilIt->getStructureIndex()<selectedStructure.getStructureIndex())
                effectorJointIndex=jointIndex;
            
            /* Add joints for all residues in this coil region: */
            bbIt=coilIt->endBackbone();
            while(bbIt!=coilIt->beginBackbone())
            {
                --bbIt;
                
                /* Add one rotation joint for each rotatable bond in the coil: */
                if(bbIt.isRotatable())
                {
                    MD::Protein::BondRotator rotator(protein,bbIt);
                    MD::Point rotPoint=rotator.getEnd();
                    MD::Point offset=rotPoint-currentOffset;
                    currentOffset=rotPoint-MD::Point::origin;
                    ik->setJoint(jointIndex,offset,-rotator.getAxis());
                    ik->setJointAngle(jointIndex,IK::Scalar(0));
                    jointMap[jointIndex]=std::make_pair(bbIt.getResidueIndex()-firstResidueIndex,bbIt.getBondType());
                    ++jointIndex;
                }
            }
        }
        
        /* Set leaf offset: */
        MD::Point offset=bbIt.getAtom2()->getPosition()-currentOffset;
        ik->setLeafOffset(offset);
    }
    
    /* Initialize IK parameters: */
    IK::Scalar orientationWeights[3]={10,10,10};
    ik->setOrientationWeights(orientationWeights);
    ik->setMaxResidual(IK::Scalar(1.0e-12));
    ik->setMaxStepError(IK::Scalar(1.0e-3));
    effectorStepSize=IK::Scalar(1.0e-4);
    correctionStepSize=IK::Scalar(1.0e-4);

    /* Initialize angle arrays: */
    for(int i=0;i<numResidues;++i)
    {
        currentPhis[i]=MD::Scalar(0);
        currentPsis[i]=MD::Scalar(0);
        deltaPhis[i]=MD::Scalar(0);
        deltaPsis[i]=MD::Scalar(0);
    }

    /* Calculate the initial effector transformation: */
    ik->setEffectorJointIndex(effectorJointIndex);
    initialTransformation=ik->getEffectorTransformation();
    ik->setEffectorJointIndex(numJoints);
    correctionTransformation=ik->getEffectorTransformation();

    /* Open undo buffer: */
    undoBuffer->startInteraction(protein,firstResidueIndex,numResidues,ikUpdateDirection);

    return ok;
}

bool ProteinInteractor::pickDragBox(const DragBox::HTransformation& modelView,const DragBox::HTransformation& projection,const DragBox::Point& mouseClip)
{
    /* Try picking the interaction box: */
    if(!selectedStructure.isValid()||!box.pick(modelView,projection,mouseClip))
        return false;

    bool ok=startInteraction();
    if(!ok)
        box.release();
    return ok;
}

bool ProteinInteractor::pickDragBox(const DragBox::Transformation& initialDragTransformation)
{
    /* Try picking the interaction box: */
    if(!selectedStructure.isValid()||!box.pick(initialDragTransformation))
        return false;

    bool ok=startInteraction();
    if(!ok)
        box.release();
    return ok;
}

void ProteinInteractor::drag(const ProteinInteractor::Transformation& goalTransformation)
{
    /* Calculate the difference transformation: */
    Transformation ikTransformation=goalTransformation*initialTransformation;

    /*********************************************************************
    Perform IK steps on the joint chain. First, move the dragged structure
    towards its goal position. In the case of bi-directional IK, this will
    incur (small) drift in the leaf joint in the joint chain. Therefore,
    only in the case of bi-directional IK, the leaf joint is moved back
    towards its original position in a second step.
    There are two step sizes involved here: The first is the step size for
    moving the dragged structure, the second is the step size for moving
    the leaf joint. Since these IK operations work on different scales, we
    maintain the two step sizes independently.
    *********************************************************************/

    /* Move the dragged structure towards its goal position: */
    ik->setEffectorJointIndex(effectorJointIndex); // Select the dragged structure for IK
    ik->setStepSize(effectorStepSize); // Set the last used IK step size
    IK::Scalar finalResidual=ik->performSteps(ikTransformation,2000); // Perform the IK steps and save the residual
    effectorStepSize=ik->getStepSize(); // Save the IK step size

    /* Check for bi-directional IK: */
    if(effectorJointIndex<numJoints)
    {
        /* Move the chain's leaf joint back towards its original position: */
        ik->setEffectorJointIndex(numJoints); // Select the chain's leaf joint
        ik->setStepSize(correctionStepSize); // Set the last used IK step size
        ik->performSteps(correctionTransformation,2000); // Perform the IK steps
        correctionStepSize=ik->getStepSize(); // Save the IK step size
    }

    /* Update the coils' dihedral angles according to the joint tree: */
    for(int i=0;i<numJoints;++i)
    {
        IK::Scalar newAngle=ik->getJointAngle(i);
        int angleIndex=jointMap[i].first;
        if(jointMap[i].second==MD::Protein::BackboneIterator::PHI)
        {
            deltaPhis[angleIndex]+=newAngle-currentPhis[angleIndex];
            currentPhis[angleIndex]=newAngle;
        }
        else
        {
            deltaPsis[angleIndex]+=newAngle-currentPsis[angleIndex];
            currentPsis[angleIndex]=newAngle;
        }
    }

    /* Get the last residual for BuildBeta: */
    lastIKResidual = finalResidual;
}

void ProteinInteractor::applyChanges(void)
{
    /* Apply dihedral angle changes to the protein: */
    protein->changeDihedralAngles(firstResidueIndex,numResidues,deltaPhis,deltaPsis,ikUpdateDirection,effectorJointIndex<numJoints);

    /* Reset dihedral angle differences: */
    for(int i=0;i<numResidues;++i)
    {
        deltaPhis[i]=IK::Scalar(0);
        deltaPsis[i]=IK::Scalar(0);
    }

    /* Update related data structures: */
    renderer->updateProtein();
    if(energyCalculator!=0)
        energyCalculator->updateProtein();
}

void ProteinInteractor::finishInteraction(void)
{
    /* Close undo buffer: */
    undoBuffer->finishInteraction();

    /* Delete IK data structures: */
    delete ik;
    ik=0;
    delete[] jointMap;
    jointMap=0;
    delete[] currentPhis;
    currentPhis=0;
    delete[] currentPsis;
    currentPsis=0;
    delete[] deltaPhis;
    deltaPhis=0;
    delete[] deltaPsis;
    deltaPsis=0;
    
    /* Snap drag box back to structure: */
    resetDragBox();
}

void ProteinInteractor::glRenderAction(GLContextData& contextData) const
{
    if(selectedStructure.isValid())
    {
        glColor4f(0.0f,0.5f,0.0f,0.33f);
        glLineWidth(3.0);
        glPointSize(5.0);
        glDisable(GL_CULL_FACE);
        box.draw();
        glEnable(GL_CULL_FACE);
    }
}

void ProteinInteractor::glTransformProtein(GLContextData& contextData) const
{
}

bool ProteinInteractor::readCoilRegionAngles(const char* angleFilename)
{
    /* Open input file: */
    FILE* angleFile=fopen(angleFilename,"rt");
    if(angleFile==0)
        return false;
    strcpy(inputfilename, angleFilename);

    /* Save previous protein state: */
    undoBuffer->startInteraction(protein);
    
    /* Read all coil regions from file: */
    int numAnimationSteps=100;
    std::list<CoilRegionAngles*> coilRegionAngleDeltas;
    while(true)
    {
        /* Read coil description from file: */
        int structureIndex,firstResidueIndex,numResidues;
        if(fscanf(angleFile,"%d %d %d",&structureIndex,&firstResidueIndex,&numResidues)!=3)
            break;
        
        /* Get appropriate structure selector and check consistency: */
        MD::Protein::StructureSelector structure=protein->pickStructure(structureIndex);
        if(structure.getFirstResidueIndex()==firstResidueIndex&&structure.getNumResidues()==numResidues)
        {
            /* Save current angles in undo buffer: */
            undoBuffer->addResidueSequence(firstResidueIndex,numResidues,updateDirection==LEFT?-1:1);
            
            /* Retrieve current angles from structure: */
            CoilRegionAngles currentAngles(structure);
            currentAngles.getCurrentAngles();
            
            /* Read target dihedral angles: */
            CoilRegionAngles targetAngles(structure);
            for(int i=0;i<numResidues;++i)
            {
                double phi,psi;
                fscanf(angleFile,"%lf %lf",&phi,&psi);
                targetAngles.phis[i]=Math::rad(MD::Scalar(phi));
                targetAngles.psis[i]=Math::rad(MD::Scalar(psi));
            }
            
            /* Calculate and save angle deltas: */
            CoilRegionAngles* angleDeltas=new CoilRegionAngles(structure);
            for(int i=0;i<angleDeltas->getNumResidues();++i)
            {
                double dPhi=(targetAngles.phis[i]-currentAngles.phis[i]);
                if(dPhi<=-Math::Constants<double>::pi)
                    dPhi+=2.0*Math::Constants<double>::pi;
                else if(dPhi>Math::Constants<double>::pi)
                    dPhi-=2.0*Math::Constants<double>::pi;
                angleDeltas->phis[i]=MD::Scalar(dPhi/double(numAnimationSteps));
                double dPsi=(targetAngles.psis[i]-currentAngles.psis[i]);
                if(dPsi<=-Math::Constants<double>::pi)
                    dPsi+=2.0*Math::Constants<double>::pi;
                else if(dPsi>Math::Constants<double>::pi)
                    dPsi-=2.0*Math::Constants<double>::pi;
                angleDeltas->psis[i]=MD::Scalar(dPsi/double(numAnimationSteps));
            }
            coilRegionAngleDeltas.push_back(angleDeltas);
        }
    }
    
    fclose(angleFile);
    
    /* Morph protein from current coil region angles to target angles: */
    for(int i=0;i<numAnimationSteps;++i)
    {
        for(std::list<CoilRegionAngles*>::iterator craIt=coilRegionAngleDeltas.begin();craIt!=coilRegionAngleDeltas.end();++craIt)
            (*craIt)->changeDihedralAngles(updateDirection);
        updateProteinNow();
    }
    if(energyCalculator!=0)
        energyCalculator->updateProtein();
    
    /* Delete angle deltas: */
    for(std::list<CoilRegionAngles*>::iterator craIt=coilRegionAngleDeltas.begin();craIt!=coilRegionAngleDeltas.end();++craIt)
        delete *craIt;
    coilRegionAngleDeltas.clear();
    
    /* Finish interaction: */
    undoBuffer->finishInteraction();
    
    /* Snap drag box back to structure: */
    resetDragBox();
    
    return true;
}

bool ProteinInteractor::writeCoilRegionAngles(const char* angleFilename) const
{
    /* Open output file: */
    FILE* angleFile=fopen(angleFilename,"wt");
    if(angleFile==0)
        return false;
    
    /* Iterate through all active coils: */
    for(std::list<MD::Protein::StructureSelector>::const_iterator cIt=activeCoils.begin();cIt!=activeCoils.end();++cIt)
    {
        /* Write coil region header: */
        fprintf(angleFile,"%d %d %d\n",cIt->getStructureIndex(),cIt->getFirstResidueIndex(),cIt->getNumResidues());
        
        /* Retrieve dihedral angles: */
        MD::Scalar* phis=new MD::Scalar[cIt->getNumResidues()];
        MD::Scalar* psis=new MD::Scalar[cIt->getNumResidues()];
        cIt->getDihedralAngles(phis,psis);
        
        /* Write dihedral angles: */
        for(int i=0;i<cIt->getNumResidues();++i)
            fprintf(angleFile,"%10.4f %10.4f\n",double(Math::deg(phis[i])),double(Math::deg(psis[i])));
        delete[] phis;
        delete[] psis;
    }
    
    fclose(angleFile);
    return true;
}
