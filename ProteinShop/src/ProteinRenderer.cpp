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
ProteinRenderer - Class to visualize protein structures.

***********************************************************************/

#include <Math/Constants.h>
#include <GLGeometry.h>

#include "ProteinRenderer.h"
#include "ProteinGui.h"

#include "GLSphere.cpp"

namespace MD {
extern GLText *text;

/************************************************
Methods of class ProteinRenderer::StructureFlags:
************************************************/

ProteinRenderer::StructureFlags::StructureFlags(const ConfigurationFile::SectionIterator& configFileSection,const Protein::SecondaryStructure* sStructure)
    :structure(sStructure),
     drawAtoms(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawAtoms",true)),
     drawBonds(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawBonds",true)),
     drawCPK(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawCPK",true)),
     drawTube(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawTube",true)),
     drawLine(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawLine",true)),
     bondWidth(configFileSection.retrieveValue<float,ValueCoder<float> >("./bondWidth",2.0f)),
     bondColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./bondColor",Color(1.0,0.0,0.0))),
     drawBackbone(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawBackbone",true)),
     backboneWidth(configFileSection.retrieveValue<float,ValueCoder<float> >("./backboneWidth",4.0f)),
     drawBackboneRibbon(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawBackboneRibbon",true)),
     backboneRibbonWidth(configFileSection.retrieveValue<Scalar,ValueCoder<Scalar> >("./backboneRibbonWidth",Scalar(2))),
     drawCartoon(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawCartoon",true)),
     cartoonSpline1(0),cartoonSpline2(0),cartoonP(0),cartoonX(0),cartoonY(0),cartoonZ(0),
     drawHydrogenBondSites(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawHydrogenBondSites",structure->getStructureType()==Protein::SecondaryStructure::BETA_STRAND)),
     hydrogenBondSiteDiameter(configFileSection.retrieveValue<float,ValueCoder<float> >("./hydrogenBondSiteDiameter",5.0f)),
     hydrogenBondSiteWidth(configFileSection.retrieveValue<float,ValueCoder<float> >("./hydrogenBondSiteWidth",1.0f)),
     amideColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./amideColor",Color(0.5,0.5,1.0))),
     carboxylColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./carboxylColor",Color(1.0,1.0,0.0))),
     drawHydrogenCages(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawHydrogenCages",structure->getStructureType()==Protein::SecondaryStructure::BETA_STRAND)),
     hydrogenCageWidth(configFileSection.retrieveValue<float,ValueCoder<float> >("./hydrogenCageWidth",1.5f)),
     hydrogenCageColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./hydrogenCageColor",Color(0.6,0.5,0.0))),
     hydrogenCageLarge(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./hydrogenCageLarge",structure->getStructureType()==Protein::SecondaryStructure::BETA_STRAND))
{
}

ProteinRenderer::StructureFlags::~StructureFlags(void)
{
    delete cartoonSpline1;
    delete cartoonSpline2;
    delete[] cartoonP;
    delete[] cartoonX;
    delete[] cartoonY;
    delete[] cartoonZ;
}

/****************************************
Static elements of class ProteinRenderer:
****************************************/

const double ProteinRenderer::maxAtomRadius=1.8;
const ProteinRenderer::Color ProteinRenderer::elementColors[118]={
    Color(0.8f,0.8f,0.8f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.2f,0.2f,0.2f),Color(0.5f,0.5f,1.0f),Color(1.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,1.0f,0.0f),Color(1.0f,1.0f,0.0f),
    Color(0.0f,1.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(1.0f,0.5f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f),
    Color(0.0f,0.0f,0.0f),Color(0.0f,0.0f,0.0f)
};

/********************************
Methods of class ProteinRenderer:
********************************/

void ProteinRenderer::createBackboneRibbonSpline(void)
{
    /* Count the total number of backbone atoms: */
    int numBackboneAtoms=0;
    for(const Protein::Residue* rPtr=protein->residues;rPtr!=0;rPtr=rPtr->succ)
    {
        if(backboneRibbonUseAllAtoms)
            numBackboneAtoms+=rPtr->backbone.size();
        else
            numBackboneAtoms+=1;
    }
    
    /* Create knot values for the backbone ribbon spline: */
    int numKnots=numBackboneAtoms+backboneRibbonDegree-1;
    Scalar* knots=new Scalar[numKnots];
    for(int i=0;i<backboneRibbonDegree;++i)
        knots[i]=Scalar(0);
    for(int i=backboneRibbonDegree;i<numBackboneAtoms-1;++i)
        knots[i]=Scalar(i-(backboneRibbonDegree-1));
    for(int i=numBackboneAtoms-1;i<numKnots;++i)
        knots[i]=Scalar(numBackboneAtoms-backboneRibbonDegree);
    
    /* Create control points for the backbone ribbon spline: */
    Point* points=new Point[numBackboneAtoms];
    Point* pPtr=points;
    int firstControlPointIndex=0;
    Scalar boundaryParameterValue=Scalar(0);
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        const Protein::SecondaryStructure* sPtr=sfIt->structure;
        
        /* Iterate through all residues in this structure: */
        int numControlPoints=0;
        for(const Protein::Residue* rPtr=sPtr->residueBegin;rPtr!=sPtr->residueEnd;rPtr=rPtr->succ)
        {
            /* Create control point(s) for this residue: */
            if(backboneRibbonUseAllAtoms)
            {
                for(std::vector<Protein::ChainAtom*>::const_iterator aIt=rPtr->backbone.begin();aIt!=rPtr->backbone.end();++aIt,++pPtr)
                    *pPtr=(*aIt)->getPosition();
                numControlPoints+=rPtr->backbone.size();
            }
            else
            {
                *pPtr=rPtr->backbone[1]->getPosition();
                ++pPtr;
                numControlPoints+=1;
            }
        }
        
        /* Calculate spline parameters for this structure's boundaries: */
        sfIt->backboneRibbonNumSamples=numControlPoints*backboneRibbonSampleDensity;
        sfIt->backboneRibbonParameterMin=boundaryParameterValue;
        firstControlPointIndex+=numControlPoints;
        if(firstControlPointIndex<numBackboneAtoms)
        {
            boundaryParameterValue=knots[firstControlPointIndex-1]+knots[firstControlPointIndex+backboneRibbonDegree-1];
            for(int i=0;i<backboneRibbonDegree-1;++i)
                boundaryParameterValue+=knots[firstControlPointIndex+i]*Scalar(2);
            boundaryParameterValue/=Scalar(backboneRibbonDegree*2);
        }
        else
            boundaryParameterValue=Scalar(numBackboneAtoms-backboneRibbonDegree);
        sfIt->backboneRibbonParameterMax=boundaryParameterValue;
    }
    
    /* Create a spline object containing all backbone atoms: */
    delete backboneRibbonSpline;
    backboneRibbonSpline=new Geometry::SplineCurve<Scalar,3>(backboneRibbonDegree,numBackboneAtoms,knots,points);
    delete[] knots;
    delete[] points;
}

void ProteinRenderer::createCartoonSpline(ProteinRenderer::StructureFlags& sf)
{
    /* Create knot vector for the cartoon rendering splines: */
    int numPoints=sf.structure->getNumResidues()+2;
    int degree=cartoonDegree;
    if(degree>sf.structure->getNumResidues()+1)
        degree=sf.structure->getNumResidues()+1;
    int numKnots=numPoints+(degree-1);
    Scalar* knots=new Scalar[numKnots];
    for(int i=0;i<degree;++i)
        knots[i]=Scalar(0);
    for(int i=degree;i<numPoints-1;++i)
        knots[i]=Scalar(i-(degree-1));
    for(int i=numPoints-1;i<numKnots;++i)
        knots[i]=Scalar(numPoints-degree);

    /* Create control points for the major cartoon rendering spline: */
    Point* points=new Point[numPoints];
    delete sf.cartoonSpline1;
    sf.cartoonSpline1=new Geometry::SplineCurve<Scalar,3>(degree,numPoints,knots,points);
    
    /* Create control points for the minor beta strand cartoon rendering spline: */
    delete sf.cartoonSpline2;
    if(sf.structure->getStructureType()==Protein::SecondaryStructure::BETA_STRAND)
        sf.cartoonSpline2=new Geometry::SplineCurve<Scalar,3>(degree,numPoints,knots,points);
    else
        sf.cartoonSpline2=0;
    
    delete[] knots;
    delete[] points;
    
    /* Create evaluation arrays for the cartoon rendering splines: */
    sf.numCartoonSplineSamples=sf.cartoonSpline1->getNumSegments()*cartoonSampleDensity;
    sf.maxCartoonSplineParameter=sf.cartoonSpline1->getUMax();
    delete sf.cartoonP;
    sf.cartoonP=new Point[sf.numCartoonSplineSamples+1];
    delete sf.cartoonX;
    sf.cartoonX=new Vector[sf.numCartoonSplineSamples+1];
    delete sf.cartoonY;
    sf.cartoonY=new Vector[sf.numCartoonSplineSamples+1];
    delete sf.cartoonZ;
    sf.cartoonZ=new Vector[sf.numCartoonSplineSamples+1];
}

void ProteinRenderer::updateCartoonSpline(ProteinRenderer::StructureFlags& sf)
{
    if(sf.drawCartoon || sf.drawTube)
    {
        /* Update major cartoon spline's control points: */
        const Protein::Residue* rPtr=sf.structure->residueBegin;
        Point p=rPtr->backbone[0]->getPosition();
        if(rPtr->pred!=0)
            p=Geometry::mid(p,rPtr->pred->backbone[rPtr->pred->backbone.size()-1]->getPosition());
        sf.cartoonSpline1->setPoint(0,p);
        int index;
        for(rPtr=sf.structure->residueBegin,index=1;rPtr!=sf.structure->residueEnd;rPtr=rPtr->succ,++index)
        {
            if(rPtr->getType()!=Protein::Residue::UNK || rPtr->getType()!=Protein::Residue::NME)
            {
                // Program received signal SIGSEGV, Segmentation fault.
                // (gdb) print rPtr->backbone[2]
                // $4 = (MD::Protein::ChainAtom * const&) @0x81d1a40: 0x0
                // original line of code is redacted and substituted:
                //p=Geometry::mid(rPtr->backbone[0]->getPosition(),rPtr->backbone[2]->getPosition());

				// check the size to avoid the ghost backbone[2] in NME
                if ( rPtr->backbone.size() > 2 )
                    p= Geometry::mid (
                        rPtr->backbone[0]->getPosition(),
                        rPtr->backbone[2]->getPosition()
                    );
                else
                    p = rPtr->backbone[0]->getPosition();
                p=Geometry::mid(p,rPtr->backbone[1]->getPosition());
            }
            else if (rPtr->getType()!=Protein::Residue::NME)
                p=Geometry::mid(rPtr->backbone[0]->getPosition(),rPtr->backbone[rPtr->backbone.size()-1]->getPosition());
				else
				p=rPtr->backbone[rPtr->backbone.size()-1]->getPosition();

            sf.cartoonSpline1->setPoint(index,p);
        }
        for(rPtr=sf.structure->residueBegin;rPtr->succ!=sf.structure->residueEnd;rPtr=rPtr->succ)
            ;
        p=rPtr->backbone[rPtr->backbone.size()-1]->getPosition();
        if(rPtr->succ!=0)
            p=Geometry::mid(p,rPtr->succ->backbone[0]->getPosition());
        sf.cartoonSpline1->setPoint(index,p);

        /* Update minor beta strand cartoon spline's control points: */
        if(sf.structure->getStructureType()==Protein::SecondaryStructure::BETA_STRAND)
        {
            Scalar crossFactor=Scalar(1);
            for(rPtr=sf.structure->residueBegin,index=1;rPtr!=sf.structure->residueEnd;rPtr=rPtr->succ,++index)
            {
                std::vector<Atom*>::const_iterator aIt;
                for(aIt=rPtr->backbone[0]->getBonds().begin();aIt!=rPtr->backbone[0]->getBonds().end();++aIt)
                    if(static_cast<const Protein::ChainAtom*>(*aIt)->backboneIndex<0)
                        break;
                Vector d1=(*aIt)->getPosition()-rPtr->backbone[0]->getPosition();
                
				Vector d2;
				if(rPtr->backbone.size() > 2)
				{
				for(aIt=rPtr->backbone[2]->getBonds().begin();aIt!=rPtr->backbone[2]->getBonds().end();++aIt)
                    if(static_cast<const Protein::ChainAtom*>(*aIt)->backboneIndex<0)
                        break;
				d2=(*aIt)->getPosition()-rPtr->backbone[2]->getPosition();
				}
				else
				{
				for(aIt=rPtr->backbone[1]->getBonds().begin();aIt!=rPtr->backbone[1]->getBonds().end();++aIt)
                    if(static_cast<const Protein::ChainAtom*>(*aIt)->backboneIndex<0)
                        break;
                d2=(*aIt)->getPosition()-rPtr->backbone[1]->getPosition();
                }
				Vector c=d1+d2;
                // Vector c=Geometry::cross(*pPtr1-rPtr->backbone[0]->getPosition(),rPtr->backbone[2]->getPosition()-*pPtr1);
                c.normalize();
                c*=crossFactor;
                crossFactor=-crossFactor;
                sf.cartoonSpline2->setPoint(index,sf.cartoonSpline1->getPoint(index)+c);
            }
            sf.cartoonSpline2->setPoint(0,(sf.cartoonSpline2->getPoint(1)-sf.cartoonSpline1->getPoint(1))+sf.cartoonSpline1->getPoint(0));
            sf.cartoonSpline2->setPoint(index,(sf.cartoonSpline2->getPoint(index-1)-sf.cartoonSpline1->getPoint(index-1))+sf.cartoonSpline1->getPoint(index));
        }
        
        /* Evaluate cartoon splines: */
        Geometry::SplineCurve<Scalar,3>::EvaluationCache* c1=sf.cartoonSpline1->createEvaluationCache();
        switch(sf.structure->getStructureType())
        {
            case Protein::SecondaryStructure::COIL:
            {
                Vector n;
                for(int i=0;i<=sf.numCartoonSplineSamples;++i)
                {
                    /* Evaluate spline: */
                    Scalar u=Scalar(i)*sf.maxCartoonSplineParameter/Scalar(sf.numCartoonSplineSamples);
                    sf.cartoonP[i]=sf.cartoonSpline1->evaluate(u,c1,sf.cartoonY[i]);

                    /* Construct coordinate frame around cartoonP[i]: */
                    sf.cartoonY[i].normalize();
                    if(i==0)
                        n=Geometry::normal(sf.cartoonY[0]);
                    else
                        n-=(sf.cartoonY[i]*n)*sf.cartoonY[i];
                    n.normalize();
                    sf.cartoonX[i]=n;
                    sf.cartoonZ[i]=Geometry::cross(sf.cartoonX[i],sf.cartoonY[i]);
                    // sf.cartoonZ[i].normalize();
                }
                break;
            }

            case Protein::SecondaryStructure::ALPHA_HELIX:
                for(int i=0;i<=sf.numCartoonSplineSamples;++i)
                {
                    /* Evaluate spline: */
                    Scalar u=Scalar(i)*sf.maxCartoonSplineParameter/Scalar(sf.numCartoonSplineSamples);
                    sf.cartoonP[i]=sf.cartoonSpline1->evaluate(u,c1,sf.cartoonY[i],sf.cartoonZ[i]);

                    /* Construct coordinate frame around cartoonP[i]: */
                    sf.cartoonY[i].normalize();
                    sf.cartoonX[i]=Geometry::cross(sf.cartoonY[i],sf.cartoonZ[i]);
                    sf.cartoonX[i].normalize();
                    sf.cartoonZ[i]=Geometry::cross(sf.cartoonX[i],sf.cartoonY[i]);
                    // sf.cartoonZ[i].normalize();
                }
                break;

            case Protein::SecondaryStructure::BETA_STRAND:
            {
                Geometry::SplineCurve<Scalar,3>::EvaluationCache* c2=sf.cartoonSpline2->createEvaluationCache();
                for(int i=0;i<=sf.numCartoonSplineSamples;++i)
                {
                    /* Evaluate both splines: */
                    Scalar u=Scalar(i)*sf.maxCartoonSplineParameter/Scalar(sf.numCartoonSplineSamples);
                    sf.cartoonP[i]=sf.cartoonSpline1->evaluate(u,c1,sf.cartoonY[i]);
                    Point p2=sf.cartoonSpline2->evaluate(u,c2);

                    /* Construct coordinate frame around cartoonP[i]: */
                    sf.cartoonY[i].normalize();
                    sf.cartoonZ[i]=Geometry::cross(p2-sf.cartoonP[i],sf.cartoonY[i]);
                    sf.cartoonZ[i].normalize();
                    sf.cartoonX[i]=Geometry::cross(sf.cartoonY[i],sf.cartoonZ[i]);
                    // sf.cartoonX[i].normalize();
                }
                delete c2;
                break;
            }

	    default:
	        break;
        }
        delete c1;
    }
}

void ProteinRenderer::glDrawAtoms(GLContextData& contextData) const
{
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);

    /* Set atom rendering parameters: */
    glMaterial(GL_FRONT,atomMaterial);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT,GL_DIFFUSE);
    glPushMatrix();
    
    /* Render each secondary structure independently: */
    Point currentOrigin=Point::origin;
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawAtoms)
        {
            /* Iterate through all residues in this secondary structure: */
            for(const Protein::Residue* rPtr=sfIt->structure->residueBegin;rPtr!=sfIt->structure->residueEnd;rPtr=rPtr->succ)
            {
                /* Push residue (PDB) index onto the selection name stack: */
                glPushName(rPtr->residueIndex);
                
                /* Iterate through all atoms in this residue: */
                for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
                {
                    Vector displacement=aPtr->getPosition()-currentOrigin;
                    glTranslate(displacement);
                    if(mapAtomValues)
                        glColor(atomColorMap.mapColor(atomValues[aPtr->getAtomIndex()]));
                    else
                        glColor(elementColors[aPtr->getType()]);
                    glCallList(dataItem->atomSphereDisplayListBaseId+aPtr->getType());
                    currentOrigin=aPtr->getPosition();
                }
                
                glPopName();
            }
        }
        
        ++structureIndex;
        glPopName();
    }
    
    /* Reset OpenGL state: */
    glPopMatrix();
    glDisable(GL_COLOR_MATERIAL);
}

void ProteinRenderer::glDrawCPK(GLContextData& contextData) const
{
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    
    /* Set atom rendering parameters: */
	glMaterial(GL_FRONT,bondMaterial);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT,GL_DIFFUSE);

    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawCPK)
        {
       		/* Iterate through all residues in this secondary structure: */
    		for(const Protein::Residue* rPtr=sfIt->structure->residueBegin;rPtr!=sfIt->structure->residueEnd;rPtr=rPtr->succ)
    		{
    		/* Push residue (PDB) index onto the selection name stack: */
    		glPushName(rPtr->residueIndex);

			 /* Iterate through all atoms in this residue: */
        	 for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
        	 {
        		 /* Iterate through all bonds for this atom: */
            	 for(std::vector<Atom*>::const_iterator bIt=aPtr->getBonds().begin();bIt!=aPtr->getBonds().end();++bIt)
            	 {
             		 if(static_cast<const Atom*>(aPtr)<*bIt)
                	 {
                        	 /* Draw a bicolor capped cylinder around the bond axis: */
                        	 Vector bondAxis=(*bIt)->getPosition()-aPtr->getPosition();
                        	 double bondAxisLen=Geometry::mag(bondAxis);
                        	 bondAxis/=Scalar(bondAxisLen);

                        	 for(int i=0;i<=numLineVertices;++i)
                        	 {
                            	 dataItem->lineVertices[1][i][2]=GLfloat(bondAxisLen*0.5);
                        		 dataItem->lineVertices[2][i][2]=GLfloat(bondAxisLen);
                        	 }

                        	 /* Set OpenGL modelview matrix to align z axis with bond axis: */
                        	 glPushMatrix();
                        	 glTranslate(aPtr->getPosition()-Point::origin);
                        	 Vector rotAxis(-bondAxis[1],bondAxis[0],Scalar(0));
                        	 Scalar rotAngle=Math::deg(Math::acos(bondAxis[2]));
                        	 glRotate(rotAngle,rotAxis);

                        	 /* Draw root atom cap: */
                        	 glBegin(GL_POLYGON);

							 if (rPtr->getResType()==Protein::Residue::PHOBIC && showHydrophobic)
								 glColor(phobicColor);
							 else if (rPtr->getResType()==Protein::Residue::PHILIC && showHydrophilic) 
								 glColor(philicColor);
							 else if (rPtr->getResType()==Protein::Residue::DISULFIDE && showDisulfide) 
								 glColor(disulfideColor);
							 else
                           		 glColor(elementColors[aPtr->getType()]);

							 glNormal(0.0,0.0,-1.0);
							 for(int i=numLineVertices;i>=0;--i)
                        		 glVertex(dataItem->lineVertices[0][i]);
                        	 glEnd();

                    		/* Draw cylinder at root atom: */
                        	 glBegin(GL_QUAD_STRIP);
                        	 for(int i=0;i<=numLineVertices;++i)
                        	 {
                            	 glNormal(dataItem->lineNormals[i]);
                            	 glVertex(dataItem->lineVertices[1][i]);
                        		 glVertex(dataItem->lineVertices[0][i]);
                        	 }
                        	 glEnd();

							 // draw end sphere
							 glCallList(dataItem->cpkSphereDisplayListId);

                        	 /* Draw cylinder at other atom: */
                        	 glBegin(GL_QUAD_STRIP);
							 if (rPtr->getResType()==Protein::Residue::PHOBIC && showHydrophobic)
								 glColor(phobicColor);
							 else if (rPtr->getResType()==Protein::Residue::PHILIC && showHydrophilic) 
								 glColor(philicColor);
							 else if (rPtr->getResType()==Protein::Residue::DISULFIDE && showDisulfide) 
								 glColor(disulfideColor);
							 else
	                       		 glColor(elementColors[(*bIt)->getType()]);

                        	 for(int i=0;i<=numLineVertices;++i)
                        	 {
                                	 glNormal(dataItem->lineNormals[i]);
                                	 glVertex(dataItem->lineVertices[2][i]);
                                	 glVertex(dataItem->lineVertices[1][i]);
                        	 }
                        	 glEnd();

							 /* Draw other atom cap: */
                        	 glBegin(GL_POLYGON);
                        	 glNormal(1.0,0.0,0.0);
                        	 for(int i=0;i<=numLineVertices;++i)
                        		 glVertex(dataItem->lineVertices[2][i]);
                        	 glEnd();

                        	 // draw end sphere
							 glPushMatrix();
                        	 glTranslate(0.0,0.0,bondAxisLen);
                        	 glCallList(dataItem->cpkSphereDisplayListId);
							 glPopMatrix();

                        	 /* Reset OpenGL modelview matrix: */
                    		 glPopMatrix();
                		 }
            		 }	
        		 } // for each atom
				 glPopName();
			 } // for each residue
        } // if 
    	++structureIndex;
        glPopName();
		} // for each structure
    /* Reset OpenGL state: */
    glDisable(GL_COLOR_MATERIAL);
}

void ProteinRenderer::glDrawTube(GLContextData& contextData) const
{
    /* Set cartoon rendering parameters: */
    glMaterial(GL_FRONT,cartoonMaterial);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT,GL_DIFFUSE);

    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawTube)
        {
            switch(sfIt->structure->getStructureType())
            {
                case Protein::SecondaryStructure::COIL:
                    /* Render coil region: */
                    if (grayoutCartoon)
						glColor(grayOutColor);
					else
						glColor(sfIt->backboneColor);
                    for(int i=0;i<numCoilVertices;++i)
                    {
                        Scalar angle1=Scalar(i)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac1=Math::cos(angle1);
                        Scalar zFac1=Math::sin(angle1);
                        Scalar angle2=Scalar(i+1)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac2=Math::cos(angle2);
                        Scalar zFac2=Math::sin(angle2);
                        glBegin(GL_QUAD_STRIP);
                        for(int j=0;j<=sfIt->numCartoonSplineSamples;++j)
                        {
                            Vector v2=sfIt->cartoonX[j]*xFac2+sfIt->cartoonZ[j]*zFac2;
                            glNormal(v2);
                            glVertex(sfIt->cartoonP[j]+v2*coilRadius);
                            Vector v1=sfIt->cartoonX[j]*xFac1+sfIt->cartoonZ[j]*zFac1;
                            glNormal(v1);
                            glVertex(sfIt->cartoonP[j]+v1*coilRadius);
                        }
                        glEnd();
                    }
                    break;
                
                case Protein::SecondaryStructure::ALPHA_HELIX:
                {
                    /* Render alpha helix: */
                    if (grayoutCartoon)
						glColor(grayOutColor);
					else
						glColor(sfIt->backboneColor);

                    for(int i=0;i<numCoilVertices;++i)
                    {
                        Scalar angle1=Scalar(i)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac1=Math::cos(angle1);
                        Scalar zFac1=Math::sin(angle1);
                        Scalar angle2=Scalar(i+1)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac2=Math::cos(angle2);
                        Scalar zFac2=Math::sin(angle2);
                        glBegin(GL_QUAD_STRIP);
                        for(int j=0;j<=sfIt->numCartoonSplineSamples;++j)
                        {
                            Vector v2=sfIt->cartoonX[j]*xFac2+sfIt->cartoonZ[j]*zFac2;
                            glNormal(v2);
                            glVertex(sfIt->cartoonP[j]+v2*coilRadius);
                            Vector v1=sfIt->cartoonX[j]*xFac1+sfIt->cartoonZ[j]*zFac1;
                            glNormal(v1);
                            glVertex(sfIt->cartoonP[j]+v1*coilRadius);
                        }
                        glEnd();
                    }
                    break;
                }
                
                case Protein::SecondaryStructure::BETA_STRAND:
                {
                    /* Render beta strand: */
                    if (grayoutCartoon)
						glColor(grayOutColor);
					else
						glColor(sfIt->backboneColor);
               		int tailEnd=sfIt->numCartoonSplineSamples-4;
                         for(int i=0;i<numCoilVertices;++i)
                    {
                        Scalar angle1=Scalar(i)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac1=Math::cos(angle1);
                        Scalar zFac1=Math::sin(angle1);
                        Scalar angle2=Scalar(i+1)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac2=Math::cos(angle2);
                        Scalar zFac2=Math::sin(angle2);
                        Scalar angle3=Scalar(i)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac3=Math::cos(angle3)*Scalar(2);
                        Scalar zFac3=Math::sin(angle3)*Scalar(2);
                        glBegin(GL_QUAD_STRIP);
                      
						for(int j=0;j<=tailEnd;++j)
                        {
                            Vector v2=sfIt->cartoonX[j]*xFac2+sfIt->cartoonZ[j]*zFac2;
                            glNormal(v2);
                            glVertex(sfIt->cartoonP[j]+v2*coilRadius);
                            Vector v1=sfIt->cartoonX[j]*xFac1+sfIt->cartoonZ[j]*zFac1;
                            glNormal(v1);
                            glVertex(sfIt->cartoonP[j]+v1*coilRadius);
                        }
						glEnd();
                        glBegin(GL_QUAD_STRIP);
                        for(int j=tailEnd;j<=sfIt->numCartoonSplineSamples;++j)
                        {
                            int k = sfIt->numCartoonSplineSamples - j +1;
							Scalar sk = Scalar(k/2.0);
							Vector v2=sfIt->cartoonX[j]*xFac2+sfIt->cartoonZ[j]*zFac2;
                            glNormal(v2);
                            glVertex(sfIt->cartoonP[j]+v2*coilRadius*sk);
                            Vector v1=sfIt->cartoonX[j]*xFac1+sfIt->cartoonZ[j]*zFac1;
                            glNormal(v1);
                            glVertex(sfIt->cartoonP[j]+v1*coilRadius*sk);
                        }
						glEnd();
                    }
                    break;
                }
			default:
		    break;
            }
        }
        ++structureIndex;
        glPopName();
    }
    /* Reset OpenGL state: */
    glDisable(GL_COLOR_MATERIAL);
}

void ProteinRenderer::glDrawLines(GLContextData& contextData) const
{
    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawLine)
        {
            /* Set rendering parameters for this structure: */
            glLineWidth(sfIt->bondWidth);
            glColor(sfIt->bondColor);
            
            /* Iterate through all residues in this secondary structure: */
            for(const Protein::Residue* rPtr=sfIt->structure->residueBegin;rPtr!=sfIt->structure->residueEnd;rPtr=rPtr->succ)
            {
                /* Push residue (PDB) index onto the selection name stack: */
                glPushName(rPtr->residueIndex);
                
                /* Iterate through all atoms in this residue: */
                glBegin(GL_LINES);
                for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
                {
                    /* Iterate through all bonds for this atom: */
                    for(std::vector<Atom*>::const_iterator bIt=aPtr->getBonds().begin();bIt!=aPtr->getBonds().end();++bIt)
                    {
						if(static_cast<const Atom*>(aPtr)<*bIt)
                        {
						if (rPtr->getResType()==Protein::Residue::PHOBIC && showHydrophobic)
							glColor(phobicColor);
						else if (rPtr->getResType()==Protein::Residue::PHILIC && showHydrophilic) 
							glColor(philicColor);
						else if (rPtr->getResType()==Protein::Residue::DISULFIDE && showDisulfide) 
							glColor(disulfideColor);
						else
                           	glColor(elementColors[(*bIt)->getType()]);
                           	//glColor(elementColors[aPtr->getType()]);

                            /* Calculate bond axis for illumination: */
                            Vector bondAxis=(*bIt)->getPosition()-aPtr->getPosition();
                            glTexCoord(bondAxis.normalize());
                            glVertex(aPtr->getPosition());
                            glVertex((*bIt)->getPosition());
                        }
                    }
                }
                glEnd();
                
                glPopName();
            }
            ++structureIndex;
            glPopName();
        }
    }
}

void ProteinRenderer::glDrawBondCylinders(GLContextData& contextData) const
{
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    
    /* Set atom rendering parameters: */
    glMaterial(GL_FRONT,bondMaterial);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT,GL_DIFFUSE);
    
    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawBonds)
        {
            /* Iterate through all residues in this secondary structure: */
            for(const Protein::Residue* rPtr=sfIt->structure->residueBegin;rPtr!=sfIt->structure->residueEnd;rPtr=rPtr->succ)
            {
                /* Push residue (PDB) index onto the selection name stack: */
                glPushName(rPtr->residueIndex);

                /* Iterate through all atoms in this residue: */
                for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
                {
                    /* Iterate through all bonds for this atom: */
                    for(std::vector<Atom*>::const_iterator bIt=aPtr->getBonds().begin();bIt!=aPtr->getBonds().end();++bIt)
                    {
                        if(static_cast<const Atom*>(aPtr)<*bIt)
                        {
                            /* Draw a bicolor capped cylinder around the bond axis: */
                            Vector bondAxis=(*bIt)->getPosition()-aPtr->getPosition();
                            double bondAxisLen=Geometry::mag(bondAxis);
                            bondAxis/=Scalar(bondAxisLen);
                            
                            for(int i=0;i<=numBondVertices;++i)
                            {
                                dataItem->bondVertices[1][i][2]=GLfloat(bondAxisLen*0.5);
                                dataItem->bondVertices[2][i][2]=GLfloat(bondAxisLen);
                            }
                            
                            /* Set OpenGL modelview matrix to align z axis with bond axis: */
                            glPushMatrix();
                            glTranslate(aPtr->getPosition()-Point::origin);
                            Vector rotAxis(-bondAxis[1],bondAxis[0],Scalar(0));
                            Scalar rotAngle=Math::deg(Math::acos(bondAxis[2]));
                            glRotate(rotAngle,rotAxis);
                            
                            /* Draw root atom cap: */
                            glBegin(GL_POLYGON);

							if (rPtr->getResType()==Protein::Residue::PHOBIC && showHydrophobic)
								glColor(phobicColor);
							else if (rPtr->getResType()==Protein::Residue::PHILIC && showHydrophilic) 
								glColor(philicColor);
							else if (rPtr->getResType()==Protein::Residue::DISULFIDE && showDisulfide) 
								glColor(disulfideColor);
							else
                            	glColor(elementColors[aPtr->getType()]);

                            glNormal(0.0,0.0,-1.0);
                            for(int i=numBondVertices;i>=0;--i)
                                glVertex(dataItem->bondVertices[0][i]);
                            glEnd();
                            
                            /* Draw cylinder at root atom: */
                            glBegin(GL_QUAD_STRIP);
                            for(int i=0;i<=numBondVertices;++i)
                            {
                                glNormal(dataItem->bondNormals[i]);
                                glVertex(dataItem->bondVertices[1][i]);
                                glVertex(dataItem->bondVertices[0][i]);
                            }
                            glEnd();
                            
                            /* Draw cylinder at other atom: */
                            glBegin(GL_QUAD_STRIP);
							if (rPtr->getResType()==Protein::Residue::PHOBIC && showHydrophobic)
								glColor(phobicColor);
							else if (rPtr->getResType()==Protein::Residue::PHILIC && showHydrophilic) 
								glColor(philicColor);
							else if (rPtr->getResType()==Protein::Residue::DISULFIDE && showDisulfide) 
								glColor(disulfideColor);
							else
								glColor(elementColors[(*bIt)->getType()]);
                            for(int i=0;i<=numBondVertices;++i)
                            {
                                glNormal(dataItem->bondNormals[i]);
                                glVertex(dataItem->bondVertices[2][i]);
                                glVertex(dataItem->bondVertices[1][i]);
                            }
                            glEnd();
                            
                            /* Draw other atom cap: */
                            glBegin(GL_POLYGON);
                            glNormal(0.0,0.0,1.0);
                            for(int i=0;i<=numBondVertices;++i)
                                glVertex(dataItem->bondVertices[2][i]);
                            glEnd();
                            
                            /* Reset OpenGL modelview matrix: */
                            glPopMatrix();
                        }
                    }
                }
                
                glPopName();
            }
        }
        
        ++structureIndex;
        glPopName();
    }
    
    /* Reset OpenGL state: */
    glDisable(GL_COLOR_MATERIAL);
}

void ProteinRenderer::glDrawBackbone(GLContextData& contextData) const
{
    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawBackbone)
        {
            /* Set rendering parameters for this structure: */
            glLineWidth(sfIt->backboneWidth);
            //glColor(sfIt->backboneColor);
            
            /* Render residues in this structure independently: */
            for(const Protein::Residue* rPtr=sfIt->structure->residueBegin;rPtr!=sfIt->structure->residueEnd;rPtr=rPtr->succ)
            {
                /* Push residue (PDB) index onto the selection name stack: */
                glPushName(rPtr->residueIndex);
				
				if (rPtr->getResType()==Protein::Residue::PHOBIC && showHydrophobic)
					glColor(phobicColor);
				else if (rPtr->getResType()==Protein::Residue::PHILIC && showHydrophilic) 
					glColor(philicColor);
				else
            		glColor(sfIt->backboneColor);
                
                glBegin(GL_LINE_STRIP);
                if(rPtr->pred!=0)
                {
                    glVertex(Geometry::mid(rPtr->pred->backbone[rPtr->pred->backbone.size()-1]->getPosition(),rPtr->backbone[0]->getPosition()));
                    #if 0
                    /* Draw midpoint between last atom in previous residue and first atom in this one: */
                    Point mid=Point::origin;
                    mid.affineAdd(rPtr->pred->backbone[rPtr->pred->backbone.size()-1]->getPosition());
                    mid.affineAdd(rPtr->backbone[0]->getPosition());
                    mid.affineDivide(2);
                    glVertex(mid);
                    #endif
                }
                for(std::vector<Protein::ChainAtom*>::const_iterator bbIt=rPtr->backbone.begin();bbIt!=rPtr->backbone.end();++bbIt)
                    glVertex((*bbIt)->getPosition());
                if(rPtr->succ!=0)
                {
                    glVertex(Geometry::mid(rPtr->backbone[rPtr->backbone.size()-1]->getPosition(),rPtr->succ->backbone[0]->getPosition()));
                    #if 0
                    /* Draw midpoint between last atom in this residue and first atom in the next one: */
                    Point mid=Point::origin;
                    mid.affineAdd(rPtr->backbone[rPtr->backbone.size()-1]->getPosition());
                    mid.affineAdd(rPtr->succ->backbone[0]->getPosition());
                    mid.affineDivide(2);
                    glVertex(mid);
                    #endif
                }
                glEnd();
                
                glPopName();
            }
        }
        
        ++structureIndex;
        glPopName();
    }
    
    /* Reset OpenGL state: */
}

void ProteinRenderer::glDrawBackboneRibbon(GLContextData& contextData) const
{
    /* Set ribbon rendering parameters: */
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
    glMaterial(GL_FRONT_AND_BACK,backboneRibbonMaterial);
        
    /* Create a spline evaluation cache: */
    Geometry::SplineCurve<Scalar,3>::EvaluationCache* cache=backboneRibbonSpline->createEvaluationCache();
    
    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawBackboneRibbon)
        {
            /* Set rendering parameters for this structure: */
            glColor(sfIt->backboneColor);
            glBegin(GL_QUAD_STRIP);
            
            /* Render all samples: */
            Scalar sideFactor=sfIt->backboneRibbonWidth*Scalar(0.5);
            for(int i=0;i<=sfIt->backboneRibbonNumSamples;++i)
            {
                /* Evaluate the spline curve: */
                Scalar u=Scalar(i)*(sfIt->backboneRibbonParameterMax-sfIt->backboneRibbonParameterMin)/Scalar(sfIt->backboneRibbonNumSamples)+sfIt->backboneRibbonParameterMin;
                Vector deriv1,deriv2;
                Point point=backboneRibbonSpline->evaluate(u,cache,deriv1,deriv2);

                /* Calculate a Frenet frame centered at the evaluated point: */
                deriv1.normalize();
                Vector side=Geometry::cross(deriv1,deriv2);
                side.normalize();
                deriv2=Geometry::cross(side,deriv1);
                glNormal(deriv2);
                side*=sideFactor;
                glVertex(point-side);
                glVertex(point+side);
            }
            
            glEnd();
        }
        
        ++structureIndex;
        glPopName();
    }

    /* Clean up: */
    delete cache;
    
    /* Reset OpenGL state: */
    glEnable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
    glDisable(GL_COLOR_MATERIAL);
}

void ProteinRenderer::glDrawCartoon(GLContextData& contextData) const
{
    /* Set cartoon rendering parameters: */
    glMaterial(GL_FRONT,cartoonMaterial);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT,GL_DIFFUSE);
    
    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawCartoon)
        {
            switch(sfIt->structure->getStructureType())
            {
                case Protein::SecondaryStructure::COIL:
                    /* Render coil region: */
                    if (grayoutCartoon)
						glColor(grayOutColor);
					else
                    	glColor(sfIt->backboneColor);
                    for(int i=0;i<numCoilVertices;++i)
                    {
                        Scalar angle1=Scalar(i)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac1=Math::cos(angle1);
                        Scalar zFac1=Math::sin(angle1);
                        Scalar angle2=Scalar(i+1)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(numCoilVertices);
                        Scalar xFac2=Math::cos(angle2);
                        Scalar zFac2=Math::sin(angle2);
                        glBegin(GL_QUAD_STRIP);
                        for(int j=0;j<=sfIt->numCartoonSplineSamples;++j)
                        {
                            Vector v2=sfIt->cartoonX[j]*xFac2+sfIt->cartoonZ[j]*zFac2;
                            glNormal(v2);
                            glVertex(sfIt->cartoonP[j]+v2*coilRadius);
                            Vector v1=sfIt->cartoonX[j]*xFac1+sfIt->cartoonZ[j]*zFac1;
                            glNormal(v1);
                            glVertex(sfIt->cartoonP[j]+v1*coilRadius);
                        }
                        glEnd();
                    }
                    break;

                case Protein::SecondaryStructure::ALPHA_HELIX:
                {
                    /* Render alpha helix: */
                    if (grayoutCartoon)
						glColor(grayOutColor);
					else
                    	glColor(sfIt->backboneColor);
                    glBegin(GL_QUADS);
                    Vector x0=sfIt->cartoonX[0]*alphaHelixWidth;
                    Vector z0=sfIt->cartoonZ[0]*alphaHelixThickness;
                    glNormal(-sfIt->cartoonY[0]);
                    glVertex(sfIt->cartoonP[0]-z0-x0);
                    glVertex(sfIt->cartoonP[0]-z0+x0);
                    glVertex(sfIt->cartoonP[0]+z0+x0);
                    glVertex(sfIt->cartoonP[0]+z0-x0);
                    Vector xn=sfIt->cartoonX[sfIt->numCartoonSplineSamples]*alphaHelixWidth;
                    Vector zn=sfIt->cartoonZ[sfIt->numCartoonSplineSamples]*alphaHelixThickness;
                    glNormal(sfIt->cartoonY[sfIt->numCartoonSplineSamples]);
                    glVertex(sfIt->cartoonP[sfIt->numCartoonSplineSamples]-zn-xn);
                    glVertex(sfIt->cartoonP[sfIt->numCartoonSplineSamples]+zn-xn);
                    glVertex(sfIt->cartoonP[sfIt->numCartoonSplineSamples]+zn+xn);
                    glVertex(sfIt->cartoonP[sfIt->numCartoonSplineSamples]-zn+xn);
                    glEnd();
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0;i<=sfIt->numCartoonSplineSamples;++i)
                    {
                        Vector x=sfIt->cartoonX[i]*alphaHelixWidth;
                        Vector z=sfIt->cartoonZ[i]*alphaHelixThickness;
                        glNormal(sfIt->cartoonZ[i]);
                        glVertex(sfIt->cartoonP[i]+z-x);
                        glVertex(sfIt->cartoonP[i]+z+x);
                    }
                    glEnd();
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0;i<=sfIt->numCartoonSplineSamples;++i)
                    {
                        Vector x=sfIt->cartoonX[i]*alphaHelixWidth;
                        Vector z=sfIt->cartoonZ[i]*alphaHelixThickness;
                        glNormal(-sfIt->cartoonZ[i]);
                        glVertex(sfIt->cartoonP[i]-z+x);
                        glVertex(sfIt->cartoonP[i]-z-x);
                    }
                    glEnd();
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0;i<=sfIt->numCartoonSplineSamples;++i)
                    {
                        Vector x=sfIt->cartoonX[i]*alphaHelixWidth;
                        Vector z=sfIt->cartoonZ[i]*alphaHelixThickness;
                        glNormal(sfIt->cartoonX[i]);
                        glVertex(sfIt->cartoonP[i]+z+x);
                        glVertex(sfIt->cartoonP[i]-z+x);
                    }
                    glEnd();
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0;i<=sfIt->numCartoonSplineSamples;++i)
                    {
                        Vector x=sfIt->cartoonX[i]*alphaHelixWidth;
                        Vector z=sfIt->cartoonZ[i]*alphaHelixThickness;
                        glNormal(-sfIt->cartoonX[i]);
                        glVertex(sfIt->cartoonP[i]-z-x);
                        glVertex(sfIt->cartoonP[i]+z-x);
                    }
                    glEnd();
                    break;
                }
                
                case Protein::SecondaryStructure::BETA_STRAND:
                {
                    /* Render beta strand: */
                    if (grayoutCartoon)
						glColor(grayOutColor);
					else
                    	glColor(sfIt->backboneColor);
                    glBegin(GL_QUADS);
                    Vector x0=sfIt->cartoonX[0]*betaStrandWidth;
                    Vector z0=sfIt->cartoonZ[0]*betaStrandThickness;
                    glNormal(-sfIt->cartoonY[0]);
                    glVertex(sfIt->cartoonP[0]-z0-x0);
                    glVertex(sfIt->cartoonP[0]-z0+x0);
                    glVertex(sfIt->cartoonP[0]+z0+x0);
                    glVertex(sfIt->cartoonP[0]+z0-x0);
                    glEnd();
                    int tailEnd=sfIt->numCartoonSplineSamples-cartoonSampleDensity/2;
                    Scalar widthFactor=betaStrandHeadWidth/Scalar(cartoonSampleDensity/2);
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0;i<=tailEnd;++i)
                    {
                        Vector x=sfIt->cartoonX[i]*betaStrandWidth;
                        Vector z=sfIt->cartoonZ[i]*betaStrandThickness;
                        glNormal(sfIt->cartoonZ[i]);
                        glVertex(sfIt->cartoonP[i]+z-x);
                        glVertex(sfIt->cartoonP[i]+z+x);
                    }
                    for(int i=tailEnd;i<=sfIt->numCartoonSplineSamples;++i)
                    {
                        Scalar width=betaStrandWidth*Scalar(sfIt->numCartoonSplineSamples-i)*widthFactor;
                        Vector x=sfIt->cartoonX[i]*width;
                        Vector z=sfIt->cartoonZ[i]*betaStrandThickness;
                        glNormal(sfIt->cartoonZ[i]);
                        glVertex(sfIt->cartoonP[i]+z-x);
                        glVertex(sfIt->cartoonP[i]+z+x);
                    }
                    glEnd();
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0;i<=tailEnd;++i)
                    {
                        Vector x=sfIt->cartoonX[i]*betaStrandWidth;
                        Vector z=sfIt->cartoonZ[i]*betaStrandThickness;
                        glNormal(-sfIt->cartoonZ[i]);
                        glVertex(sfIt->cartoonP[i]-z+x);
                        glVertex(sfIt->cartoonP[i]-z-x);
                    }
                    for(int i=tailEnd;i<=sfIt->numCartoonSplineSamples;++i)
                    {
                        Scalar width=betaStrandWidth*Scalar(sfIt->numCartoonSplineSamples-i)*widthFactor;
                        Vector x=sfIt->cartoonX[i]*width;
                        Vector z=sfIt->cartoonZ[i]*betaStrandThickness;
                        glNormal(-sfIt->cartoonZ[i]);
                        glVertex(sfIt->cartoonP[i]-z+x);
                        glVertex(sfIt->cartoonP[i]-z-x);
                    }
                    glEnd();
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0;i<=tailEnd;++i)
                    {
                        Vector x=sfIt->cartoonX[i]*betaStrandWidth;
                        Vector z=sfIt->cartoonZ[i]*betaStrandThickness;
                        glNormal(sfIt->cartoonX[i]);
                        glVertex(sfIt->cartoonP[i]+z+x);
                        glVertex(sfIt->cartoonP[i]-z+x);
                        if(i==tailEnd)
                        {
                            glNormal(-sfIt->cartoonY[i]);
                            glVertex(sfIt->cartoonP[i]+z+x);
                            glVertex(sfIt->cartoonP[i]-z+x);
                        }
                    }
                    for(int i=tailEnd;i<=sfIt->numCartoonSplineSamples;++i)
                    {
                        Scalar width=betaStrandWidth*Scalar(sfIt->numCartoonSplineSamples-i)*widthFactor;
                        Vector x=sfIt->cartoonX[i]*width;
                        Vector z=sfIt->cartoonZ[i]*betaStrandThickness;
                        if(i==tailEnd)
                        {
                            glNormal(-sfIt->cartoonY[i]);
                            glVertex(sfIt->cartoonP[i]+z+x);
                            glVertex(sfIt->cartoonP[i]-z+x);
                        }
                        glNormal(sfIt->cartoonY[i]+sfIt->cartoonX[i]);
                        glVertex(sfIt->cartoonP[i]+z+x);
                        glVertex(sfIt->cartoonP[i]-z+x);
                    }
                    glEnd();
                    glBegin(GL_QUAD_STRIP);
                    for(int i=0;i<=tailEnd;++i)
                    {
                        Vector x=sfIt->cartoonX[i]*betaStrandWidth;
                        Vector z=sfIt->cartoonZ[i]*betaStrandThickness;
                        glNormal(-sfIt->cartoonX[i]);
                        glVertex(sfIt->cartoonP[i]-z-x);
                        glVertex(sfIt->cartoonP[i]+z-x);
                        if(i==tailEnd)
                        {
                            glNormal(-sfIt->cartoonY[i]);
                            glVertex(sfIt->cartoonP[i]-z-x);
                            glVertex(sfIt->cartoonP[i]+z-x);
                        }
                    }
                    for(int i=tailEnd;i<=sfIt->numCartoonSplineSamples;++i)
                    {
                        Scalar width=betaStrandWidth*Scalar(sfIt->numCartoonSplineSamples-i)*widthFactor;
                        Vector x=sfIt->cartoonX[i]*width;
                        Vector z=sfIt->cartoonZ[i]*betaStrandThickness;
                        if(i==tailEnd)
                        {
                            glNormal(-sfIt->cartoonY[i]);
                            glVertex(sfIt->cartoonP[i]-z-x);
                            glVertex(sfIt->cartoonP[i]+z-x);
                        }
                        glNormal(sfIt->cartoonY[i]-sfIt->cartoonX[i]);
                        glVertex(sfIt->cartoonP[i]-z-x);
                        glVertex(sfIt->cartoonP[i]+z-x);
                    }
                    glEnd();
                    break;
                }
		
		default:
		    break;
            }
        }
        
        ++structureIndex;
        glPopName();
    }

    /* Reset OpenGL state: */
    glDisable(GL_COLOR_MATERIAL);
}

void ProteinRenderer::glDrawHydrogenBonds(GLContextData& contextData) const
{
    /* Set hydrogen bond rendering parameters: */
    glLineWidth(hydrogenBondWidth);
    glColor(hydrogenBondColor);
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(4,0xAAAA);
    glBegin(GL_LINES);
    
    /* Iterate through all possible combinations of opposite dipoles: */
    for(std::vector<Protein::Dipole>::const_iterator nh=protein->amides.begin();nh!=protein->amides.end();++nh)
        for(std::vector<Protein::Dipole>::const_iterator co=protein->carboxyls.begin();co!=protein->carboxyls.end();++co)
        {
            if(formHydrogenBond(*nh,*co))
            {
                glVertex(nh->getMinorAtom()->getPosition());
                glVertex(co->getMinorAtom()->getPosition());
            }
        }

    glEnd();
    glDisable(GL_LINE_STIPPLE);
}

void ProteinRenderer::glDrawHydrogenBondSites(GLContextData& contextData) const
{
    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawHydrogenBondSites)
        {
            #if 1
            /* Set rendering parameters for this structure: */
            glPointSize(sfIt->hydrogenBondSiteDiameter);
            glLineWidth(sfIt->hydrogenBondSiteWidth);
            
            /* Iterate through all amide groups in this secondary structure: */
            glColor(sfIt->amideColor);
            for(std::vector<Protein::Dipole>::const_iterator nh=sfIt->structure->amidesBegin;nh!=sfIt->structure->amidesEnd;++nh)
            {
                glPushName(nh->getMajorAtom()->residue->residueIndex);
                Point bondSite=nh->getBondSite();
                glBegin(GL_POINTS);
                glVertex(bondSite);
                glEnd();
                glBegin(GL_LINES);
                glVertex(nh->getMajorAtom()->getPosition());
                glVertex(bondSite);
                glEnd();
                glPopName();
            }
            
            /* Iterate through all carboxyl groups in this secondary structure: */
            glColor(sfIt->carboxylColor);
            for(std::vector<Protein::Dipole>::const_iterator co=sfIt->structure->carboxylsBegin;co!=sfIt->structure->carboxylsEnd;++co)
            {
                glPushName(co->getMajorAtom()->residue->residueIndex);
                Point bondSite=co->getBondSite();
                glBegin(GL_POINTS);
                glVertex(bondSite);
                glEnd();
                glBegin(GL_LINES);
                glVertex(co->getMajorAtom()->getPosition());
                glVertex(bondSite);
                glEnd();
                glPopName();
            }
            #else
            /* Set point rendering parameters for this structure: */
            glPointSize(sfIt->hydrogenBondSiteDiameter);
            glBegin(GL_POINTS);

            /* Iterate through all carboxyl groups in this secondary structure: */
            glColor(sfIt->carboxylColor);
            for(std::vector<Protein::Dipole>::const_iterator co=sfIt->structure->carboxylsBegin;co!=sfIt->structure->carboxylsEnd;++co)
                glVertex(co->getBondSite());

            glEnd();

            /* Set line rendering parameters for this structure: */
            glLineWidth(sfIt->hydrogenBondSiteWidth);
            glBegin(GL_LINES);

            /* Iterate through all amide groups in this secondary structure: */
            glColor(sfIt->amideColor);
            for(std::vector<Protein::Dipole>::const_iterator nh=sfIt->structure->amidesBegin;nh!=sfIt->structure->amidesEnd;++nh)
            {
                glVertex(nh->getMajorAtom()->getPosition());
                glVertex(nh->getBondSite());
            }

            /* Iterate through all carboxyl groups in this secondary structure: */
            glColor(sfIt->carboxylColor);
            for(std::vector<Protein::Dipole>::const_iterator co=sfIt->structure->carboxylsBegin;co!=sfIt->structure->carboxylsEnd;++co)
            {
                glVertex(co->getMajorAtom()->getPosition());
                glVertex(co->getBondSite());
            }

            glEnd();
            #endif
        }
        
        ++structureIndex;
        glPopName();
    }
}

void ProteinRenderer::glDrawHydrogenCages(GLContextData& contextData) const
{
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    
    /* Render each secondary structure independently: */
    GLuint structureIndex=0;
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        /* Push structure index onto the selection name stack: */
        glPushName(structureIndex);
        
        if(sfIt->drawHydrogenCages)
        {
            /* Set point rendering parameters for this structure: */
            glPointSize(sfIt->hydrogenBondSiteDiameter);
            glBegin(GL_POINTS);

            /* Iterate through all carboxyl groups in this secondary structure: */
            glColor(sfIt->carboxylColor);
            for(std::vector<Protein::Dipole>::const_iterator co=sfIt->structure->carboxylsBegin;co!=sfIt->structure->carboxylsEnd;++co)
                glVertex(co->getMinorAtom()->getPosition());

            glEnd();
            
            /* Set line rendering parameters for this structure: */
            glLineWidth(sfIt->hydrogenBondSiteWidth);
            glBegin(GL_LINES);

            /* Iterate through all carboxyl groups in this secondary structure: */
            glColor(sfIt->carboxylColor);
            for(std::vector<Protein::Dipole>::const_iterator co=sfIt->structure->carboxylsBegin;co!=sfIt->structure->carboxylsEnd;++co)
            {
                glVertex(co->getMajorAtom()->getPosition());
                glVertex(co->getMinorAtom()->getPosition());
            }

            glEnd();
            
            /* Set hydrogen cage rendering parameters for this structure: */
            glLineWidth(sfIt->hydrogenCageWidth);
            glColor(sfIt->hydrogenCageColor);
            for(std::vector<Protein::Dipole>::const_iterator nh=sfIt->structure->amidesBegin;nh!=sfIt->structure->amidesEnd;++nh)
            {
                /* Update modelview matrix to position cage: */
                glPushMatrix();
                Scalar matrix[4][4]; // Column-major transformation matrix
                Vector y=nh->getMinorAtom()->getPosition()-nh->getMajorAtom()->getPosition();
                y.normalize();
                Vector x=nh->getAlphaCarbon()->getPosition()-nh->getMajorAtom()->getPosition();
                Vector z=Geometry::cross(x,y);
                z.normalize();
                x=Geometry::cross(y,z);
                for(int i=0;i<3;++i)
                {
                    matrix[0][i]=x[i];
                    matrix[1][i]=y[i];
                    matrix[2][i]=z[i];
                    matrix[3][i]=nh->getMinorAtom()->getPosition()[i];
                }
                matrix[0][3]=matrix[1][3]=matrix[2][3]=Scalar(0);
                matrix[3][3]=Scalar(1);
                glMultMatrix(&matrix[0][0]);
                
                /* Render the cage: */
                if(sfIt->hydrogenCageLarge)
                    glCallList(dataItem->hydrogenCageLargeDisplayListId);
                else
                    glCallList(dataItem->hydrogenCageSmallDisplayListId);
                
                glPopMatrix();
            }
        }
        
        ++structureIndex;
        glPopName();
    }
}

void ProteinRenderer::glDrawCollisions(GLContextData& contextData) const
{
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);

    /* Set collision sphere rendering parameters: */
    glMaterial(GL_FRONT,collisionSphereMaterial);
    
    /* Calculate offsets for the 27 neighbours of a cell: */
    int cellOffsets[27];
    int* coPtr=cellOffsets;
    for(int z=-1;z<=1;++z)
        for(int y=-1;y<=1;++y)
            for(int x=-1;x<=1;++x,++coPtr)
                *coPtr=(z*numCollisionCells[1]+y)*numCollisionCells[0]+x;

    /* Iterate through all atoms: */
    const AtomListItem* aliPtr=atomListItems;
    for(int i=0;i<protein->getNumAtoms();++i,++aliPtr)
    {
        const Protein::ChainAtom* aPtr1=aliPtr->atom;
        /* Check all 27 neighbouring cells: */
        for(int j=0;j<27;++j)
        {
            for(AtomListItem* cellAliPtr=aliPtr->cell[cellOffsets[j]].atoms;cellAliPtr!=0;cellAliPtr=cellAliPtr->cellSucc)
            {
                if(cellAliPtr>aliPtr)
                {
                    const Protein::ChainAtom* aPtr2=cellAliPtr->atom;
                    double radiusSum=(aPtr1->getCovalentRadius()+aPtr2->getCovalentRadius())*0.75;
                    if(Geometry::sqrDist(aPtr1->getPosition(),aPtr2->getPosition())<Math::sqr(radiusSum))
                    {
                        /* Draw a big red sphere at the collision's center point: */
                        double penetrationDepth=(radiusSum-Geometry::dist(aPtr1->getPosition(),aPtr2->getPosition()))*2.0;
                        glPushMatrix();
                        glTranslate(Geometry::mid(aPtr1->getPosition(),aPtr2->getPosition())-Position::origin);
                        glScale(penetrationDepth,penetrationDepth,penetrationDepth);
                        glCallList(dataItem->collisionSphereDisplayListId);
                        glPopMatrix();
                    }
                }
            }
        }
    }
}

ProteinRenderer::ProteinRenderer(const ConfigurationFile::SectionIterator& sConfigFileSection,const Protein* sProtein)
    :configFileSection(sConfigFileSection),protein(sProtein),
     atomListItems(new AtomListItem[protein->getNumAtoms()]),collisionCells(0),
     boundingBox(Geometry::Box<Scalar,3>::empty),
     drawAtoms(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawAtoms",false)),
     atomMaterial(configFileSection.retrieveValue<GLMaterial,ValueCoder<GLMaterial> >("./atomMaterial",GLMaterial(Color(0.1,0.1,0.1),Color(1.0,1.0,1.0),Color(1.0,1.0,1.0),25.0))),
     atomTesselation(configFileSection.retrieveValue<int,ValueCoder<int> >("./atomTesselation",4)),
     mapAtomValues(false),
     atomValues(new float[protein->getNumAtoms()]),
     drawBonds(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawBonds",false)),
     bondMaterial(configFileSection.retrieveValue<GLMaterial,ValueCoder<GLMaterial> >("./bondMaterial",GLMaterial(Color(0.1,0.1,0.1),Color(1.0,0.0,0.0),Color(1.0,1.0,1.0),25.0))),
     bondMaterialVersion(0),
     numLineVertices(configFileSection.retrieveValue<int,ValueCoder<int> >("./numLineVertices",6)),
     lineRadius(configFileSection.retrieveValue<float,ValueCoder<float> >("./lineRadius",0.15f)),
 	 cpkSphereRadius(configFileSection.retrieveValue<float,ValueCoder<float> >("./cpkSphereRadius",0.3f)),
     drawCPK(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawCPK",false)),
     drawTube(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawTube",false)),
     drawLine(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawLine",false)),
     numBondVertices(configFileSection.retrieveValue<int,ValueCoder<int> >("./numBondVertices",24)),
     bondRadius(configFileSection.retrieveValue<float,ValueCoder<float> >("./bondRadius",0.25f)),
     alphaHelixColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./alphaHelixColor",Color(1.0,0.0,0.7))),
     betaStrandColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./betaStrandColor",Color(0.7,0.0,1.0))),
     coilColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./coilColor",Color(1.0,0.0,1.0))),
     highlightColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./highlightColor",Color(1.0,1.0,0.0))),
     drawBackbone(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawBackbone",true)),
     drawBackboneRibbon(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawBackboneRibbon",false)),
     backboneRibbonMaterial(configFileSection.retrieveValue<GLMaterial,ValueCoder<GLMaterial> >("./backboneRibbonMaterial",GLMaterial(Color(0.05,0.05,0.05),Color(1.0,0.0,1.0),Color(1.0,1.0,1.0),50.0))),
     backboneRibbonUseAllAtoms(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./backboneRibbonUseAllAtoms",false)),
     backboneRibbonDegree(configFileSection.retrieveValue<int,ValueCoder<int> >("./backboneRibbonDegree",3)),
     backboneRibbonSampleDensity(configFileSection.retrieveValue<int,ValueCoder<int> >("./backboneRibbonSampleDensity",16)),
     backboneRibbonSpline(0),
     drawCartoon(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawCartoon",true)),
     grayoutCartoon(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./grayOutCartoon",false)),
     grayOutColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./grayOutColor",Color(0.4,0.4,0.4))),
     cartoonMaterial(configFileSection.retrieveValue<GLMaterial,ValueCoder<GLMaterial> >("./cartoonMaterial",GLMaterial(Color(0.05,0.05,0.05),Color(1.0,0.0,1.0),Color(1.0,1.0,1.0),50.0))),
     cartoonDegree(configFileSection.retrieveValue<int,ValueCoder<int> >("./cartoonDegree",3)),
     cartoonSampleDensity(configFileSection.retrieveValue<int,ValueCoder<int> >("./cartoonSampleDensity",16)),
     alphaHelixWidth(configFileSection.retrieveValue<Scalar,ValueCoder<Scalar> >("./alphaHelixWidth",Scalar(1.5))),
     alphaHelixThickness(configFileSection.retrieveValue<Scalar,ValueCoder<Scalar> >("./alphaHelixThickness",Scalar(0.5))),
     betaStrandWidth(configFileSection.retrieveValue<Scalar,ValueCoder<Scalar> >("./betaStrandWidth",Scalar(1.5))),
     betaStrandThickness(configFileSection.retrieveValue<Scalar,ValueCoder<Scalar> >("./betaStrandThickness",Scalar(0.5))),
     betaStrandHeadWidth(configFileSection.retrieveValue<Scalar,ValueCoder<Scalar> >("./betaStrandHeadWidth",Scalar(1.333))),
     numCoilVertices(configFileSection.retrieveValue<int,ValueCoder<int> >("./numCoilVertices",12)),
     coilRadius(configFileSection.retrieveValue<Scalar,ValueCoder<Scalar> >("./coilRadius",Scalar(0.333))),
     drawHydrogenBonds(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawHydrogenBonds",true)),
     hydrogenBondWidth(configFileSection.retrieveValue<float,ValueCoder<float> >("./hydrogenBondWidth",3.0f)),
     hydrogenBondColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./hydrogenBondColor",Color(1.0,0.9,0.0))),
     drawHydrogenBondSites(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawHydrogenBondSites",false)),
     drawHydrogenCages(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawHydrogenCages",false)),
     drawCollisions(configFileSection.retrieveValue<bool,ValueCoder<bool> >("./drawCollisions",false)),
     collisionSphereMaterial(configFileSection.retrieveValue<GLMaterial,ValueCoder<GLMaterial> >("./collisionSphereMaterial",GLMaterial(Color(0.1,0.1,0.1),Color(1.0,0.0,0.0),Color(1.0,1.0,1.0),25.0))),
     collisionSphereTesselation(configFileSection.retrieveValue<int,ValueCoder<int> >("./collisionSphereTesselation",4)),
     phobicColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./hydrophobicColor",Color(1.0,1.0,1.0))),
     philicColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./hydrophilicColor",Color(0.5,0.5,1.0))),
     disulfideColor(configFileSection.retrieveValue<Color,ValueCoder<Color> >("./disulfideColor",Color(.75,0.5,0.0))),
	 showHydrophobic(false),showHydrophilic(false),showDisulfide(false),
	 drawResidueName(false),drawResidueAngles(false),drawAtomNames(false),
	 updateEnergyValue(false),settime(false),energyValue(0.0)
{
    /* Initialize atom list: */
    AtomListItem* aliPtr=atomListItems;
    for(const Protein::ChainAtom* aPtr=protein->atoms;aPtr!=0;aPtr=aPtr->succ,++aliPtr)
        aliPtr->atom=aPtr;

    /* Initialize atom mapping values: */
    for(int i=0;i<protein->getNumAtoms();++i)
        atomValues[i]=0.0f;

    /* Create atom color map: */
    static GLColorMap::Color mapColors[3]={GLColorMap::Color(0.0,1.0,0.0),GLColorMap::Color(1.0,1.0,0.0),GLColorMap::Color(1.0,0.0,0.0)};
    atomColorMap=GLColorMap(3,mapColors,0.0,1.0);

    /* Initialize per-structure rendering parameters: */
    updateStructureFlags();

    /* Initialize bond illuminator: */
    bondIlluminator.setSceneCenter(GLVector<GLfloat,3>(protein->calcCentroid().getComponents()));
    bondIlluminator.enableAutoView();
    bondIlluminator.enableAutoLight(GL_LIGHT0);
	text = new GLText;
   	text->setWinVars(renderWindow->w(), renderWindow->h());
}

ProteinRenderer::~ProteinRenderer(void)
{
    delete[] atomListItems;
    delete[] collisionCells;
    delete backboneRibbonSpline;
	delete text;
}


void ProteinRenderer::clearContext(GLContextData &contextData)
{
    /* Get a pointer to the context data: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);

    /* At this point we should also remove the context data entry itself from the context! */
    contextData.removeDataItem(this);

    /* Delete all associated display lists and array data: */
    delete dataItem;
}


bool ProteinRenderer::isContextInitialized(GLContextData &contextData) const
{
    return ( contextData.retrieveDataItem<DataItem>(this) );
}


void ProteinRenderer::initContext(GLContextData& contextData)
{
    /* Create a new context entry: */
    DataItem* dataItem=new DataItem(numBondVertices,numLineVertices);
    contextData.addDataItem(this,dataItem);

    /* Create a sphere for collision visualization: */
    glNewList(dataItem->collisionSphereDisplayListId,GL_COMPILE);
    drawSphere(1.0,collisionSphereTesselation);
    glEndList();
	
	/* Create a small sphere for cpk visualization: */
    glNewList(dataItem->cpkSphereDisplayListId,GL_COMPILE);
    drawSphere(cpkSphereRadius, atomTesselation);
    glEndList();

    /* Create atom spheres for all common elements: */
    const int numElements=8;
    const Atom::Element elements[numElements]={Atom::H,Atom::C,Atom::N,Atom::O,Atom::P,Atom::S,Atom::Cl,Atom::Zn};
    for(int i=0;i<numElements;++i)
    {
        glNewList(dataItem->atomSphereDisplayListBaseId+elements[i],GL_COMPILE);
        drawSphere(Atom::getVanDerWaalsRadius(elements[i]),atomTesselation);
        glEndList();
    }

    #if USE_BONDCYLINDER_VERTEXARRAY
    /* Create bond rendering cylinder vertices: */
    for(int ring=0;ring<3;++ring)
        for(int i=0;i<numBondVertices;++i)
        {
            GLfloat angle=2.0f*GLfloat(Math::pi)*GLfloat(i)/GLfloat(numBondVertices);
            dataItem->bondVertices[ring*numBondVertices+i].normal[0]=Math::cos(angle);
            dataItem->bondVertices[ring*numBondVertices+i].normal[1]=Math::sin(angle);
            dataItem->bondVertices[ring*numBondVertices+i].normal[2]=0.0f;
            dataItem->bondVertices[ring*numBondVertices+i].vertex[0]=Math::cos(angle);
            dataItem->bondVertices[ring*numBondVertices+i].vertex[1]=Math::sin(angle);
            dataItem->bondVertices[ring*numBondVertices+i].vertex[2]=GLfloat(ring);
        }

    /* Create bond rendering cylinder index arrays: */
    for(int i=0;i<numBondVertices;++i)
        dataItem->bondIndices[0][i]=numBondVertices-1-i;
    for(int i=0;i<numBondVertices;++i)
    {
        dataItem->bondIndices[1][2*i+0]=1*numBondVertices+i;
        dataItem->bondIndices[1][2*i+1]=0*numBondVertices+i;
    }
    dataItem->bondIndices[1][2*numBondVertices+0]=1*numBondVertices;
    dataItem->bondIndices[1][2*numBondVertices+1]=0*numBondVertices;
    for(int i=0;i<numBondVertices;++i)
    {
        dataItem->bondIndices[2][2*i+0]=2*numBondVertices+i;
        dataItem->bondIndices[2][2*i+1]=1*numBondVertices+i;
    }
    dataItem->bondIndices[2][2*numBondVertices+0]=2*numBondVertices;
    dataItem->bondIndices[2][2*numBondVertices+1]=1*numBondVertices;
    for(int i=0;i<numBondVertices;++i)
        dataItem->bondIndices[3][i]=2*numBondVertices+i;
    #else
    /* Create bond rendering cylinders: */
    for(int i=0;i<=numBondVertices;++i)
    {
        GLfloat angle=2.0f*Math::Constants<GLfloat>::pi*GLfloat(i)/GLfloat(numBondVertices);
        dataItem->bondNormals[i][0]=Math::cos(angle);
        dataItem->bondNormals[i][1]=Math::sin(angle);
        dataItem->bondNormals[i][2]=0.0f;
    }
    for(int ring=0;ring<3;++ring)
    {
        for(int i=0;i<=numBondVertices;++i)
        {
            GLfloat angle=2.0f*Math::Constants<GLfloat>::pi*GLfloat(i)/GLfloat(numBondVertices);
            dataItem->bondVertices[ring][i][0]=Math::cos(angle)*bondRadius;
            dataItem->bondVertices[ring][i][1]=Math::sin(angle)*bondRadius;
            dataItem->bondVertices[ring][i][2]=0.0f;
        }
    }
    #endif

    /* Create line rendering segmentss: */
    for(int i=0;i<=numLineVertices;++i)
    {
        GLfloat
		angle=2.0f*Math::Constants<GLfloat>::pi*GLfloat(i)/GLfloat(numLineVertices);
        dataItem->lineNormals[i][0]=Math::cos(angle);
        dataItem->lineNormals[i][1]=Math::sin(angle);
        dataItem->lineNormals[i][2]=0.0f;
    }
    for(int ring=0;ring<3;++ring)
    {
        for(int i=0;i<=numLineVertices;++i)
        {
            GLfloat
			angle=2.0f*Math::Constants<GLfloat>::pi*GLfloat(i)/GLfloat(numLineVertices);
            dataItem->lineVertices[ring][i][0]=Math::cos(angle)*lineRadius;
            dataItem->lineVertices[ring][i][1]=Math::sin(angle)*lineRadius;
            dataItem->lineVertices[ring][i][2]=0.0f;
        }
    }

    /* Create hydrogen cage template: */
    Scalar a=Protein::Dipole::NHdist;
    Scalar b=Protein::Dipole::hydrogenBondRangeMinor.getMax();
    Scalar c=Protein::Dipole::hydrogenBondRangeMajor.getMin();
    const int nb=6;
    Scalar x[nb],y[nb];
    x[0]=Scalar(0);
    y[0]=Protein::Dipole::hydrogenBondRangeMinor.getMin();
    Scalar tan2AngleMax=(Scalar(1)-Math::sqr(Protein::Dipole::cosAngleMin))/Math::sqr(Protein::Dipole::cosAngleMin);
    Scalar denom=Scalar(1)+tan2AngleMax;
    y[1]=(-a+Math::sqrt(Math::sqr(c)*denom-Math::sqr(a)*tan2AngleMax))/denom;
    x[1]=y[1]*Math::sqrt(tan2AngleMax);
    x[nb-1]=Scalar(0);
    y[nb-1]=b;
    Scalar angleMax=Math::acos(Protein::Dipole::cosAngleMin);
    for(int i=2;i<nb-1;++i)
    {
        Scalar angle=Scalar(nb-i)*angleMax/Scalar(nb-2);
        x[i]=b*Math::sin(angle);
        y[i]=b*Math::cos(angle);
    }
    const int na=8;
    Scalar sj[na],cj[na];
    for(int j=0;j<na;++j)
    {
        Scalar angle=Scalar(j)*Scalar(2)*Math::Constants<Scalar>::pi/Scalar(na);
        sj[j]=Math::sin(angle);
        cj[j]=Math::cos(angle);
    }

    /* Create a display list for small hydrogen cages: */
    glNewList(dataItem->hydrogenCageSmallDisplayListId,GL_COMPILE);
    for(int i=2;i<nb-1;++i)
    {
        glBegin(GL_LINE_LOOP);
        for(int j=0;j<na;++j)
            glVertex(x[i]*cj[j],y[i],x[i]*sj[j]);
        glEnd();
    }
    for(int j=0;j<na/2;++j)
    {
        glBegin(GL_LINE_STRIP);
        for(int i=2;i<nb;++i)
            glVertex(x[i]*cj[j],y[i],x[i]*sj[j]);
        for(int i=nb-2;i>=2;--i)
            glVertex(x[i]*cj[j+na/2],y[i],x[i]*sj[j+na/2]);
        glEnd();
    }
    glEndList();

    /* Create a display list for large hydrogen cages: */
    glNewList(dataItem->hydrogenCageLargeDisplayListId,GL_COMPILE);
    for(int i=1;i<nb-1;++i)
    {
        glBegin(GL_LINE_LOOP);
        for(int j=0;j<na;++j)
            glVertex(x[i]*cj[j],y[i],x[i]*sj[j]);
        glEnd();
    }
    for(int j=0;j<na/2;++j)
    {
        glBegin(GL_LINE_LOOP);
        for(int i=1;i<nb;++i)
            glVertex(x[i]*cj[j],y[i],x[i]*sj[j]);
        for(int i=nb-2;i>=0;--i)
            glVertex(x[i]*cj[j+na/2],y[i],x[i]*sj[j+na/2]);
        glEnd();
    }
    glEndList();

    /* Initialize the line illuminator: */
    bondIlluminator.init(contextData);
    bondIlluminator.setMaterial(contextData,bondMaterial);
}

void ProteinRenderer::updateStructureFlags(void)
{
    /* Synchronize structure parameter list with protein's secondary structure sequence: */
    structureFlags.clear();
    for(const Protein::SecondaryStructure* sPtr=protein->secondaryStructures;sPtr!=0;sPtr=sPtr->succ)
    {
        ConfigurationFile::SectionIterator structureSection=configFileSection;
        switch(sPtr->getStructureType())
        {
            case Protein::SecondaryStructure::COIL:
                structureSection.setSection("./Coil");
                break;
            
            case Protein::SecondaryStructure::ALPHA_HELIX:
                structureSection.setSection("./AlphaHelix");
                break;

            case Protein::SecondaryStructure::BETA_STRAND:
                structureSection.setSection("./BetaStrand");
                break;
        }
        structureFlags.push_back(StructureFlags(structureSection,sPtr));
    }
    resetAllBackboneColor();
    
    /* Create backbone ribbon spline: */
    createBackboneRibbonSpline();
    
    /* Create cartoon splines: */
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
        createCartoonSpline(*sfIt);
    
    /* Fill in protein data structures: */
    updateProtein();
}

void ProteinRenderer::updateProtein(void)
{
    /* Update the protein's bounding box: */
    Point min,max;
    min=max=protein->atoms->getPosition();
    for(const Protein::ChainAtom* aPtr=protein->atoms->succ;aPtr!=0;aPtr=aPtr->succ)
    {
        for(int i=0;i<3;++i)
        {
            if(min[i]>aPtr->getPosition()[i])
                min[i]=aPtr->getPosition()[i];
            else if(max[i]<aPtr->getPosition()[i])
                max[i]=aPtr->getPosition()[i];
        }
    }
    for(int i=0;i<3;++i)
    {
        min[i]-=maxAtomRadius*2.0;
        max[i]+=maxAtomRadius*2.0;
    }
    boundingBox=Geometry::Box<Scalar,3>(min,max);
    
    if(drawCartoon || drawTube)
    {
        /* Create and evaluate cartoon splines: */
        for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
            updateCartoonSpline(*sfIt);
    }
    
    if(drawCollisions)
    {
        /* Create the collision detection data structure: */
        double optCollisionCellSize=Math::pow(boundingBox.getSize(0)*boundingBox.getSize(1)*boundingBox.getSize(2)/double(protein->getNumAtoms()),1.0/3.0);
        double maxCollisionRadius=0.0;
        for(int i=0;i<118;++i)
        {
            double collisionRadius=Atom::getCovalentRadius((Atom::Element)i)*0.75;
            if(maxCollisionRadius<collisionRadius)
                maxCollisionRadius=collisionRadius;
        }
        if(optCollisionCellSize<2.0*maxCollisionRadius)
            optCollisionCellSize=2.0*maxCollisionRadius;
        for(int i=0;i<3;++i)
        {
            numCollisionCells[i]=int(Math::floor(boundingBox.getSize(i)/optCollisionCellSize));
            collisionCellSize[i]=boundingBox.getSize(i)/double(numCollisionCells[i]);
            numCollisionCells[i]+=2;
        }
        delete[] collisionCells;
        collisionCells=new CollisionCell[numCollisionCells[2]*numCollisionCells[1]*numCollisionCells[0]];
        collisionCellBase=&collisionCells[(1*numCollisionCells[1]+1)*numCollisionCells[0]+1];

        /* Insert all atoms into the collision detection structure: */
        AtomListItem* aliPtr=atomListItems;
        for(int i=0;i<protein->getNumAtoms();++i,++aliPtr)
        {
            int cellIndex[3];
            for(int i=0;i<3;++i)
                cellIndex[i]=int(Math::floor((aliPtr->atom->getPosition()[i]-boundingBox.getOrigin()[i])/collisionCellSize[i]));
            aliPtr->cell=&collisionCellBase[(cellIndex[2]*numCollisionCells[1]+cellIndex[1])*numCollisionCells[0]+cellIndex[0]];
            aliPtr->cellSucc=aliPtr->cell->atoms;
            aliPtr->cell->atoms=aliPtr;
        }
    }
    
    /* Copy all backbone atom positions into spline control point array: */
    int pointIndex=0;
    for(const Protein::Residue* rPtr=protein->residues;rPtr!=0;rPtr=rPtr->succ)
    {
        if(backboneRibbonUseAllAtoms)
        {
            for(std::vector<Protein::ChainAtom*>::const_iterator aIt=rPtr->backbone.begin();aIt!=rPtr->backbone.end();++aIt,++pointIndex)
                backboneRibbonSpline->setPoint(pointIndex,(*aIt)->getPosition());
        }
        else
        {
            backboneRibbonSpline->setPoint(pointIndex,rPtr->backbone[1]->getPosition());
            ++pointIndex;
        }
    }
}

void ProteinRenderer::setViewDirection(const GLLineIlluminator::Vector& viewDirection)
{
    bondIlluminator.disableAutoView();
    bondIlluminator.setViewDirection(viewDirection);
}

void ProteinRenderer::setLightDirection(const GLLineIlluminator::Vector& lightDirection)
{
    bondIlluminator.disableAutoLight();
    bondIlluminator.setLightDirection(lightDirection);
}

void ProteinRenderer::glRenderAction(GLContextData& contextData) const
{
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);

    /* Upload new bond material if necessary: */
    if(dataItem->bondMaterialVersion!=bondMaterialVersion)
    {
        bondIlluminator.setMaterial(contextData,bondMaterial);
        dataItem->bondMaterialVersion=bondMaterialVersion;
    }
    
    /* Render line parts of protein visualization: */
    glDisable(GL_LIGHTING);
    // if(drawBonds)
      //glDrawBonds(contextData);
    if(drawLine)
        glDrawLines(contextData);
    if(drawBackbone)
        glDrawBackbone(contextData);
    if(drawHydrogenBonds)
        glDrawHydrogenBonds(contextData);
    if(drawHydrogenBondSites)
        glDrawHydrogenBondSites(contextData);
    if(drawHydrogenCages)
        glDrawHydrogenCages(contextData);
    
    /* Render polygon parts of protein visualization: */
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    if(drawBonds)
        glDrawBondCylinders(contextData);
    if(drawCPK)
        glDrawCPK(contextData);
    if(drawTube)
        glDrawTube(contextData);
    if(drawAtoms)
        glDrawAtoms(contextData);
    if(drawBackboneRibbon)
        glDrawBackboneRibbon(contextData);
    if(drawCartoon)
        glDrawCartoon(contextData);
    if(drawCollisions)
        glDrawCollisions(contextData);
}

void ProteinRenderer::highlightResidue(GLContextData& contextData,const Protein::Residue* rPtr) const
{
    #if 0
    /* Draw residue's bonds as thicker unlit yellow sticks: */
    glDisable(GL_LIGHTING);
    glLineWidth(4.0f);
    glColor(1.0,1.0,0.0);

    /* Iterate through all atoms in this residue: */
    glBegin(GL_LINES);
    for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
    {
        /* Iterate through all bonds for this atom: */
        for(std::vector<Atom*>::const_iterator bIt=aPtr->getBonds().begin();bIt!=aPtr->getBonds().end();++bIt)
        {
            glVertex(aPtr->getPosition());
            glVertex(Geometry::mid(aPtr->getPosition(),(*bIt)->getPosition()));
        }
    }
    glEnd();
    #elif 0
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    
    /* Set atom rendering parameters: */
    glMaterial(GL_FRONT,atomMaterial);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT,GL_DIFFUSE);
    glPushMatrix();
    
    /* Iterate through all atoms in this residue: */
    Point currentOrigin=Point::origin;
    for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
    {
        Vector displacement=aPtr->getPosition()-currentOrigin;
        glTranslate(displacement);
        glCallList(dataItem->atomSphereDisplayListBaseId+aPtr->getType());
        currentOrigin=aPtr->getPosition();
    }
    
    /* Reset OpenGL state: */
    glPopMatrix();
    glDisable(GL_COLOR_MATERIAL);
    #else
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    
    /* Set atom rendering parameters: */
    glMaterial(GL_FRONT,bondMaterial);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT,GL_DIFFUSE);
    
    /* Iterate through all atoms in this residue: */
    for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
    {
        /* Iterate through all bonds for this atom: */
        for(std::vector<Atom*>::const_iterator bIt=aPtr->getBonds().begin();bIt!=aPtr->getBonds().end();++bIt)
        {
            const Protein::ChainAtom* oPtr=static_cast<const Protein::ChainAtom*>(*bIt);
            if(oPtr->getResidue()!=rPtr||aPtr<oPtr)
            {
                /* Draw a bicolor capped cylinder around the bond axis: */
                Vector bondAxis=oPtr->getPosition()-aPtr->getPosition();
                double bondAxisLen=Geometry::mag(bondAxis);
                bondAxis/=Scalar(bondAxisLen);

                for(int i=0;i<=numBondVertices;++i)
                {
                    dataItem->bondVertices[1][i][2]=GLfloat(bondAxisLen*0.5);
                    dataItem->bondVertices[2][i][2]=GLfloat(bondAxisLen);
                }

                /* Set OpenGL modelview matrix to align z axis with bond axis: */
                glPushMatrix();
                glTranslate(aPtr->getPosition()-Point::origin);
                Vector rotAxis(-bondAxis[1],bondAxis[0],Scalar(0));
                Scalar rotAngle=Math::deg(Math::acos(bondAxis[2]));
                glRotate(rotAngle,rotAxis);

                /* Draw root atom cap: */
                glBegin(GL_POLYGON);
                glColor(elementColors[aPtr->getType()]);
                glNormal(0.0,0.0,-1.0);
                for(int i=numBondVertices;i>=0;--i)
                    glVertex(dataItem->bondVertices[0][i]);
                glEnd();

                /* Draw cylinder at root atom: */
                glBegin(GL_QUAD_STRIP);
                for(int i=0;i<=numBondVertices;++i)
                {
                    glNormal(dataItem->bondNormals[i]);
                    glVertex(dataItem->bondVertices[1][i]);
                    glVertex(dataItem->bondVertices[0][i]);
                }
                glEnd();

                /* Draw cylinder at other atom: */
                glBegin(GL_QUAD_STRIP);
                glColor(elementColors[(*bIt)->getType()]);
                for(int i=0;i<=numBondVertices;++i)
                {
                    glNormal(dataItem->bondNormals[i]);
                    glVertex(dataItem->bondVertices[2][i]);
                    glVertex(dataItem->bondVertices[1][i]);
                }
                glEnd();

                /* Draw other atom cap: */
                glBegin(GL_POLYGON);
                glNormal(0.0,0.0,1.0);
                for(int i=0;i<=numBondVertices;++i)
                    glVertex(dataItem->bondVertices[2][i]);
                glEnd();

                /* Reset OpenGL modelview matrix: */
                glPopMatrix();
            }
        }
    }
    
    /* Reset OpenGL state: */
    glDisable(GL_COLOR_MATERIAL);
    #endif
	labelResidue(rPtr); 
	labelAtom(rPtr);
	showResidueAngles(rPtr);
}

void ProteinRenderer::labelResidue(const Protein::Residue* rPtr) const
{
    if(!drawResidueName || rPtr == NULL) return;

    Point p=rPtr->backbone[0]->getPosition();
    p=Geometry::mid(rPtr->backbone[0]->getPosition(),rPtr->backbone[2]->getPosition());
    p=Geometry::mid(p,rPtr->backbone[1]->getPosition());

	text->resetFontSize();
	text->setFontColor(1.0, 0.0, 0.0);
	text->setFontSize(36);

	glDisable (GL_LIGHTING);
	/* Set OpenGL modelview matrix to align z axis with bond axis: */
    glPushMatrix();
	glTranslate(p[0]+3.0, p[1]+3.0, p[2]);
	char *buf = new char[10];
	buf[0]='\0';

	if(rPtr->getModelId() > 0) 
		sprintf(buf,"%s(%d:%d)",rPtr->getPdbResidueName(),rPtr->getModelId(),rPtr->getPdbResidueIndex());
	else if(strncmp(rPtr->getChainId(),"",1))
		sprintf(buf,"%s(%s:%d)",rPtr->getPdbResidueName(),rPtr->getChainId(),rPtr->getPdbResidueIndex());
	else
		sprintf(buf,"%s(%d)",rPtr->getPdbResidueName(),rPtr->getPdbResidueIndex());

	text->drawFreeString(buf);
	glColor3f(1.0f, 0.0f, 0.0f);
	glPopMatrix();
   	glEnable(GL_LIGHTING);
	delete [] buf;
}

void ProteinRenderer::labelAtom(const Protein::Residue* rPtr) const
{
	if(!drawAtomNames || rPtr == NULL)	return;
	
    /* Get a pointer to the context entry: */
    Point p=rPtr->backbone[0]->getPosition();
	p=Geometry::mid(rPtr->backbone[0]->getPosition(),rPtr->backbone[2]->getPosition());
	p=Geometry::mid(p,rPtr->backbone[1]->getPosition());

	text->resetFontSize();
	text->setFontSize(14);
	text->setFontColor(1.0, 0.0, 0.0);
    glDisable(GL_LIGHTING);

    /* Iterate through all atoms in this residue: */
    for(const Protein::ChainAtom* aPtr=rPtr->beginAtom;aPtr!=rPtr->endAtom;aPtr=aPtr->succ)
    {
		/* Set OpenGL modelview matrix to align z axis with bond axis: */
		glPushMatrix();
		glTranslate(aPtr->getPosition()-Point::origin);
		text->setFontColor(1.0, 0.0, 0.0);
		text->drawFreeString(aPtr->getAtomName());
		glPopMatrix();
	}
	glEnable(GL_LIGHT0);
}

void ProteinRenderer::showResidueAngles(const Protein::Residue* rPtr) const
{
    if(!drawResidueAngles || rPtr == NULL) return;

    MD::Scalar* phis=new MD::Scalar[1];
    MD::Scalar* psis=new MD::Scalar[1];
    protein->getDihedralAngles (
	    protein->getResidueIndex(rPtr),1,phis,psis );

	/* Get a pointer to the context entry: */
    Point p=rPtr->backbone[0]->getPosition();
	p=Geometry::mid(rPtr->backbone[0]->getPosition(),rPtr->backbone[2]->getPosition());
	p=Geometry::mid(p,rPtr->backbone[1]->getPosition());

	text->setFontColor(0.0, 0.0, 1.0);
	text->resetFontSize();
	text->setFontSize(14);
    glDisable(GL_LIGHTING);

	sprintf(text->fontbuffer, "Phi: %f", Math::deg(phis[0]) );
	text->drawFlatString(text->fontbuffer, 10, 20);

	sprintf(text->fontbuffer, "Psi: %f", Math::deg(psis[0]) );
	text->drawFlatString(text->fontbuffer, 10, 45);

    glEnable(GL_LIGHT0);

	delete[] phis;
    delete[] psis;

}

void ProteinRenderer::updateEnergy(double value)
{
    if ( !isnan(value) )
        energyValue = value;
	sprintf(text->fontbuffer, "Energy: %g", energyValue);
	text->resetFontSize();
	text->setFontSize(14);
	text->setFontColor(1.0, 0.0, 0.0);
	text->drawFlatString(text->fontbuffer, renderWindow->w()-200, 10);
}

void ProteinRenderer::showRMSD(double value)
{
	if (value <= -1)
		sprintf(text->fontbuffer, "RMSD:");
	else
		sprintf(text->fontbuffer, "RMSD: %g", value);
	
	text->resetFontSize();
	text->setFontSize(14);
	text->setFontColor(1.0, 0.0, 0.0);
	text->drawFlatString(text->fontbuffer, 10, renderWindow->h()-20);
}

int ProteinRenderer::glPick(GLContextData& contextData,int which) const
{
    /* Get a pointer to the context entry: */
    DataItem* dataItem=contextData.retrieveDataItem<DataItem>(this);
    
    /* Upload new bond material if necessary: */
    if(dataItem->bondMaterialVersion!=bondMaterialVersion)
    {
        bondIlluminator.setMaterial(contextData,bondMaterial);
        dataItem->bondMaterialVersion=bondMaterialVersion;
    }
    
    /* Start OpenGL selection mode: */
    GLsizei hitBufferSize=structureFlags.size()*7*5; // Enough entries for every structure hitting on every rendering pass with 2 names on the stack
    GLuint* hitBuffer=new GLuint[hitBufferSize];
    glSelectBuffer(hitBufferSize,hitBuffer);
    glRenderMode(GL_SELECT);
    glInitNames();
    
    /* Render line parts of protein visualization: */
    glDisable(GL_LIGHTING);
    // if(drawBonds)
      //glDrawBonds(contextData);
    if(drawLine)
        glDrawLines(contextData);
    if(drawBackbone)
        glDrawBackbone(contextData);
    if(drawHydrogenBonds)
        glDrawHydrogenBonds(contextData);
    if(drawHydrogenBondSites)
        glDrawHydrogenBondSites(contextData);
    if(drawHydrogenCages)
        glDrawHydrogenCages(contextData);
    
    /* Render polygon parts of protein visualization: */
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    if(drawBonds)
        glDrawBondCylinders(contextData);
    if(drawCPK)
        glDrawCPK(contextData);
    if(drawTube)
        glDrawTube(contextData);
    if(drawAtoms)
        glDrawAtoms(contextData);
    if(drawBackboneRibbon)
        glDrawBackboneRibbon(contextData);
    if(drawCartoon)
        glDrawCartoon(contextData);
    if(drawCollisions)
        glDrawCollisions(contextData);
    
    /* Process hit records: */
    int numHits=glRenderMode(GL_RENDER);
    GLuint* hPtr=hitBuffer;
    GLuint minZ=0xffffffffU;
    int minStructure=-1;
    for(int i=0;i<numHits;++i)
    {
        GLuint numNames=hPtr[0];
        GLuint z1=hPtr[1];
        hPtr+=3;
        if(numNames>which&&z1<minZ)
        {
            minStructure=int(hPtr[which]);
            minZ=z1;
        }
        hPtr+=numNames;
    }
    delete[] hitBuffer;
    
    return minStructure;
}

ProteinRenderer::Color ProteinRenderer::getResidueBackboneColor(const Protein::Residue* rPtr) const
{
    /* Find structure flags for the residue's secondary structure: */
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==rPtr->secondaryStructure)
            return sfIt->backboneColor;
    }
    
    return Color(0,0,0);
}

bool ProteinRenderer::getDrawAtoms(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawAtoms;
    }
    return false;
}

void ProteinRenderer::setDrawAtoms(const Protein::StructureSelector& selector,bool newDrawAtoms)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->drawAtoms=newDrawAtoms;
    }
}

bool ProteinRenderer::getDrawBonds(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawBonds;
    }
    return false;
}

void ProteinRenderer::setDrawBonds(const Protein::StructureSelector& selector,bool newDrawBonds)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->drawBonds=newDrawBonds;
    }
}

bool ProteinRenderer::getDrawCPK(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawCPK;
    }
    return false;
}

void ProteinRenderer::setDrawCPK(const Protein::StructureSelector& selector,bool newDrawCPK)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->drawCPK=newDrawCPK;
    }
}

bool ProteinRenderer::getDrawTube(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawTube;
    }
    return false;
}

void ProteinRenderer::setDrawTube(const Protein::StructureSelector& selector,bool newDrawTube)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
        {
		    sfIt->drawTube=newDrawTube;
            if(drawTube&&sfIt->drawTube)
                updateCartoonSpline(*sfIt);
		}		

	}
}

bool ProteinRenderer::getDrawLine(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawLine;
    }
    return false;
}

void ProteinRenderer::setDrawLine(const Protein::StructureSelector& selector,bool newDrawLine)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->drawLine=newDrawLine;
    }
}

bool ProteinRenderer::getDrawBackbone(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawBackbone;
    }
    return false;
}

void ProteinRenderer::setDrawBackbone(const Protein::StructureSelector& selector,bool newDrawBackbone)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->drawBackbone=newDrawBackbone;
    }
}

bool ProteinRenderer::getDrawBackboneRibbon(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawBackboneRibbon;
    }
    return false;
}

void ProteinRenderer::setDrawBackboneRibbon(const Protein::StructureSelector& selector,bool newDrawBackboneRibbon)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->drawBackboneRibbon=newDrawBackboneRibbon;
    }
}

bool ProteinRenderer::getDrawCartoon(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawCartoon;
    }
    return false;
}

void ProteinRenderer::setDrawCartoon(const Protein::StructureSelector& selector,bool newDrawCartoon)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
        {
            sfIt->drawCartoon=newDrawCartoon;
            if(drawCartoon&&sfIt->drawCartoon)
                updateCartoonSpline(*sfIt);
        }
    }
}

void ProteinRenderer::resetAllBackboneColor(void)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        switch(sfIt->structure->structureType)
        {
            case Protein::SecondaryStructure::COIL:
                sfIt->backboneColor=coilColor;
                break;

            case Protein::SecondaryStructure::ALPHA_HELIX:
                sfIt->backboneColor=alphaHelixColor;
                break;

            case Protein::SecondaryStructure::BETA_STRAND:
                sfIt->backboneColor=betaStrandColor;
                break;
        }
    }
}

void ProteinRenderer::resetBackboneColor(const Protein::StructureSelector& selector)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
        {
            switch(sfIt->structure->structureType)
            {
                case Protein::SecondaryStructure::COIL:
                    sfIt->backboneColor=coilColor;
                    break;

                case Protein::SecondaryStructure::ALPHA_HELIX:
                    sfIt->backboneColor=alphaHelixColor;
                    break;

                case Protein::SecondaryStructure::BETA_STRAND:
                    sfIt->backboneColor=betaStrandColor;
                    break;
	       
	       default:
	    	    break;
	    }
        }
    }
}

void ProteinRenderer::setBackboneColor(const Protein::StructureSelector& selector,const Color& newBackboneColor)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->backboneColor=newBackboneColor;
    }
}

bool ProteinRenderer::getDrawHydrogenBondSites(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            return sfIt->drawHydrogenBondSites;
    }
    return false;
}

void ProteinRenderer::setAllDrawHydrogenBondSites(bool newDrawHydrogenBondSites)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
        sfIt->drawHydrogenBondSites=newDrawHydrogenBondSites;
}

void ProteinRenderer::setDrawHydrogenBondSites(const Protein::StructureSelector& selector,bool newDrawHydrogenBondSites)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->drawHydrogenBondSites=newDrawHydrogenBondSites;
    }
}

bool ProteinRenderer::getDrawHydrogenCages(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->drawHydrogenCages;
    }
    return false;
}

void ProteinRenderer::setDrawHydrogenCages(const Protein::StructureSelector& selector,bool newDrawHydrogenCages)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->drawHydrogenCages=newDrawHydrogenCages;
    }
}

bool ProteinRenderer::getDrawLargeHydrogenCages(const Protein::StructureSelector& selector) const
{
    for(std::vector<StructureFlags>::const_iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(sfIt->structure==selector.structure)
            return sfIt->hydrogenCageLarge;
    }
    return false;
}

void ProteinRenderer::setDrawLargeHydrogenCages(const Protein::StructureSelector& selector,bool newDrawLargeHydrogenCages)
{
    for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
    {
        if(selector.structure==sfIt->structure)
            sfIt->hydrogenCageLarge=newDrawLargeHydrogenCages;
    }
}

}
