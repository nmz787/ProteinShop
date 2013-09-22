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


namespace MD { class ProteinRenderer; }


#ifndef PROTEINRENDERER_INCLUDED
#define PROTEINRENDERER_INCLUDED

#define USE_BONDCYLINDER_VERTEXARRAY 0

#include <utility>
#include <ConfigurationFile.h>
#include <Geometry/Box.h>
#include <Geometry/SplineCurve.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <GLVertex.h>
#include <GLColorMap.h>
#include <GLMaterial.h>
#include <GLContextData.h>
#include <GLLineIlluminator.h>
#include <GLText.h>
#include <ctime>

#include "Protein.h"

namespace MD {

class ProteinRenderer
{
    /* Embedded classes: */
    public:
    typedef GLMaterial::Color Color; // Data type for colors
    
    private:
    struct CollisionCell;
    
    struct AtomListItem // Structure associating atoms with collision detection cells
    {
        /* Elements: */
        public:
        const Protein::ChainAtom* atom; // Pointer to atom
        CollisionCell* cell; // Pointer to the cell containing the atom
        AtomListItem* cellSucc; // Pointer to next atom in the same cell
    };
    
    struct CollisionCell // Structure for collision detection cells
    {
        /* Elements: */
        public:
        AtomListItem* atoms; // Pointer to list of atoms inside this cell
        
        /* Constructors and destructors: */
        CollisionCell(void) // Creates empty cell
            :atoms(0)
        {
        };
    };
    
    struct StructureFlags // Structure holding rendering flags for each secondary structure
    {
        /* Elements: */
        public:
        const Protein::SecondaryStructure* structure; // The secondary structure
        bool drawAtoms;
        bool drawBonds;
        bool drawCPK;
        bool drawTube;
        bool drawLine;
        GLfloat bondWidth; // Cosmetic bond (line) width
        Color bondColor;
        bool drawBackbone;
        GLfloat backboneWidth;
        Color backboneColor;
        bool drawBackboneRibbon;
        Scalar backboneRibbonWidth;
        int backboneRibbonNumSamples;
        Scalar backboneRibbonParameterMin,backboneRibbonParameterMax;
        bool drawCartoon;
        Geometry::SplineCurve<Scalar,3>* cartoonSpline1; // Pointer to major spline used by cartoon renderer
        Geometry::SplineCurve<Scalar,3>* cartoonSpline2; // Pointer to minor spline used by beta strand cartoon renderer
        int numCartoonSplineSamples;
        Scalar maxCartoonSplineParameter;
        Point* cartoonP;
        Vector* cartoonX;
        Vector* cartoonY;
        Vector* cartoonZ;
        bool drawHydrogenBondSites;
        GLfloat hydrogenBondSiteDiameter;
        GLfloat hydrogenBondSiteWidth;
        Color amideColor,carboxylColor;
        bool drawHydrogenCages;
        GLfloat hydrogenCageWidth;
        Color hydrogenCageColor;
        bool hydrogenCageLarge;
        
        /* Constructors and destructors: */
        StructureFlags(const ConfigurationFile::SectionIterator& configFileSection,const Protein::SecondaryStructure* sStructure); // Creates default rendering flags
        ~StructureFlags(void);
    };
	
    struct DataItem:public GLContextData::DataItem
    {
        /* Elements: */
        public:
        GLuint collisionSphereDisplayListId; // Display list for collision visualization sphere
        GLuint atomSphereDisplayListBaseId; // Base display list for atom spheres
		GLuint cpkSphereDisplayListId; // Display list for cpk sphere
        unsigned int bondMaterialVersion; // Version of bond material currently uploaded into the context
        #if USE_BONDCYLINDER_VERTEXARRAY
        GLVertex<void,0,void,0,GLfloat,GLfloat,3>* bondVertices; // Vertex array to render stick model cylinders
        GLuint* bondIndices[4]; // Four index arrays to render stick model cylinders
        #else
        GLVector<GLfloat,3>* bondNormals; // Ring of normals around fancy stick model cylinders
        GLVector<GLfloat,3>* bondVertices[3]; // Three rings of vertices around fancy stick model cylinders
        GLVector<GLfloat,3>* lineNormals; // Ring of normals around lines
        GLVector<GLfloat,3>* lineVertices[3]; // Three rings of vertices lines
        #endif
        GLuint hydrogenCageSmallDisplayListId; // Display list for small hydrogen cages
        GLuint hydrogenCageLargeDisplayListId; // Display list for large hydrogen cages
        
        /* Constructors and destructors: */
        DataItem(int numBondVertices, int numLineVertices)
            :collisionSphereDisplayListId(glGenLists(1)),
             atomSphereDisplayListBaseId(glGenLists(118)),
             cpkSphereDisplayListId(glGenLists(1)),
			 bondMaterialVersion(0),
             #if USE_BONDCYLINDER_VERTEXARRAY
             bondVertices(new GLVertex<void,0,void,0,GLfloat,GLfloat,3>[numBondVertices*3]),
             bondIndices(0),
             #else
             bondNormals(new GLVector<GLfloat,3>[numBondVertices+1]),
             lineNormals(new GLVector<GLfloat,3>[numLineVertices+1]),
 			 #endif
             hydrogenCageSmallDisplayListId(glGenLists(2)),
             hydrogenCageLargeDisplayListId(hydrogenCageSmallDisplayListId+1)
        {
            #if USE_BONDCYLINDER_VERTEXARRAY
            bondIndices[0]=new GLuint[numBondVertices];
            bondIndices[1]=new GLuint[(numBondVertices+1)*2];
            bondIndices[2]=new GLuint[(numBondVertices+1)*2];
            bondIndices[3]=new GLuint[numBondVertices];
            #else
            for(int ring=0;ring<3;++ring)
			{
                bondVertices[ring]=new GLVector<GLfloat,3>[numBondVertices+1];
                lineVertices[ring]=new GLVector<GLfloat,3>[numLineVertices+1];
			}
			#endif
        };
        ~DataItem(void)
        {
            glDeleteLists(collisionSphereDisplayListId,1);
            glDeleteLists(atomSphereDisplayListBaseId,118);
            glDeleteLists(cpkSphereDisplayListId,1);
            #if USE_BONDCYLINDER_VERTEXARRAY
            delete[] bondVertices;
            for(int part=0;part<4;++part)
                delete[] bondIndices[part];
            #else
            delete[] bondNormals;
            delete[] lineNormals;
            for(int ring=0;ring<3;++ring)
			{
                delete[] bondVertices[ring];
                delete[] lineVertices[ring];
            }
			#endif
            glDeleteLists(hydrogenCageSmallDisplayListId,2);
        };
    };

    /* Elements: */
    ConfigurationFile::SectionIterator configFileSection; // Configuration file section containing rendering parameters
    const Protein* protein; // Pointer to the visualized protein
    AtomListItem* atomListItems; // C-style array of atom list items; stays valid throughout lifetime of renderer
    double collisionCellSize[3];
    int numCollisionCells[3],collisionCellsWidth[3];
    CollisionCell* collisionCells;
    CollisionCell* collisionCellBase;
    static const double maxAtomRadius; // Maximum radius of any atom sphere
    Geometry::Box<Scalar,3> boundingBox; // Protein's current bounding box
    std::vector<StructureFlags> structureFlags; // Rendering flags for each secondary structure
    bool drawAtoms; // Global flag for atom rendering
    GLMaterial atomMaterial; // Material for atom rendering
    int atomTesselation; // Number of subdivisions for atom spheres
    static const Color elementColors[118]; // Standard colors to render atoms
    bool mapAtomValues; // Flag to enable value mapping for atoms
    GLColorMap atomColorMap; // Color map for mapping atom values
    float* atomValues; // Array of values that can be used to colormap atoms
    bool drawBonds; // Global flag for bond rendering
    GLLineIlluminator bondIlluminator;
    GLMaterial bondMaterial;
	GLText*  text;
    double energyValue;
    unsigned int bondMaterialVersion;
    int numBondVertices; // Number of vertices for fancy stick model rendering
    GLfloat bondRadius; // Cylinder radius for fancy stick model rendering
	GLfloat lineRadius; // Line radius for CPK model rendering
	GLfloat cpkSphereRadius; // Sphere radius for CPK model rendering
	int numLineVertices; // Number of vertices for CPK model rendering
    bool drawCPK; // Global flag for CPK rendering
    bool drawTube; // Global flag for CPK rendering
    bool drawLine; // Global flag for CPK rendering
    Color alphaHelixColor; // Color for alpha helices
    Color betaStrandColor; // Color for beta strands
    Color coilColor; // Color for coils
    Color highlightColor; // Color for selected structures
    Color phobicColor;		// Color for hydrophobic
	Color philicColor;		// Color for hydrophilic
	Color disulfideColor;	// Color for disulfide
    bool drawBackbone; // Global flag for backbone rendering
    bool drawBackboneRibbon; // Global flag for backbone ribbon rendering
    bool showHydrophobic; // Global flag for Hydrophobic rendering
	bool showHydrophilic; // Global flag for Hydrophilic rdering
	bool showDisulfide; // Global flag for disulfide bond
	bool drawAtomNames; // flag for atom name rendering
    bool drawResidueName; // flag for residue name rendering
    bool drawResidueAngles; // flag for residue angles rendering
	bool updateEnergyValue; // flag for energy value
    bool settime;
	time_t oldtime;
    GLMaterial backboneRibbonMaterial;
    bool backboneRibbonUseAllAtoms; // Flag if backbone spline uses all backbone atoms or only one per residue as control points
    int backboneRibbonDegree; // B-spline degree of backbone ribbon
    int backboneRibbonSampleDensity; // Sampling density for rendering the backbone ribbon
    Geometry::SplineCurve<Scalar,3>* backboneRibbonSpline; // Spline curve used to render the backbone
    bool drawCartoon; // Global flag for cartoon rendering
	bool grayoutCartoon; // Global flag for cartoon grayout when align two proteins
    Color grayOutColor;  // cartoon grayout color
    GLMaterial cartoonMaterial;
    int cartoonDegree; // B-spline degree of cartoon splines
    int cartoonSampleDensity; // Sampling density for rendering cartoon splines
    Scalar alphaHelixWidth,alphaHelixThickness; // Width and thickness of alpha helix strip
    Scalar betaStrandWidth,betaStrandThickness; // Width and thickness of beta strand arrow
    Scalar betaStrandHeadWidth; // Scale factor for width of beta strand arrow heads
    int numCoilVertices; // Number of vertices for coil tubes
    Scalar coilRadius; // Radius of coil tubes
    bool drawHydrogenBonds; // Global flag for hydrogen bond rendering
    GLfloat hydrogenBondWidth;
    Color hydrogenBondColor;
    bool drawHydrogenBondSites; // Global flag for hydrogen bond site rendering
    bool drawHydrogenCages; // Global flag for hydrogen cage rendering
    bool drawCollisions; // Global flag for collision sphere rendering
    GLMaterial collisionSphereMaterial; // Material for collision spheres
    int collisionSphereTesselation; // Number of subdivisions for collision spheres

    /* Private methods: */
    void createBackboneRibbonSpline(void);
    void createCartoonSpline(StructureFlags& sf);
    void updateCartoonSpline(StructureFlags& sf);
    void glDrawAtoms(GLContextData& contextData) const;
    void glDrawBondCylinders(GLContextData& contextData) const;
    void glDrawBackbone(GLContextData& contextData) const;
    void glDrawBackboneRibbon(GLContextData& contextData) const;
    void glDrawCartoon(GLContextData& contextData) const;
    void glDrawHydrogenBonds(GLContextData& contextData) const;
    void glDrawHydrogenBondSites(GLContextData& contextData) const;
    void glDrawHydrogenCages(GLContextData& contextData) const;
    void glDrawCollisions(GLContextData& contextData) const;
	void glDrawCPK(GLContextData& contextData) const;
	void glDrawTube(GLContextData& contextData) const;
	void glDrawLines(GLContextData& contextData) const;

    /* Constructors and destructors: */
    public:
    ProteinRenderer(const ConfigurationFile::SectionIterator& sConfigFileSection,const Protein* sProtein); // Creates a renderer for the given protein
    ~ProteinRenderer(void);

    /* Methods: */
    void clearContext(GLContextData& contextData); // Delete per-context renderer data
    void initContext(GLContextData& contextData); // Initializes per-context renderer data
    bool isContextInitialized(GLContextData &contextData) const;
    void updateStructureFlags(void);  // Updates structure rendering flags after change to protein's secondary structure sequence
    void updateProtein(void); // Updates internal representation of the protein structure after changes
    template <class TransformationParam>
    std::pair<double,double> calcDepthRange(const TransformationParam& modelView) const // Returns depth value range for given modelview transformation
    {
        double frontDist,backDist;
        frontDist=backDist=-modelView.transform(boundingBox.getVertex(0))[2];
        for(int i=1;i<8;++i)
        {
            double dist=-modelView.transform(boundingBox.getVertex(i))[2];
            if(frontDist>dist)
                frontDist=dist;
            else if(backDist<dist)
                backDist=dist;
        }
        return std::pair<double,double>(frontDist-maxAtomRadius,backDist+maxAtomRadius);
    };
    void setViewDirection(const GLLineIlluminator::Vector& viewDirection); // Sets the view direction for rendering
    void setLightDirection(const GLLineIlluminator::Vector& lightDirection); // Sets the light direction for rendering
    void glRenderAction(GLContextData& contextData) const; // Renders the protein
    void highlightResidue(GLContextData& contextData,const Protein::Residue* rPtr) const; // Highlights a single residue
    void labelResidue(const Protein::Residue* rPtr) const; // label a single residue
	void labelAtom(const Protein::Residue* rPtr) const;
	void showResidueAngles(const Protein::Residue* rPtr) const;
	void updateEnergy(double value); 
	void showRMSD(double value); 
    int glPick(GLContextData& contextData,int which) const; // Uses OpenGL picking to return the index of the secondary structure (which==0) or residue (which==1) hit by the query point; -1 if none found
    Color getResidueBackboneColor(const Protein::Residue* rPtr) const; // Returns the color a residue's ba

    /* Parameter access methods: */
    bool getDrawAtoms(void) const
    {
        return drawAtoms;
    };
    void setDrawAtoms(bool newDrawAtoms)
    {
        drawAtoms=newDrawAtoms;
    };
    bool getMapAtomValues(void) const
    {
        return mapAtomValues;
    };
    void setMapAtomValues(bool newMapAtomValues)
    {
        mapAtomValues=newMapAtomValues;
    };
    void setMapAtomValueRange(float newMinValue,float newMaxValue) // Sets mapping range for atom values
    {
        atomColorMap.setScalarRange(newMinValue,newMaxValue);
    };
    void setAtomValue(int atomIndex,float value) // Sets colormapping source value for an atom
    {
        atomValues[atomIndex]=value;
    };
    bool getDrawBonds(void) const
    {
        return drawBonds;
    };
    void setDrawBonds(bool newDrawBonds)
    {
        drawBonds=newDrawBonds;
    };
    bool getDrawCPK(void) const
    {
        return drawCPK;
    };
    void setDrawCPK(bool newDrawCPK)
    {
        drawCPK=newDrawCPK;
    };
    bool getDrawTube(void) const
    {
        return drawTube;
    };
    void setDrawTube(bool newDrawTube)
    {
        drawTube=newDrawTube;
        if(drawTube)
        {
            /* Create and evaluate cartoon splines: */
            for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
                updateCartoonSpline(*sfIt);
        }
    };
    bool getDrawLine(void) const
    {
        return drawLine;
    };
    void setDrawLine(bool newDrawLine)
    {
        drawLine=newDrawLine;
    };
    bool getDrawBackbone(void) const
    {
        return drawBackbone;
    };
    void setDrawBackbone(bool newDrawBackbone)
    {
        drawBackbone=newDrawBackbone;
    };
    bool getDrawBackboneRibbon(void) const
    {
        return drawBackboneRibbon;
    };
    void setDrawBackboneRibbon(bool newDrawBackboneRibbon)
    {
        drawBackboneRibbon=newDrawBackboneRibbon;
    };
    bool getDrawCartoon(void) const
    {
        return drawCartoon;
    };
    void setDrawCartoon(bool newDrawCartoon)
    {
        drawCartoon=newDrawCartoon;
        if(drawCartoon)
        {
            /* Create and evaluate cartoon splines: */
            for(std::vector<StructureFlags>::iterator sfIt=structureFlags.begin();sfIt!=structureFlags.end();++sfIt)
                updateCartoonSpline(*sfIt);
        }
    };
    bool getDrawHydrogenBonds(void) const
    {
        return drawHydrogenBonds;
    };
    void setDrawHydrogenBonds(bool newDrawHydrogenBonds)
    {
        drawHydrogenBonds=newDrawHydrogenBonds;
    };
    bool getDrawHydrogenBondSites(void) const
    {
        return drawHydrogenBondSites;
    };
    void setDrawHydrogenBondSites(bool newDrawHydrogenBondSites)
    {
        drawHydrogenBondSites=newDrawHydrogenBondSites;
    };
    bool getDrawHydrogenCages(void) const
    {
        return drawHydrogenCages;
    };
    void setDrawHydrogenCages(bool newDrawHydrogenCages)
    {
        drawHydrogenCages=newDrawHydrogenCages;
    };
    bool getDrawCollisions(void) const
    {
        return drawCollisions;
    };
    void setDrawCollisions(bool newDrawCollisions)
    {
        drawCollisions=newDrawCollisions;
        if(drawCollisions)
            updateProtein(); // Might have to re-calculate collision data structure
    };

    /* Per-structure parameter access methods: */
    bool getDrawAtoms(const Protein::StructureSelector& selector) const;
    void setDrawAtoms(const Protein::StructureSelector& selector,bool newDrawAtoms);
    bool getDrawBonds(const Protein::StructureSelector& selector) const;
    void setDrawBonds(const Protein::StructureSelector& selector,bool newDrawBonds);
    bool getDrawCPK(const Protein::StructureSelector& selector) const;
    void setDrawCPK(const Protein::StructureSelector& selector,bool newDrawCPK);
    bool getDrawTube(const Protein::StructureSelector& selector) const;
    void setDrawTube(const Protein::StructureSelector& selector,bool newDrawTube);
    bool getDrawLine(const Protein::StructureSelector& selector) const;
    void setDrawLine(const Protein::StructureSelector& selector,bool newDrawLine);
    bool getDrawBackbone(const Protein::StructureSelector& selector) const;
    void setDrawBackbone(const Protein::StructureSelector& selector,bool newDrawBackbone);
    bool getDrawBackboneRibbon(const Protein::StructureSelector& selector) const;
    void setDrawBackboneRibbon(const Protein::StructureSelector& selector,bool newDrawBackboneRibbon);
    bool getDrawCartoon(const Protein::StructureSelector& selector) const;
    void setDrawCartoon(const Protein::StructureSelector& selector,bool newDrawCartoon);
    void grayOutCartoon(bool value) { grayoutCartoon=value; }
    void setShowHydrophobic(bool value) { showHydrophobic=value; }
    void setShowHydrophilic(bool value) { showHydrophilic=value; }
    void setShowDisulfide(bool value) { showDisulfide=value; }
    void setDrawAtomNames(bool value)		{ drawAtomNames=value; }
    void setDrawResidueName(bool value)		{ drawResidueName=value; }
    void setDrawResidueAngles(bool value)	{ drawResidueAngles=value; }
	void setEnergyUpdate(bool value)		{ updateEnergyValue=value; }
    bool getShowHydrophobic(void) 			{ return showHydrophobic; }
    bool getShowHydrophilic(void) 			{ return showHydrophilic; }
    bool getShowDisulfide(void) 			{ return showDisulfide; }
    bool getDrawAtomNames(void)				{ return drawAtomNames; }
    bool getDrawResidueName(void)			{ return drawResidueName; }
    bool getDrawResidueAngles(void)			{ return drawResidueAngles; }
	bool getEnergyUpdate(void)	   			{ return updateEnergyValue; }
    bool getGrayOutCartoon(void) 			{ return grayoutCartoon; }
    void resetAllBackboneColor(void);
    void resetBackboneColor(const Protein::StructureSelector& selector);
    void setBackboneColor(const Protein::StructureSelector& selector,const Color& newBackboneColor);
    void selectStructure(const Protein::StructureSelector& selector)
    {
        setBackboneColor(selector,highlightColor);
    };
    void deselectStructure(const Protein::StructureSelector& selector)
    {
        resetBackboneColor(selector);
    };
    bool getDrawHydrogenBondSites(const Protein::StructureSelector& selector) const;
    void setAllDrawHydrogenBondSites(bool newDrawHydrogenBondSites);
    void setDrawHydrogenBondSites(const Protein::StructureSelector& selector,bool newDrawHydrogenBondSites);
    bool getDrawHydrogenCages(const Protein::StructureSelector& selector) const;
    void setDrawHydrogenCages(const Protein::StructureSelector& selector,bool newDrawHydrogenCages);
    bool getDrawLargeHydrogenCages(const Protein::StructureSelector& selector) const;
    void setDrawLargeHydrogenCages(const Protein::StructureSelector& selector,bool newDrawLargeHydrogenCages);
};

}

#endif
