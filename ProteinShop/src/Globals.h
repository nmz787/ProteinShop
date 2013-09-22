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
Globals - Global variables for interactive protein manipulation program.

***********************************************************************/


class ProteinState;


#ifndef GLOBALS_INCLUDED
#define GLOBALS_INCLUDED

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
#include <pthread.h>
#include <ConfigurationFile.h>
#include <Geometry/Box.h>
#include <GLTypes.h>
#include "DragBox.h"
#include "Protein.h"
#include "ProteinClient.h"
#include "ProteinInteractor.h"
#include "ProteinRenderer.h"
#include "EnergyRenderer.h"
#include "EnergyAPI.h"
#include "UndoBuffer.h"
#include "Message.h"

/* Forward declarations: */
class GLContextData;

/*****************************************************
Class to decode color values from configuration files:
*****************************************************/

template <>
class ValueCoder<GLColor<GLfloat,4> >
{
    /* Methods: */
    public:
    static std::string encode(const GLColor<GLfloat,4>& value)
    {
        char buffer[80];
        if(value[3]==1.0f)
            sprintf(buffer,"(%f, %f, %f)",value[0],value[1],value[2]);
        else
            sprintf(buffer,"(%f, %f, %f, %f)",value[0],value[1],value[2],value[3]);
        return std::string(buffer);
    };
    static GLColor<GLfloat,4> decode(const std::string& text)
    {
        GLColor<GLfloat,4> result;
        int numComponents=sscanf(text.c_str(),"(%f, %f, %f, %f)",&result[0],&result[1],&result[2],&result[3]);
        if(numComponents==3)
            result[3]=1.0f;
        else if(numComponents!=4)
            throw DecodingError(std::string("Unable to convert \"")+text+std::string("\" to GLColor<GLfloat,4>"));
        return result;
    };
};

template <>
class ValueCoder<GLMaterial>
{
    /* Methods: */
    public:
    static std::string encode(const GLMaterial& value)
    {
        char buffer[256];
        char* bufPtr=buffer;
        bufPtr+=sprintf(bufPtr,"(%f, %f, %f), ",value.getAmbient()[0],value.getAmbient()[1],value.getAmbient()[2]);
        bufPtr+=sprintf(bufPtr,"(%f, %f, %f), ",value.getDiffuse()[0],value.getDiffuse()[1],value.getDiffuse()[2]);
        bufPtr+=sprintf(bufPtr,"(%f, %f, %f), ",value.getSpecular()[0],value.getSpecular()[1],value.getSpecular()[2]);
        bufPtr+=sprintf(bufPtr,"%f, ",value.getShininess());
        return std::string(buffer);
    };
    static GLMaterial decode(const std::string& text)
    {
        GLMaterial::Color ambient,diffuse,specular,emission;
        GLfloat shininess;
        if(sscanf(text.c_str(),"(%f, %f, %f), (%f, %f, %f), (%f, %f, %f), %f",&ambient[0],&ambient[1],&ambient[2],&diffuse[0],&diffuse[1],&diffuse[2],&specular[0],&specular[1],&specular[2],&shininess)!=10)
            throw DecodingError(std::string("Unable to convert \"")+text+std::string("\" to GLMaterial"));
        return GLMaterial(ambient,diffuse,specular,shininess);
    };
};

// proteinMutex guards or regulates the variables in this block except ikUpdatePosted
extern pthread_mutex_t proteinMutex;
extern pthread_cond_t ikUpdateRequestedCond;
extern volatile bool ikUpdateRequested;
extern pthread_cond_t ikUpdateDoneCond;
extern volatile bool ikUpdateDone;
extern volatile bool ikUpdatePosted;
extern pthread_t ikUpdateThreadId;

extern ConfigurationFile* configFile; // Configuration file containing program parameters
extern UndoBuffer undoBuffer;
extern EnergyLibrary* energyLibrary; // Object representing a dynamically linked energy calculation library
extern PShopMessage msg;

extern void updateDihedralAngles(void);
extern void recalculateEnergy(void);
extern void updateAtomEnergies(void);
extern double getEnergyOutput(void);
extern void updateGui(void);
extern int zipAntiParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float plane[3], float center[3]);
extern int zipParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float plane[3], float center[3]);
extern int altZipAntiParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float plane[3], float center[3]);
extern int altZipParallel(MD::Protein::Residue* manipulatedResidue,MD::Protein::Residue* anchorResidue, float plane[3], float center[3]);
extern void redrawRenderWindows(void);
extern void updateProteinNow(void);


struct ProteinState
{
    MD::Protein* protein; // The visualized protein
    char *name; // name to use when displaying a list of proteins

    ProteinClient* client; // A protein client talking to a global optimization server
    uint proteinId; // Unique ID of current protein in server's optimization tree

    MD::Protein::Residue* selectedResidue; // The selected residue
    ProteinInteractor* interactor; // Interaction object

    MD::ProteinRenderer* proteinRenderer; // Renderer attached to the protein
    EnergyRenderer *energyRenderer; // energy renderer attached to the protein

    EnergyCalculator* energyCalculator; // Energy calculator attached to the protein
    bool visualizeEnergy; // Flag to toggle energy visualization
    double visualizeEnergyMinRange,visualizeEnergyMaxRange, engUpdateRate; // Value range of energy visualization color map

    ProteinState (const char *filename);
    ~ProteinState();
};

extern void checkDrawToggles(ProteinState* state);
extern void checkShowToggles(ProteinState* state);

// ProteinState managment
extern ProteinState *createProtein (const char *filename);
extern ProteinState *curProtein();
extern void deleteAllProteins();
extern bool deleteProtein (ProteinState *state);
extern ProteinState *getProtein (uint index);
extern uint numProteins();
extern bool setCurProtein (ProteinState *state);

/// Translate the errno variable into a string.
extern const char *errnoString();

/// Get the next larger power of two.
extern uint ceilingPowerOfTwo (uint number);

/// Clamp a scalar variable to an interval.
template <typename Scalar>
void clamp (Scalar &toClamp, const Scalar &lowerBound, const Scalar &upperBound)
{
    if ( toClamp < lowerBound ) toClamp = lowerBound;
    if ( toClamp > upperBound ) toClamp = upperBound;
}


#endif
