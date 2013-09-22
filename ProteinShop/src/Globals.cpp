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

#include <cerrno>
#include <cstdio>
#include <vector>

#include "Globals.h"
#include "ProteinGui.h"


/****************
Global variables:
****************/

ConfigurationFile* configFile=0;
UndoBuffer undoBuffer;
EnergyLibrary* energyLibrary=0;

pthread_mutex_t proteinMutex;
pthread_cond_t ikUpdateRequestedCond;
volatile bool ikUpdateRequested;
pthread_cond_t ikUpdateDoneCond;
volatile bool ikUpdateDone;
volatile bool ikUpdatePosted;
pthread_t ikUpdateThreadId;
PShopMessage msg;


/////////////////////////////
// global state management //
/////////////////////////////


using namespace std;


static vector<ProteinState*> s_proteins;
static uint s_curProtein = 0;


ProteinState::ProteinState (const char *filename) :
    protein (0),
    name (0),
    client (0),
    proteinId (0),
    selectedResidue (0),
    interactor (0),
    proteinRenderer (0),
    energyRenderer (0),
    energyCalculator (0),
    visualizeEnergy (false),
    visualizeEnergyMinRange (0.0),
    visualizeEnergyMaxRange (0.0),
    engUpdateRate(0.0)
{
    const char *basename = strrchr (filename, '/');
    if ( basename )
        ++basename;
    else
        basename = filename;
    char const *ext = strrchr (basename, '.');
    size_t len = strlen (basename);
    name = new char[len + 10];
    if ( !ext )
        strcpy (name, basename);
    else
    {
        len = (size_t) (ext - basename);
        strncpy (name, basename, len);
        name[len] = 0;
    }
    uint numDups = 0;
    for ( uint i = 0; i < s_proteins.size(); ++i )
    {
        if ( s_proteins[i]->proteinId >= proteinId )
            proteinId = (s_proteins[i]->proteinId + 1) * 2;
        if ( !strncmp(s_proteins[i]->name, name, len) )
            ++numDups;
    }
    if ( numDups )
    {
        char buf[9];
        snprintf (buf, 9, " <%d>", numDups + 1);
        strcat (name, buf);
    }
}


ProteinState::~ProteinState()
{
    using namespace MD;

    // free any OpenGL resources associated with this protein
    renderWindow->clearContext (this);

    // order of destruction should observe dependencies
    undoBuffer.deleteProtein (protein);
    delete name;
    delete client;
    delete interactor;
    delete proteinRenderer;
    delete energyRenderer;
    delete energyCalculator;
    delete protein;
}


ProteinState *createProtein (const char *filename)
{
    s_proteins.push_back (new ProteinState(filename));
    s_curProtein = s_proteins.size() - 1;
    return s_proteins.back();
}


ProteinState *curProtein()
{
    if ( s_curProtein < s_proteins.size() )
        return s_proteins[s_curProtein];
    else
        return 0;
}


void deleteAllProteins()
{
    for ( uint i = 0; i < s_proteins.size(); ++i )
    {
        // do not allow access back to state object during destruction
        ProteinState *state = s_proteins[i];
        s_proteins[i] = 0;
        delete state;
    }
    s_proteins.clear();
}


bool deleteProtein (ProteinState *state)
{
    if ( !state ) return false;
    for ( uint i = 0; i < s_proteins.size(); ++i )
    {
        if ( state == s_proteins[i] )
        {
            // do not allow access back to state object during destruction
            s_proteins[i] = 0;
            delete state;
            s_proteins.erase (s_proteins.begin() + i);
            if ( s_curProtein && s_curProtein == s_proteins.size() )
                --s_curProtein;
            return true;
        }
    }
    return false;
}


ProteinState *getProtein (uint index)
{
    if ( index < s_proteins.size() )
        return s_proteins[index];
    else
        return 0;
}


uint numProteins()
{
    return s_proteins.size();
}


bool setCurProtein (ProteinState *state)
{
    if ( !state ) return false;
    for ( uint i = 0; i < s_proteins.size(); ++i )
    {
        if ( s_proteins[i] == state )
        {
            s_curProtein = i;
            return true;
        }
    }
    return false;
}


const char *errnoString()
{
    const char *errorCode = 0;
    switch ( errno )
    {
        case EACCES:        errorCode = "EACCES";       break;
        case EFAULT:        errorCode = "EFAULT";       break;
        case EINVAL:        errorCode = "EINVAL";       break;
        case EPERM:         errorCode = "EPERM";        break;
        case ENAMETOOLONG:  errorCode = "ENAMETOOLONG"; break;
        case ENOENT:        errorCode = "ENOENT";       break;
        case ENOTDIR:       errorCode = "ENOTDIR";      break;
        case ENOMEM:        errorCode = "ENOMEM";       break;
        case ERANGE:        errorCode = "ERANGE";       break;
        case EROFS:         errorCode = "EROFS";        break;
        default: break;
    }
    return errorCode;
}


uint ceilingPowerOfTwo (uint number)
{
    uint result = 1;
    while ( result < number )
        result = result << 1;
    return result;
}

