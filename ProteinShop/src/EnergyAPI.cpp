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
EnergyAPI - Classes to dynamically load energy computation libraries at
run-time.
***********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dlfcn.h>
#include <string>
#include <stdexcept>

#include "EnergyAPI.h"


/******************************
Methods of class EnergyLibrary:
******************************/

EnergyLibrary::EnergyLibrary (void *sDsoHandle) :
    dsoHandle (sDsoHandle)
{
}

EnergyLibrary* EnergyLibrary::load(const char* dsoName)
{
    /* Open DSO containing class: */
    void* dsoHandle=dlopen(dsoName,RTLD_LAZY|RTLD_GLOBAL);
    if(dsoHandle==0) // Check for errors during DSO loading
        throw std::runtime_error(std::string("EnergyLibrary::load: Unable to open energy calculation DSO \"")+std::string(dsoName)+std::string("\" due to ")+std::string(dlerror()));

    /* Get address to library object creation function: */
    CreateEnergyLibraryFunction createEnergyLibrary=(CreateEnergyLibraryFunction)dlsym(dsoHandle,"createEnergyLibrary");
    if(createEnergyLibrary==0) // Check for errors during function lookup
    {
        dlclose(dsoHandle);
        throw std::runtime_error(std::string("EnergyLibrary::load: Unable to retrieve energy calculator creation function from DSO due to ")+std::string(dlerror()));
    }

    /* Create and return library object: */
    return createEnergyLibrary(dsoHandle);
}

EnergyLibrary::~EnergyLibrary(void)
{
    // this call causes a segfault for some unknown reason
    // the segfault does not happen on this calling stack
    // yet removing this call seems to avert it
    // it was not caused by the cleanup code in the Fortran energy library
//  dlclose(dsoHandle);
}


