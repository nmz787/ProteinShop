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


class EnergyOptimizer;
class EnergyCalculator;
class EnergyLibrary;


#ifndef ENERGYAPI_INCLUDED
#define ENERGYAPI_INCLUDED


#include <cstdlib>

#include "ProteinRecord.h"


/* Forward declarations: */
namespace MD { class Protein; }


/** Class representing specific energy optimizer for a single protein.  There
    should be only one instance of this class for each EnergyCalculator, and
    that instance should be destroyed when the EnergyCalculator that created it
    is destroyed. */
class EnergyOptimizer
{
public:

    /** Virtual destructor providing linkage to implementation destructor. */
    virtual ~EnergyOptimizer() {}

    /** Cancel an ongoing optimization in progress.  If no optimization is in
        progress, this call has no effect. */
    virtual void cancel() = 0;

    /** Release the optimizer thread mutex after a call to isUpdateAvailable()
        that returned true. */
    virtual void finishUpdate() = 0;

    /// Return true if minimization states are being recorded in a file.
    virtual bool isRecording() = 0;

    /// Return true if the recording of minimization states is supported.
    virtual bool isRecordingSupported() = 0;

    /** Return true if the optimizer is currently running on this instance. */
    virtual bool isRunning() = 0;

    /** Grab the optimizer thread mutex and check to see if an updated set of
        atom coordinates are available.  If an update is available, call
        *Coordinates() to get the atom coordinate arrays, then call
        finishUpdate() to release the mutex.  If no update is available, the
        mutex is released before this call returns.
        @return
            True if updated coordinates are available, false if not.  Any
            subsequent call to this method will return false until another
            update has become available. */
    virtual bool isUpdateAvailable() = 0;

    /** Start the optimizer for this protein.  If the optimizer is already
        running for another protein, this call has no effect.
        @return
            True if the optimization started, false if not. */
    virtual bool optimize() = 0;

    /** Specify the recording file, or disable recording if filename is null,
        and do nothing if recording is not supported. */
    virtual void record (const char *fileName) = 0;

    /** Get the updated atom coordinates.  This call is unsynchronized, but
        should be bracketed by calls to isUpdateAvailable() and finishUpdate().
        @return
            The X-coordinates of the optimized protein atoms. */
    virtual const double *xCoordinates() = 0;

    /** Get the updated atom coordinates.  This call is unsynchronized, but
        should be bracketed by calls to isUpdateAvailable() and finishUpdate().
        @return
            The Y-coordinates of the optimized protein atoms. */
    virtual const double *yCoordinates() = 0;

    /** Get the updated atom coordinates.  This call is unsynchronized, but
        should be bracketed by calls to isUpdateAvailable() and finishUpdate().
        @return
            The Z-coordinates of the optimized protein atoms. */
    virtual const double *zCoordinates() = 0;
};


/** Class representing specific energy calculator for a single protein.  If an
    optimizer is implemented, this class won't do much of anything while it is
    running.  However, if the optimizer is locked (i.e. not running) during an
    update interval, this class will work.  Energy will be in a constant state
    of flux while the optimizer is running anyway. */
class EnergyCalculator
{
public:

    /// Destroys energy calculator and all resources it allocated.
    virtual ~EnergyCalculator() {}

    /** Load the state variables of the calculator from a record.  Do nothing if
        this function is not supported, or if a minimization is underway. */
    virtual void apply (const ProteinRecord &state) = 0;

    /// Calculates and returns total energy value for current atom coordinates.
    virtual double calcEnergy() = 0;

    /** Returns per-atom energy of given atom, using current energy component
        states. */
    virtual double getAtomEnergy (int atomIndex) = 0;

    /// Get all atom energies for a specified component in an array.
    virtual const double *getAtomEnergies (int energyComponentIndex) = 0;

    /// Returns current state of energy component.
    virtual bool getEnergyComponentState (int energyComponentIndex) = 0;

    /// Return true if the apply() function is supported.
    virtual bool isApplySupported() = 0;

    /// Return true if gradients are available in this implementation.
    virtual bool isGradientSupported() = 0;

    /// Get the optimizer object for the protein, null if none is implemented.
    virtual EnergyOptimizer *optimizer() = 0;

    /// Enables or disables an energy component for subsequent per-atom queries.
    virtual void setEnergyComponentState (int energyComponentIndex,
                                          bool enable) = 0;

    /// Tells the energy calculator that atom coordinates have changed.
    virtual void updateProtein() = 0;

    /** Get the X component gradients from the most rencent calculation.  If
        this calculator does not support gradients, return null. */
    virtual const double *xGradients() = 0;

    /** Get the Y component gradients from the most rencent calculation.  If
        this calculator does not support gradients, return null. */
    virtual const double *yGradients() = 0;

    /** Get the Z component gradients from the most rencent calculation.  If
        this calculator does not support gradients, return null. */
    virtual const double *zGradients() = 0;
};


/// Class representing dynamically loaded energy computation libraries.
class EnergyLibrary
{
    /// Handle of dynamic shared object containing energy computation library.
    void *dsoHandle;

protected:

    /** Constructs an energy library object; called from static DSO loader
        function. */
    EnergyLibrary (void* sDsoHandle);

public:

    /// Function signature for energy library creation function inside DSO.
    typedef EnergyLibrary* (*CreateEnergyLibraryFunction)(void*);

    /// Opens a DSO containing an energy library and returns library object.
    static EnergyLibrary* load (const char* dsoName);

    /// Destroys library object and closes DSO.
    virtual ~EnergyLibrary();

    /// Returns number of energy components reported by energy library.
    virtual int getNumEnergyComponents() const = 0;

    /** Returns user-readable name for each energy component reported by energy
        library. */
    virtual const char* getEnergyComponentName (int index) const = 0;

    /// Creates energy calculator for given protein.
    virtual EnergyCalculator* createEnergyCalculator (const MD::Protein* p) = 0;
};


#endif
