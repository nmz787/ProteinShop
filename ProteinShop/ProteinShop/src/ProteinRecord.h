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
ProteinRecord - A record object encapsulating the state of a protein at
a point in time.

***********************************************************************/


class ProteinRecord;


#ifndef PROTEIN_RECORD_INCLUDED
#define PROTEIN_RECORD_INCLUDED


#include <cstdio>
#include <cstdlib>

#include "Buffer.h"
#include "Globals.h"


/** Class to access recorded protein state information from a binary file.  The
    ability to create this file is not provided by this class; an implementation
    is available in the AMBER plug-in.  The file format is simple:
    <ul>
        <li>4 bytes (uint):  N = number of records.
        <li>4 bytes (uint):  S = size of each record in bytes
            = 8 + 8 * A * (6 + E).
        <li>4 bytes (uint):  A = number of atoms.
        <li>4 bytes (uint):  E = number of energy components.
        <li>N * S bytes (double):  N records.  Each record consists of a single
            double-precision floating point value (the total energy) followed by
            6 + E arrays, each containing A double-precision floating point
            numbers.  The order of the arrays is:
            <ol>
                <li>Atom position X-coordinates.
                <li>Atom position Y-coordinates.
                <li>Atom position Z-coordinates.
                <li>Force field gradient X-components.
                <li>Force field gradient Y-components.
                <li>Force field gradient Z-components.
                <li>First atom energy term.
                <li>Second atom energy term.
                <li>Et-cetera.
                <li>Final (E^th) atom energy term.
            </ol>
    </ul>
    The gentle reader might be interested to note that the reason why we use a
    column-major storage format is inherited from the AmberEngine component of
    the AMBER plug-in.
*/
class ProteinRecord
{
public:

    /** Create a file record access object for a particular protein.
        @param fileName
            The name of the file containing the records.  If the name is null or
            the file does not exist, no records will be loaded and this object
            will be empty unless the snapshot() method is called.  If the file
            exists, the first record will be loaded from it. */
    ProteinRecord (const char *fileName = 0);

    /** Close the file and delete the current record, if any. */
    ~ProteinRecord();

    /** Apply the state in this object's record to the specified protein.
        @param state
            The protein to load the current record's state into.  If the number
            of atoms does not match, if the number of energy components does not
            match, or if there is no record loaded, this call has no effect.
        @return
            True if the call succeeds and both the EnergyCalculator and the
            Protein reflect the contents of the record.  False if some condition
            was not met, in which case the state is not modified. */
    void apply (ProteinState &state) const;

    /** Get the atom energies for a specified component.
        @param index
            The index of the energy component, with a value less than
            EnergyLibrary::getNumEnergyComponents().
        @return
            The requested atom energies, or null if @a index is out of range or
            the record is empty.  If there was no energy calculator for the
            protein, these values are all zero. */
    const double *atomEnergies (uint index) const
    {
        return m_record (6 + index, 0);
    }

    /** Get the total energy for the molecule.
        @return
            The optimization target value. */
    double energy() const { return m_energy; }

    /** Determine whether or not this object is empty.
        @return
            True if this object contains a record, false if not. */
    bool isEmpty() const { return ( !m_record.numElements() ); }

    /** Load a specified record from the file.
        @param index
            In the range [0,numRecords()-1].
        @return
            True if the record was loaded, false if @a index is out of range or
            if there is no file. */
    bool load (uint index);

    /** Get the number of atoms in this object.
        @return
            The number of atoms determined by the input file or snapshot(). */
    uint numAtoms() const { return m_numAtoms; }

    /** Get the number of energy components in this object.
        @return
            The number of energy term (e.g. nonbonded, dihedral, etc.) arrays in
            each record. */
    uint numEnergyComponents() const { return m_numEnergyComponents; }

    /** Get the number of records available.
        @return
            The number of states recorded in the file, 0 if there is no file. */
    uint numRecords() const { return m_numRecords; }

    /** Create a record inside this object that represents the current state of
        a specified protein.  If this object is associated with a file, that
        file is closed, and subsequent calls to load() will return false.
        @param state
            The the protein to take a snapshot of. */
    void snapshot (ProteinState &state);
    
    /** Get the atom X-coordinates for this record.
        @return
            The () X-coordinates of the atom positions in an array of numAtoms()
            elements; null if this object is empty. */
    const double *xCoordinates() const { return m_record (0, 0); }
    
    /** Get the X-components of the atom energy gradients for this record.
        @return
            The X-components of the force field gradients in an array of
            numAtoms() elements.  If the calculator did not support gradients,
            these values are zero.  If this object is empty, return null. */
    const double *xGradients() const { return m_record (3, 0); }

    /** Get the atom Y-coordinates for this record.
        @return
            The Y-coordinates of the atom positions in an array of numAtoms()
            elements; null if this object is empty. */
    const double *yCoordinates() const { return m_record (1, 0); }

    /** Get the Y-components of the atom energy gradients for this record.
        @return
            The Y-components of the force field gradients in an array of
            numAtoms() elements.  If the calculator did not support gradients,
            these values are zero.  If this object is empty, return null. */
    const double *yGradients() const { return m_record (4, 0); }

    /** Get the atom Z-coordinates for this record.
        @return
            The Z-coordinates of the atom positions in an array of numAtoms()
            elements; null if this object is empty. */
    const double *zCoordinates() const { return m_record (2, 0); }
    
    /** Get the Z-components of the atom energy gradients for this record.
        @return
            The Z-components of the force field gradients in an array of
            numAtoms() elements.  If the calculator did not support gradients,
            these values are zero.  If this object is empty, return null. */
    const double *zGradients() const { return m_record (5, 0); }

private:

    double m_energy;            // the total energy in the state record
    FILE *m_file;               // the record file, if any
    uint m_numAtoms;            // number atoms in the protein
    uint m_numEnergyComponents; // number of energy terms
    uint m_numRecords;          // number of records in the file
    Buffer<double,2> m_record;  // the protein state record
    uint m_recordSize;          // size of each record in bytes
};


#endif
