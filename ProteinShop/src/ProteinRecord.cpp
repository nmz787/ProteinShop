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


#include "stdint.h"

#include "ProteinRecord.h"
using namespace MD;


ProteinRecord::ProteinRecord (const char *filename) :
    m_file (0),
    m_numAtoms (0),
    m_numEnergyComponents (0),
    m_numRecords (0),
    m_recordSize (0)
{
    if ( !filename ) return;
    m_file = fopen (filename, "r");
    if ( !m_file ) return;
    uint32_t descriptor[4];
    if ( fread(descriptor, 4, 4, m_file) == 4 )
    {
        // sucessfully read the descriptor, try to load the first record
        m_numRecords = descriptor[0];
        m_recordSize = descriptor[1];
        m_numAtoms = descriptor[2];
        m_numEnergyComponents = descriptor[3];
        m_record.setDimSizes (6 + m_numEnergyComponents, m_numAtoms);
        if ( fread(&m_energy, sizeof(double), 1, m_file) != 1 ||
             fread(m_record, m_recordSize, 1, m_file) != 1 )
        {
            // failed to load the first record, reset object
            m_record.setDimSizes (0, 0);
            m_numAtoms = 0;
            m_numEnergyComponents = 0;
            m_numRecords = 0;
            m_recordSize = 0;
            fclose (m_file);
            m_file = 0;
        }
    }
    else
    {
        // could not read the descriptor, close the file
        fclose (m_file);
        m_file = 0;
    }
}


ProteinRecord::~ProteinRecord()
{
    if ( m_file ) fclose (m_file);
}


void ProteinRecord::apply (ProteinState &state) const
{
    if ( isEmpty() || m_numAtoms != state.protein->getNumAtoms() ) return;
    if ( state.energyCalculator &&
         m_numEnergyComponents != energyLibrary->getNumEnergyComponents() )
        return;
    state.protein->setAtomPositions (
        m_record(0, 0),
        m_record(1, 0),
        m_record(2, 0)
    );
    if ( state.energyCalculator && state.energyCalculator->isApplySupported() )
        state.energyCalculator->apply (*this);
}


bool ProteinRecord::load (uint index)
{
    return (
        m_file &&
        index < m_numRecords &&
        fseek(m_file, long(16 + index * m_recordSize), SEEK_SET) == 0 &&
        fread(&m_energy, sizeof(double), 1, m_file) == 1 &&
        fread(m_record, m_recordSize, 1, m_file) == 1
    );
}


void ProteinRecord::snapshot (ProteinState &state)
{
    if ( !state.protein ) return;
    if ( m_file )
    {
        // close the file if there was one
        fclose (m_file);
        m_file = 0;
    }
    // reinitialize the state of this object
    m_numRecords = 0;
    m_numAtoms = state.protein->getNumAtoms();
    m_numEnergyComponents = ( state.energyCalculator )
        ? energyLibrary->getNumEnergyComponents()
        : 0;
    m_record.setDimSizes (6 + m_numEnergyComponents, m_numAtoms);
    m_record.initElements (0.0);
    m_recordSize = m_record.numElements() * sizeof(double);

    // initialize the atom position coordinate arrays
    double *x = m_record;
    double *y = x + m_numAtoms;
    double *z = y + m_numAtoms;
    Position pos;
    for ( Protein::ConstAtomIterator itr = state.protein->atomsBegin();
          itr != state.protein->atomsEnd();
          ++itr )
    {
        // loop through the atoms to get the positions
        pos = itr->getPosition();
        *x = pos[0];
        *y = pos[1];
        *z = pos[2];
        ++x;
        ++y;
        ++z;
    }
    if ( state.energyCalculator )
    {
        // get the energy calculator state
        m_energy = state.energyCalculator->calcEnergy();
        if ( state.energyCalculator->isGradientSupported() )
        {
            // copy the force field gradient vectors
            x = z + m_numAtoms;
            y = x + m_numAtoms;
            z = y + m_numAtoms;
            m_record.bulkCopy (
                x, state.energyCalculator->xGradients(), m_numAtoms
            );
            m_record.bulkCopy (
                y, state.energyCalculator->yGradients(), m_numAtoms
            );
            m_record.bulkCopy (
                z, state.energyCalculator->zGradients(), m_numAtoms
            );
        }
        x = z + m_numAtoms;
        for ( uint i = 0; i < m_numEnergyComponents; ++i )
        {
            // copy the atom energy terms
            m_record.bulkCopy (
                x, state.energyCalculator->getAtomEnergies(i), m_numAtoms
            );
            x += m_numAtoms;
        }
    }
}


