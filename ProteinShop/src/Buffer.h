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

#ifndef BUFFER_INCLUDED
#define BUFFER_INCLUDED


#include <cstdarg>


/** A simple class for the management of dynamic multidimensional arrays.  It
    does not provide for the initialization of data within the buffer when it is
    created, or for the preservation of data within the buffer when it grows. */
template <typename Element, uint NumDimensions = 1> class Buffer
{
private:

    Element *m_buffer;
    uint m_bufferSize;
    uint m_dimSize[NumDimensions];
    uint m_dimStep[NumDimensions];
    uint m_numElements;
    bool m_rowMajor;

    void initStepsAndBuffer()
    {
        uint requiredBufferSize = 1;
        if ( m_rowMajor )
        {
            // recalculate steps for row-major order
            for ( uint i = 0; i < NumDimensions; ++i )
            {
                requiredBufferSize *= m_dimSize[i];
                m_dimStep[i] = 1;
                for ( uint j = i + 1; j < NumDimensions; ++j )
                    m_dimStep[i] *= m_dimSize[j];
            }
        }
        else
        {
            // recalculate steps for column-major order
            for ( uint i = 0; i < NumDimensions; ++i )
            {
                requiredBufferSize *= m_dimSize[i];
                m_dimStep[i] = 1;
                for ( uint j = 0; j < i; ++j )
                    m_dimStep[i] *= m_dimSize[j];
            }
        }
        reserve (requiredBufferSize);
        m_numElements = requiredBufferSize;
    }

public:

    Buffer() :
        m_buffer (0),
        m_bufferSize (0),
        m_numElements (0),
        m_rowMajor (true)
    {
        for ( uint i = 0; i < NumDimensions; ++i )
            m_dimStep[i] = m_dimSize[i] = 0;
    }

    Buffer (uint firstDimSize, ...) :
        m_buffer (0),
        m_bufferSize (0),
        m_numElements (0),
        m_rowMajor (true)
    {
        m_dimSize[0] = firstDimSize;
        if ( NumDimensions > 1 )
        {
            va_list sizes;
            va_start (sizes, firstDimSize);
            for ( uint i = 1; i < NumDimensions; ++i )
                m_dimSize[i] = va_arg (sizes, uint);
            va_end (sizes);
        }
        initStepsAndBuffer();
    }

    ~Buffer()
    {
        delete[] m_buffer;
        m_buffer = 0;
        m_bufferSize = 0;
        for ( uint i = 0; i < NumDimensions; ++i )
            m_dimStep[i] = m_dimSize[i] = 0;
    }

    /// Get the capacity of the storage currently allocated for the buffer.
    uint capacity() const
    {
        return m_bufferSize;
    }

    /// Copy the entire contents of this buffer from an array.
    void copyFrom (const Element *src)
    {
        bulkCopy (m_buffer, src, m_numElements);
    }

    /// Copy the entire contents of this buffer to an array.
    void copyTo (Element *dest) const
    {
        bulkCopy (dest, m_buffer, m_numElements);
    }

    /// Get the size of a specified dimension.
    uint dimSize (uint dimNumber) const
    {
        return m_dimSize[dimNumber];
    }

    /** Get the number of physical elements between consecutive index steps
        along a particular dimension within the buffer. */
    uint dimStep (uint dimNumber) const
    {
        return m_dimStep[dimNumber];
    }

    /// Calculate a linear index from a set of dimensional indices.
    uint index (uint firstDimIndex, ...) const
    {
        uint linearIndex = firstDimIndex * m_dimStep[0];
        if ( NumDimensions > 1 )
        {
            va_list indices;
            va_start (indices, firstDimIndex);
            for ( uint i = 1; i < NumDimensions; ++i )
                linearIndex += m_dimStep[i] * va_arg(indices, uint);
            va_end (indices);
        }
        return linearIndex;
    }

    /// Initialize the elements of this buffer.
    void initElements (const Element &value)
    {
        bulkSet (m_buffer, value, m_numElements);
    }

    /// Query the physical storage format.
    bool isRowMajor() const
    {
        return m_rowMajor;
    }

    /// Get the total number of elements in the buffer.
    uint numElements() const
    {
        return m_numElements;
    }

    /// Get lvalue access to a buffer element, with range checking.
    Element *operator() (uint firstDimIndex, ...)
    {
        if ( m_buffer )
        {
            uint linearIndex = firstDimIndex * m_dimStep[0];
            if ( NumDimensions > 1 )
            {
                va_list indices;
                va_start (indices, firstDimIndex);
                for ( uint i = 1; i < NumDimensions; ++i )
                    linearIndex += m_dimStep[i] * va_arg(indices, uint);
                va_end (indices);
            }
            if ( linearIndex < m_numElements )
                return ( m_buffer + linearIndex );
        }
        return 0;
    }

    /// Get rvalue access to a buffer element, with range checking.
    const Element *operator() (uint firstDimIndex, ...) const
    {
        if ( m_buffer )
        {
            uint linearIndex = firstDimIndex * m_dimStep[0];
            if ( NumDimensions > 1 )
            {
                va_list indices;
                va_start (indices, firstDimIndex);
                for ( uint i = 1; i < NumDimensions; ++i )
                    linearIndex += m_dimStep[i] * va_arg(indices, uint);
                va_end (indices);
            }
            if ( linearIndex < m_numElements )
                return ( m_buffer + linearIndex );
        }
        return 0;
    }

    /// Get direct lvalue access to the buffer.
    operator Element* ()
    {
        return m_buffer;
    }

    /// Get direct rvalue access to the buffer.
    operator const Element* () const
    {
        return m_buffer;
    }

    /** Ensure that the buffer has sufficient storage for the indicated number
        of elements. */
    void reserve (uint capacity)
    {
        if ( capacity > m_bufferSize )
        {
            // the buffer needs to get larger
            uint newBufferSize = (m_bufferSize) ? m_bufferSize : 1;
            while ( newBufferSize < capacity )
                newBufferSize *= 2;
            Element *newBuffer = new Element[newBufferSize];
            m_bufferSize = newBufferSize;
            delete[] m_buffer;
            m_buffer = newBuffer;
        }
    }

    /// Change the size of one of this buffer's dimensions.
    void setDimSize (uint dimNumber, uint size)
    {
        if ( dimNumber < NumDimensions )
        {
            m_dimSize[dimNumber] = size;
            initStepsAndBuffer();
        }
    }

    /// Change the size of this buffer's dimensions.
    void setDimSizes (uint firstDimSize, ...)
    {
        m_dimSize[0] = firstDimSize;
        if ( NumDimensions > 1 )
        {
            va_list sizes;
            va_start (sizes, firstDimSize);
            for ( uint i = 1; i < NumDimensions; ++i )
                m_dimSize[i] = va_arg (sizes, uint);
            va_end (sizes);
        }
        initStepsAndBuffer();
    }

    /// Select the physical storage format used by this buffer.
    void setRowMajor (bool rowMajor)
    {
        if ( rowMajor != m_rowMajor )
        {
            m_rowMajor = rowMajor;
            initStepsAndBuffer();
        }
    }

    /** Efficiently copy the elements of one array to another.  If the arrays
        overlap, the results are undefined. */
    static void bulkCopy (Element *dest, const Element *src, uint count)
    {
        if ( !dest || !src || !count ) return;
        while ( count >= 16 )
        {
            count -= 16;
            dest[0] = src[0];
            dest[1] = src[1];
            dest[2] = src[2];
            dest[3] = src[3];
            dest[4] = src[4];
            dest[5] = src[5];
            dest[6] = src[6];
            dest[7] = src[7];
            dest[8] = src[8];
            dest[9] = src[9];
            dest[10] = src[10];
            dest[11] = src[11];
            dest[12] = src[12];
            dest[13] = src[13];
            dest[14] = src[14];
            dest[15] = src[15];
            dest += 16;
            src += 16;
        }
        for ( uint i = 0; i < count; ++i )
            dest[i] = src[i];
    }

    /// Efficiently initialize the elements of an array.
    static void bulkSet (Element *dest, const Element &value, uint count)
    {
        if ( !dest || !count ) return;
        while ( count >= 16 )
        {
            count -= 16;
            dest[0] = value;
            dest[1] = value;
            dest[2] = value;
            dest[3] = value;
            dest[4] = value;
            dest[5] = value;
            dest[6] = value;
            dest[7] = value;
            dest[8] = value;
            dest[9] = value;
            dest[10] = value;
            dest[11] = value;
            dest[12] = value;
            dest[13] = value;
            dest[14] = value;
            dest[15] = value;
            dest += 16;
        }
        for ( uint i = 0; i < count; ++i )
            dest[i] = value;
    }
};


#endif


