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
Endianness - Helper functions to deal with different endianness of basic
data types.
***********************************************************************/

#ifndef ENDIANNESS_INCLUDED
#define ENDIANNESS_INCLUDED

#ifdef __LINUX__
#include <endian.h>
#endif

#include <stdio.h>

/****************************************
In-place endianness conversion functions:
****************************************/

inline void swapEndianness(char& x)
{
    /* Dummy function - no need to swap bytes! */
}

inline void swapEndianness(signed char& x)
{
    /* Dummy function - no need to swap bytes! */
}

inline void swapEndianness(unsigned char& x)
{
    /* Dummy function - no need to swap bytes! */
}

inline void swapEndianness(short int& x)
{
    char* xPtr=reinterpret_cast<char*>(&x);
    char temp;
    
    /* Swap all bytes: */
    temp=xPtr[0];
    xPtr[0]=xPtr[1];
    xPtr[1]=temp;
}

inline void swapEndianness(unsigned short int& x)
{
    char* xPtr=reinterpret_cast<char*>(&x);
    char temp;
    
    /* Swap all bytes: */
    temp=xPtr[0];
    xPtr[0]=xPtr[1];
    xPtr[1]=temp;
}

inline void swapEndianness(int& x)
{
    char* xPtr=reinterpret_cast<char*>(&x);
    char temp;
    
    /* Swap all bytes: */
    temp=xPtr[0];
    xPtr[0]=xPtr[3];
    xPtr[3]=temp;
    temp=xPtr[1];
    xPtr[1]=xPtr[2];
    xPtr[2]=temp;
}

inline void swapEndianness(unsigned int& x)
{
    char* xPtr=reinterpret_cast<char*>(&x);
    char temp;
    
    /* Swap all bytes: */
    temp=xPtr[0];
    xPtr[0]=xPtr[3];
    xPtr[3]=temp;
    temp=xPtr[1];
    xPtr[1]=xPtr[2];
    xPtr[2]=temp;
}

inline void swapEndianness(long& x)
{
    char* xPtr=reinterpret_cast<char*>(&x);
    char temp;
    
    /* Swap all bytes: */
    temp=xPtr[0];
    xPtr[0]=xPtr[3];
    xPtr[3]=temp;
    temp=xPtr[1];
    xPtr[1]=xPtr[2];
    xPtr[2]=temp;
}

inline void swapEndianness(unsigned long& x)
{
    char* xPtr=reinterpret_cast<char*>(&x);
    char temp;
    
    /* Swap all bytes: */
    temp=xPtr[0];
    xPtr[0]=xPtr[3];
    xPtr[3]=temp;
    temp=xPtr[1];
    xPtr[1]=xPtr[2];
    xPtr[2]=temp;
}

inline void swapEndianness(float& x)
{
    char* xPtr=reinterpret_cast<char*>(&x);
    char temp;
    
    /* Swap all bytes: */
    temp=xPtr[0];
    xPtr[0]=xPtr[3];
    xPtr[3]=temp;
    temp=xPtr[1];
    xPtr[1]=xPtr[2];
    xPtr[2]=temp;
}

inline void swapEndianness(double& x)
{
    char* xPtr=reinterpret_cast<char*>(&x);
    char temp;
    
    /* Swap all bytes: */
    temp=xPtr[0];
    xPtr[0]=xPtr[7];
    xPtr[7]=temp;
    temp=xPtr[1];
    xPtr[1]=xPtr[6];
    xPtr[6]=temp;
    temp=xPtr[2];
    xPtr[2]=xPtr[5];
    xPtr[5]=temp;
    temp=xPtr[3];
    xPtr[3]=xPtr[4];
    xPtr[4]=temp;
}

inline void swapEndianness(const unsigned short int& x)
{
	swapEndianness(const_cast<unsigned short int&>(x));
}

inline void swapEndianness(const int& x)
{
	swapEndianness(const_cast<int&>(x));
}

inline void swapEndianness(const unsigned int& x)
{
	swapEndianness(const_cast<unsigned int&>(x));
}

inline void swapEndianness(const long& x)
{
	swapEndianness(const_cast<long&>(x));
}

inline void swapEndianness(const unsigned long& x)
{
	swapEndianness(const_cast<unsigned long&>(x));
}

inline void swapEndianness(const float& x)
{
	swapEndianness(const_cast<float&>(x));
}

inline void swapEndianness(const double& x)
{
	swapEndianness(const_cast<double&>(x));
}

/*******************
Conversion routines:
*******************/

template <class Type>
Type toLE(Type value)
{
    #if __LITTLE_ENDIAN!=1234
    swapEndianness(value);
    #endif
    
    return value;
}

template <class Type>
Type toBE(Type value)
{
    #if __BIG_ENDIAN!=1234
    swapEndianness(value);
    #endif
    
    return value;
}

template <class Type>
Type fromLE(Type value)
{
    #if __LITTLE_ENDIAN!=1234
    swapEndianness(value);
    #endif
    
    return value;
}

template <class Type>
Type fromBE(Type value)
{
    #if __BIG_ENDIAN!=1234
    swapEndianness(value);
    #endif
    
    return value;
}

/*************************************************
File I/O functions for little-endian binary files:
*************************************************/

template <class Type>
void writeLE(FILE* file,Type value)
{
    #if __LITTLE_ENDIAN!=1234
    swapEndianness(value);
    #endif
    
    fwrite(&value,sizeof(Type),1,file);
}

template <class Type>
void writeLE(FILE* file,const Type* values,int numValues)
{
    #if __LITTLE_ENDIAN==1234
    fwrite(values,sizeof(Type),numValues,file);
    #else
    for(int i=0;i<numValues;++i)
    {
        Type temp=values[i];
        swapEndianness(temp);
        fwrite(&temp,sizeof(Type),1,file);
    }
    #endif
}

template <class Type>
Type readLE(FILE* file)
{
    Type result;
    fread(&result,sizeof(Type),1,file);
    
    #if __LITTLE_ENDIAN!=1234
    swapEndianness(result);
    #endif
    
    return result;
}

template <class Type>
void readLE(Type* values,int numValues,FILE* file)
{
    fread(values,sizeof(Type),numValues,file);
    
    #if __LITTLE_ENDIAN!=1234
    for(int i=0;i<numValues;++i)
        swapEndianness(values[i]);
    #endif
}

/**********************************************
File I/O functions for big-endian binary files:
**********************************************/

template <class Type>
void writeBE(FILE* file,Type value)
{
    #if __BIG_ENDIAN!=1234
    swapEndianness(value);
    #endif
    
    fwrite(&value,sizeof(Type),1,file);
}

template <class Type>
void writeBE(FILE* file,const Type* values,int numValues)
{
    #if __BIG_ENDIAN==1234
    fwrite(values,sizeof(Type),numValues,file);
    #else
    for(int i=0;i<numValues;++i)
    {
        Type temp=values[i];
        swapEndianness(temp);
        fwrite(&temp,sizeof(Type),1,file);
    }
    #endif
}

template <class Type>
Type readBE(FILE* file)
{
    Type result;
    fread(&result,sizeof(Type),1,file);
    
    #if __BIG_ENDIAN!=1234
    swapEndianness(result);
    #endif
    
    return result;
}

template <class Type>
void readBE(Type* values,int numValues,FILE* file)
{
    fread(values,sizeof(Type),numValues,file);
    
    #if __BIG_ENDIAN!=1234
    for(int i=0;i<numValues;++i)
        swapEndianness(values[i]);
    #endif
}

#endif
