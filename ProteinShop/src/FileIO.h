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
FileIO - Small set of template functions to make low-level file I/O a
bit more readable.
***********************************************************************/

#ifndef FILEIO_INCLUDED
#define FILEIO_INCLUDED

#include <stdio.h>

template <class Type> inline Type swapHiLo(const Type& source)
{
    Type result;
    const unsigned char* sPtr=(const unsigned char*)&source;
    unsigned char* rPtr=((unsigned char*)&result)+(sizeof(Type)-1);
    for(int i=0;i<sizeof(Type);++i,++sPtr,--rPtr)
        *rPtr=*sPtr;
    return result;
}

template <class Type> inline size_t fread(Type& destination,FILE* f)
{
    return fread(&destination,sizeof(Type),1,f);
}

template <class Type> inline size_t fread(Type* destination,size_t elements,FILE* f)
{
    return fread(destination,sizeof(Type),elements,f);
}

template <class Type> inline size_t fwrite(const Type& source,FILE* f)
{
    return fwrite(&source,sizeof(Type),1,f);
}

template <class Type> inline size_t fwrite(const Type* source,size_t elements,FILE* f)
{
    return fwrite(source,sizeof(Type),elements,f);
}

/************************************************
Endian-safe file input/output for integral types:
************************************************/

template <class Type> inline size_t freadLE(Type& destination,FILE* f)
{
    unsigned char buffer[sizeof(Type)];
    size_t result=fread(buffer,1,sizeof(Type),f);
    destination=0;
    for(int i=sizeof(Type)-1;i>=0;--i)
    {
        destination<<=8;
        destination|=(Type)buffer[i];
    }
    return result;
}

template <class Type> inline size_t fwriteLE(const Type& source,FILE* f)
{
    unsigned char buffer[sizeof(Type)];
    Type temp=source;
    for(int i=0;i<sizeof(Type);++i)
    {
        buffer[i]=(unsigned char)(temp&0xff);
        temp>>=8;
    }
    return fwrite(buffer,1,sizeof(Type),f);
}

#endif
