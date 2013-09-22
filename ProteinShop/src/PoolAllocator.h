/***********************************************************************
PoolAllocator - Class to quicklu allocate and release large numbers of
objects of identical size.
Copyright (c) 2003 Oliver Kreylos

This file is part of the Templatized Geometry Library (TGL).

The Templatized Geometry Library is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

The Templatized Geometry Library is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with the Templatized Geometry Library; if not, write to the Free
Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA

In addition, as a special exception, Oliver Kreylos gives permission to
The Regents of the University of California to use, copy, modify,
prepare derivative works, and sublicense under any license the source
code of the Templatized Geometry Library in binary and/or source code
forms so long as such use is in connection with software intended for
protein structure prediction.  You must obey the GNU General Public
License in all respects for all of the code used in other contexts than
software intended for protein structure prediction.  If you received
this file under the GNU General Public License and you modify this file,
you may extend this exception to your version of the file, but you are
not obligated to do so.  If you do not wish to do so, delete this
exception statement from your version.
***********************************************************************/

#ifndef POOLALLOCATOR_INCLUDED
#define POOLALLOCATOR_INCLUDED

template <class Type>
class PoolAllocator
{
    /* Embedded classes: */
    private:
    struct Link // Structure to link free memory items
    {
        /* Elements: */
        public:
        Link* succ;
    };
    struct Chunk // Structure for page-sized chunks of memory
    {
        /* Embedded classes: */
        public:
        enum {size=4*1024-16};
        /* Elements: */
        Chunk* succ; // Next chunk in list
        char mem[size]; // Uninitialized memory
    };
    
    /* Elements: */
    Chunk* chunks; // Pointer to head of chunk list
    Link* head; // Pointer to first free memory item
    
    /* Private methods: */
    PoolAllocator(const PoolAllocator& source) // Forbid copy construction
        :chunks(0),head(0)
    {
    };
    PoolAllocator& operator=(const PoolAllocator& source) // Forbid assignment
    {
        return *this;
    };
    void growPool(void) // Increases memory pool size
    {
        /* Allocate a new chunk: */
        Chunk* newChunk=new Chunk;
        newChunk->succ=chunks;
        chunks=newChunk;
        
        /* Create a linked list of free items in the new chunk: */
        const int numItems=Chunk::size/sizeof(Type);
        char* start=newChunk->mem;
        char* last=newChunk->mem+((numItems-1)*sizeof(Type));
        for(char* item=start;item!=last;item+=sizeof(Type))
            reinterpret_cast<Link*>(item)->succ=reinterpret_cast<Link*>(item+sizeof(Type));
        
        /* Connect the new free list to the old free list: */
        reinterpret_cast<Link*>(last)->succ=head;
        head=reinterpret_cast<Link*>(start);
    };
    
    /* Constructors and destructors: */
    public:
    PoolAllocator(void) // Creates an empty memory pool
        :chunks(0),head(0)
    {
    };
    ~PoolAllocator(void)
    {
        while(chunks!=0)
        {
            Chunk* succ=chunks->succ;
            delete chunks;
            chunks=succ;
        }
    };
    
    /* Methods: */
    void* allocate(void)
    {
        if(head==0) // We ran out of memory
            growPool(); // Grow the pool
        
        /* Return the first free item: */
        Link* result=head;
        head=head->succ;
        
        return reinterpret_cast<Type*>(result);
    };
    
    void free(void* item)
    {
        /* Put the new free item at the head of the list: */
        Link* l=reinterpret_cast<Link*>(item);
        l->succ=head;
        head=l;
    };
};

#endif
