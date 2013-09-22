/***********************************************************************
PriorityHeap - Implementation of a priority queue with a heap structure.
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
***********************************************************************/

#ifndef PRIORITYHEAP_INCLUDED
#define PRIORITYHEAP_INCLUDED

#include <new>

template <class Type> class StdComp // Functor to compare elements of an arbitrary type
{
    /* Methods: */
    public:
    static bool lessEqual(const Type& v1,const Type& v2)
    {
        return v1<=v2;
    };
};

template <class Content,class Comparison =StdComp<Content>,int allocStep =16> class PriorityHeap
{
    friend class Iterator;
    
    /* Embedded classes: */
    public:
    class Iterator
    {
        friend class PriorityHeap;
        
        /* Elements: */
        private:
        PriorityHeap* heap; // Pointer to heap being iterated through
        int index; // Index of element currently referenced
        
        /* Constructors and destructors: */
        public:
        Iterator(void)
            :heap(0)
        {
        };
        private:
        Iterator(PriorityHeap* sHeap,int sIndex)
            :heap(sHeap),index(sIndex)
        {
        };
        public:
        Iterator(const Iterator& source)
            :heap(source.heap),index(source.index)
        {
        };
        
        /* Methods: */
        Iterator& operator=(const Iterator& source)
        {
            heap=source.heap;
            index=source.index;
            return *this;
        };
        friend bool operator==(const Iterator& it1,const Iterator& it2)
        {
            return it1.heap==it2.heap&&it1.index==it2.index;
        };
        friend bool operator!=(const Iterator& it1,const Iterator& it2)
        {
            return it1.heap!=it2.heap||it1.index!=it2.index;
        };
        bool isFinished(void) const
        {
            return heap==0||index<0||index>=heap->numElements;
        };
        const Content& operator*(void) const
        {
            return heap->heap[index];
        };
        Content& operator*(void)
        {
            return heap->heap[index];
        };
        const Content* operator->(void) const
        {
            return heap->heap+index;
        };
        Content* operator->(void)
        {
            return heap->heap+index;
        };
        Iterator& operator++(void)
        {
            ++index;
            return *this;
        };
        Iterator operator++(int)
        {
            Iterator result(*this);
            ++index;
            return result;
        };
        Iterator& operator+=(int increment)
        {
            index+=increment;
            return *this;
        };
        Iterator operator+(int increment) const
        {
            return Iterator(heap,index+increment);
        };
    };
    
    /* Elements: */
    private:
    int allocSize; // Size of allocated heap array
    void* memChunk; // Pointer to uninitialized memory
    int numElements; // Number of elements currently in heap
    Content* heap; // Pointer to heap array
    
    /* Private methods: */
    void reallocate(int newAllocSize)
    {
        /* Allocate a new memory chunk: */
        allocSize=newAllocSize;
        void* newMemChunk=new char[allocSize*sizeof(Content)];
        Content* newHeap=static_cast<Content*>(newMemChunk);
        
        /* Copy all entries from the old heap, then delete the old entries: */
        for(int i=0;i<numElements;++i)
        {
            new(&newHeap[i]) Content(heap[i]);
            heap[i].~Content();
        }
        
        /* Delete the old heap and use the new one: */
        delete[] static_cast<char*>(memChunk);
        memChunk=newMemChunk;
        heap=newHeap;
    };
    
    /* Constructors and destructors: */
    public:
    PriorityHeap(int sAllocSize =0) // Creates empty heap
        :allocSize(sAllocSize<allocStep?allocStep:sAllocSize),memChunk(new char[allocSize*sizeof(Content)]),
         numElements(0),heap(static_cast<Content*>(memChunk))
    {
    };
    PriorityHeap(const PriorityHeap& source)
        :allocSize(source.allocSize),memChunk(new char[allocSize*sizeof(Content)]),
         numElements(source.numElements),heap(static_cast<Content*>(memChunk))
    {
        /* Copy all entries from the source heap: */
        for(int i=0;i<numElements;++i)
            new(&heap[i]) Content(source.heap[i]);
    };
    ~PriorityHeap(void)
    {
        /* Destroy all heap entries: */
        for(int i=0;i<numElements;++i)
            heap[i].~Content();
        
        delete[] static_cast<char*>(memChunk);
    };
    
    /* Methods: */
    PriorityHeap& operator=(const PriorityHeap& source)
    {
        if(this!=&source)
        {
            /* Destroy all heap entries: */
            for(int i=0;i<numElements;++i)
                heap[i].~Content();
            
            /* Allocate a new memory chunk if necessary: */
            if(source.numElements>allocSize)
            {
                delete[] memChunk;
                allocSize=source.allocSize;
                memChunk=new char[allocSize*sizeof(Content)];
                heap=static_cast<Content*>(memChunk);
            }
            
            /* Copy all elements from the source heap: */
            numElements=source.numElements;
            for(int i=0;i<numElements;++i)
                new(&heap[i]) Content(source.heap[i]);
        }
        
        return *this;
    };
    bool isEmpty(void) const
    {
        return numElements==0;
    };
    int getNumElements(void) const
    {
        return numElements;
    };
    Iterator begin(void)
    {
        return Iterator(this,0);
    };
    Iterator end(void)
    {
        return Iterator(this,numElements);
    };
    PriorityHeap& insert(const Content& newElement)
    {
        if(numElements==allocSize)
            reallocate(allocSize+allocStep);
        
        /* Find the correct insertion position for the new element: */
        int insertionPos=numElements;
        int parent=(insertionPos-1)>>1;
        if(parent<0||Comparison::lessEqual(heap[parent],newElement))
        {
            /* Copy the new element to the previous unused slot: */
            new(&heap[insertionPos]) Content(newElement);
        }
        else
        {
            /* Copy the parent to the previous unused slot: */
            new(&heap[insertionPos]) Content(heap[parent]);
            insertionPos=parent;
            while(insertionPos>0)
            {
                parent=(insertionPos-1)>>1;
                if(Comparison::lessEqual(heap[parent],newElement)) // Do the elements have to be swapped?
                    break;
                heap[insertionPos]=heap[parent];
                insertionPos=parent;
            }
            
            /* Copy the new element to the insertion spot: */
            heap[insertionPos]=newElement;
        }
        
        /* Increase the number of stored elements: */
        ++numElements;
        
        return *this;
    };
    const Content& getSmallest(void) const
    {
        return heap[0];
    };
    Content& getSmallest(void)
    {
        return heap[0];
    };
    PriorityHeap& reinsertSmallest(void)
    {
        /* Find the correct insertion position for the (changed) smallest element: */
        int insertionPos=0;
        while(true)
        {
            int child1=(insertionPos<<1)+1;
            int child2=(insertionPos<<1)+2;
            int minIndex=insertionPos;
            if(child1<numElements&&!Comparison::lessEqual(heap[minIndex],heap[child1]))
                minIndex=child1;
            if(child2<numElements&&!Comparison::lessEqual(heap[minIndex],heap[child2]))
                minIndex=child2;
            if(minIndex==insertionPos)
                break;
            Content temp=heap[insertionPos];
            heap[insertionPos]=heap[minIndex];
            heap[minIndex]=temp;
            insertionPos=minIndex;
        }
        return *this;
    };
    PriorityHeap& removeSmallest(void)
    {
        /* Find the correct position to insert the last element: */
        int insertionPos=0;
        while(true)
        {
            int child1=(insertionPos<<1)+1;
            int child2=(insertionPos<<1)+2;
            int minIndex=numElements-1;
            if(child1<numElements&&!Comparison::lessEqual(heap[minIndex],heap[child1]))
                minIndex=child1;
            if(child2<numElements&&!Comparison::lessEqual(heap[minIndex],heap[child2]))
                minIndex=child2;
            if(minIndex==numElements-1)
                break;
            heap[insertionPos]=heap[minIndex];
            insertionPos=minIndex;
        }
        --numElements;
        heap[insertionPos]=heap[numElements];
        heap[numElements].~Content();
        
        return *this;
    };
};

#endif
