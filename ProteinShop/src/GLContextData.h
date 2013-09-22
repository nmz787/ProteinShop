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
GLContextData - Class to store per-GL-context data for application
objects.

***********************************************************************/

#ifndef GLCONTEXTDATA_INCLUDED
#define GLCONTEXTDATA_INCLUDED

#include <HashTable.h>

class GLContextData
{
    /* Embedded classes: */
    public:
    struct DataItem
    {
        /* Constructors and destructors: */
        public:
        virtual ~DataItem(void)
        {
        };
    };
    
    private:
    typedef HashTable<const void*,DataItem*> ItemHash; // Class for hash table mapping pointers to data items
    
    /* Elements: */
    ItemHash context; // A hash table for the context
    
    /* Constructors and destructors: */
    public:
    GLContextData(int sTableSize,float sWaterMark =0.9f,float sGrowRate =1.7312543) // Constructs an empty context
        :context(sTableSize,sWaterMark,sGrowRate)
    {
    };
    ~GLContextData(void)
    {
        /* Delete all data items in this context: */
        for(ItemHash::Iterator it=context.begin();!it.isFinished();++it)
            delete it->getDest();
    };
    
    /* Methods: */
    bool isRealized(const void* thing) const
    {
        return context.isEntry(thing);
    };
    void addDataItem(const void* thing,DataItem* dataItem)
    {
        context.setEntry(ItemHash::Entry(thing,dataItem));
    };
    void removeDataItem(const void *thing)
    {
        context.removeEntry(thing);
    };
    template <class DataItemParam>
    DataItemParam* retrieveDataItem(const void* thing)
    {
        ItemHash::Iterator dataIt=context.findEntry(thing);
        if(dataIt.isFinished())
            return 0;
        else
            return dynamic_cast<DataItemParam*>(dataIt->getDest());
    };
};

#endif
