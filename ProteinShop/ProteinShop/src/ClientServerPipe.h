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
ClientServerPipe - Class defining the client/server protocol for
computational steering of protein optimization.
**********************************************************************/

#ifndef CLIENTSERVER_INCLUDED
#define CLIENTSERVER_INCLUDED

#include <iostream>
#include <TCPSocket2.h>

#ifdef CLIENT_SIDE
#include "Protein.h"
#endif

#include "Endianness.h" 
#include "Globals.h"

class ClientServerPipe
{
    /* Embedded classes: */
    public:
    typedef unsigned short int MessageIdType; // Network type for protocol messages
    
    enum MessageId // Enumerated type for protocol messages
    {
        NACK, // Generic negative server reply
        CONNECT_REQUEST, // Request to connect to server
        CONNECT_REPLY, // Positive connect reply
        DISCONNECT_REQUEST, // Polite request to disconnect from server
        LOCKQUEUE_REQUEST, // Request to lock the server's configuration queue
        LOCKQUEUE_REPLY, // Reply that queue has been locked
        UNLOCKQUEUE_REQUEST, // Request to unlock the server's queue
        UNLOCKQUEUE_REPLY, // Reply that queue has been unlocked
        GETNUMCONFIGS_REQUEST, // Request to send the number of configurations in the server's queue
        GETNUMCONFIGS_REPLY, // Reply with number of configurations
        GETCONFIG_REQUEST, // Request to send configuration n
        GETCONFIG_REPLY, // Reply with configuration and ID
        REMOVECONFIG_REQUEST, // Request to remove configuration n from queue
        REMOVECONFIG_REPLY, // Reply that configuration has been removed
        GETOPTIMIZATIONTREE_REQUEST, // Request to send entire optimization tree
        GETOPTIMIZATIONTREE_REPLY, // Reply with optimization tree nodes
        GETBESTCONFIG_REQUEST, // Request to send best configuration
        GETBESTCONFIG_REPLY, // Reply with configuration
        ADDCONFIG_REQUEST, // Request to add configuration to queue
        ADDCONFIG_REPLY, // Reply that configuration has been added, with ID
        GETCONFIGBYID_REQUEST, // Request to send configuration with given ID
        GETCONFIGBYID_REPLY, // Reply with configuration
        REMOVECONFIGBYID_REQUEST, // Request to remove configuration with given ID from queue
        REMOVECONFIGBYID_REPLY // Reply that configuration has been removed
    };
    
    struct OptimizationTreeNode // Structure to communicate optimization tree nodes
    {
        /* Embedded classes: */
        public:
        enum EntryStatus
        {
            PENDING,ACTIVE,RETIRED,CULLED
        };
        
        /* Elements: */
        public:
        unsigned int id; // Unique ID of configuration
        unsigned int parentId; // ID of configuration's parent
        int manipulated; // Flag if configuration was result of interactive manipulation
        double energy; // Configuration's internal energy
        int status; // Status of entry

		static void swapEndianness(OptimizationTreeNode& OptT)
        {
			::swapEndianness(OptT.id);
			::swapEndianness(OptT.parentId);
			::swapEndianness(OptT.manipulated);
			::swapEndianness(OptT.energy);
			::swapEndianness(OptT.status);
        };

    };
    
    /* Elements: */
    private:
    TCPSocket socket; // Socket representing the other end of a connection
    
    /* Constructors and destructors: */
    public:
    ClientServerPipe(const TCPSocket& sSocket) // Creates pipe wrapper for existing socket
        :socket(sSocket)
    {
    };
    
    /* Methods: */
    TCPSocket& getSocket(void) // Returns TCP socket object
    {
        return socket;
    };
    void writeMessage(MessageId messageId) // Writes a protocol message to the pipe
    {
        MessageIdType message=messageId;
        /* Need to do endianness conversion here... */
		#if __BIG_ENDIAN != 1234 // Big endian is not native format!
     	swapEndianness(message);
		msg.Debug(NETID, "writeMessage... Big endian is not native format!"); 
   		#endif
        socket.blockingWrite(&message,sizeof(MessageIdType));
    };
    MessageIdType readMessage(void) // Reads a protocol message from the pipe
    {
		MessageIdType message;
        socket.blockingRead(&message,sizeof(MessageIdType));
        /* Need to do endianness conversion here... */
		#if __BIG_ENDIAN != 1234 // Big endian is not native format!
 		swapEndianness(message);
		msg.Debug(NETID, "readMessage... Big endian is not native format!"); 
   		#endif
        return message;
    };
    template <class DataParam>
    void write(const DataParam& data) // Writes an element of the given data type to the pipe
    {
        /* Need to do endianness conversion here... */
		#if __BIG_ENDIAN != 1234 // Big endian is not native format!
        swapEndianness(data);
		#endif
        socket.blockingWrite(&data,sizeof(DataParam));
    };
    template <class DataParam>
    DataParam read(void) // Reads an element of the given data type from the pipe
    {
        DataParam result;
        socket.blockingRead(&result,sizeof(DataParam));
        /* Need to do endianness conversion here... */
		#if __BIG_ENDIAN != 1234 // Big endian is not native format!
        swapEndianness(result);
   		#endif
        return result;
    };
    template <class DataParam>
    void write(size_t numElements,const DataParam* elements) // Writes an array of elements to the pipe
    {
        /* Need to do endianness conversion here... */
		#if __BIG_ENDIAN != 1234 // Big endian is not native format!
    	/* Swap and send array items individually: */
    	for(size_t i=0;i<numElements;++i)
        {
        	/* Copy and swap the element: */
        	const DataParam* element= &elements[i];
        	swapEndianness(*element);
        	//socket.blockingWrite(&element,sizeof(DataParam));
        }
    		socket.blockingWrite(elements,numElements*sizeof(DataParam));
		#else
    	/* Send array completely and unswapped: */
    	socket.blockingWrite(elements,numElements*sizeof(DataParam));
   		#endif
    };

    template <class DataParam>
    void writeObj(size_t numElements, const DataParam* elements) // Writes an array of object elements to the pipe
    {
       /* Need to do endianness conversion here... */
		#if __BIG_ENDIAN != 1234 // big endian is not native format
    	for(size_t i=0;i<numElements;++i)
        {
        	/* Copy and swap the element: */
        	const DataParam* element= &elements[i];
 			DataParam::swapEndianness(const_cast<DataParam&>(*elements));
			socket.blockingWrite(elements,numElements*sizeof(DataParam));
        }
		#else
    	/* Send array completely and unswapped: */
    	socket.blockingWrite(elements,numElements*sizeof(DataParam));
   		#endif
    };

    template <class DataParam>
    void read(int numElements,DataParam* elements) // Reads an array of elements from the pipe
    {
        socket.blockingRead(elements,numElements*sizeof(DataParam));
        /* Need to do endianness conversion here... */
		#if __BIG_ENDIAN != 1234 // big endian is not native format
    	for(size_t i=0;i<numElements;++i)
        {
        	/* Copy and swap the element: */
        	DataParam* element=&elements[i];
        	swapEndianness(*element);
		}
		#endif
     };

    template <class DataParam>
    void readObj(int numElements,DataParam* elements) // Reads an array of elements from the pipe
    {
		socket.blockingRead(elements,numElements*sizeof(DataParam));
       	/* Need to do endianness conversion here... */
		#if __BIG_ENDIAN != 1234 // big endian is not native format
    	for(size_t i=0;i<numElements;++i)
        {
        	/* Copy and swap the element: */
        	DataParam* element=&elements[i];
 			DataParam::swapEndianness(const_cast<DataParam&>(*element));
		}
		#endif
    };
    void flush(void) // Finishes writing data of a single message to the pipe
    {
    };
    
    /* Higher-level IO methods: */
    void writeNumConfigs(int numConfigs) // Writes number of configurations in queue to pipe
    {
        std::cout<<"writeNumConfigs("<<numConfigs<<")"<<std::endl;
        write(numConfigs);
    };
    int readNumConfigs(void) // Reads number of configurations in queue from pipe
    {
        int numConfigs=read<int>();
        std::cout<<"readNumConfigs() = "<<numConfigs<<std::endl;
        return numConfigs;
    };
    void writeConfigIndex(int configIndex) // Writes a configuration's index in queue to pipe
    {
        std::cout<<"writeConfigIndex("<<configIndex<<")"<<std::endl;
        write(configIndex);
    };
    int readConfigIndex(void) // Reads a configuration's index in queue from pipe
    {
        int configIndex=read<int>();
        std::cout<<"readConfigIndex() = "<<configIndex<<std::endl;
        return configIndex;
    };
    #ifdef CLIENT_SIDE
    void writeConfiguration(const MD::Protein* protein) // Writes protein's configuration to pipe
    {
        /* Create temporary array of atom coordinate components: */
        int numAtoms=protein->getNumAtoms();
		/* the value past to write will be byteswap when endian issue occures */
		/* the numAtoms won't be correct for later use */ 
		write(protein->getNumAtoms()); // Just to make sure...

		double* atomCoords=new double[numAtoms];
        
        /* Write coordinates one component at a time: */
        const MD::Protein::ChainAtom* const* atoms=protein->getAtomPointers();
        for(int dim=0;dim<3;++dim)
        {
            for(int i=0;i<numAtoms;++i)
                atomCoords[i]=double(atoms[i]->getPosition()[dim]);
            write(numAtoms,atomCoords);
        }
        
        delete[] atomCoords;
    };
    #endif
    void writeConfiguration(int numAtoms,const double interlacedCoords[]) // Writes protein's configuration to pipe
    {
        write(numAtoms); // Just to make sure...
        for(int dim=0;dim<3;++dim)
            write(numAtoms,&interlacedCoords[dim*numAtoms]);
    };
    #ifdef CLIENT_SIDE
    bool readConfiguration(MD::Protein* protein) // Reads protein's configuration from pipe
    {
        /* Create temporary array of atom coordinate components: */
        int numAtoms=read<int>();
        double* atomCoords[3];
        for(int dim=0;dim<3;++dim)
        {
            atomCoords[dim]=new double[numAtoms];
            read(numAtoms,atomCoords[dim]);
        }
        
        /* Only update protein's state if numbers of atoms match: */
        if(numAtoms==protein->getNumAtoms())
        {
            MD::Protein::ChainAtom** atoms=protein->getAtomPointers();
            for(int i=0;i<numAtoms;++i)
                atoms[i]->setPosition(MD::Position(atomCoords[0][i],atomCoords[1][i],atomCoords[2][i]));
        }
        
        /* Delete temporary arrays: */
        for(int dim=0;dim<3;++dim)
            delete[] atomCoords[dim];
        
        return numAtoms==protein->getNumAtoms();
    };
    #endif
    bool readConfiguration(int numAtoms,double interlacedCoords[]) // Reads configuration from pipe
    {
        int pipeNumAtoms=read<int>();
        if(pipeNumAtoms==numAtoms)
        {
            /* Update server configuration: */
            for(int dim=0;dim<3;++dim)
                read(numAtoms,&interlacedCoords[dim*numAtoms]);
        }
        else
        {
            /* Read and ignore sent data: */
            double* tempCoords=new double[pipeNumAtoms];
            for(int dim=0;dim<3;++dim)
                read(pipeNumAtoms,tempCoords);
            delete[] tempCoords;
        }
        
        return numAtoms==pipeNumAtoms;
    };
    void writeNumOptimizationTreeNodes(int numTreeNodes) // Writes number of optimization tree nodes to pipe
    {
        write(numTreeNodes);
    };
    int readNumOptimizationTreeNodes(void) // Reads number of optimization tree nodes from pipe
    {
        return read<int>();
    };
    void writeOptimizationTreeNodes(int numTreeNodes,const OptimizationTreeNode nodes[]) // Writes optimization tree to pipe
    {
        writeObj(numTreeNodes,nodes);
    };
    void readOptimizationTreeNodes(int numTreeNodes,OptimizationTreeNode nodes[]) // Reads optimization tree from pipe
    {
        readObj(numTreeNodes,nodes);
    };
    void writeConfigId(unsigned int configId) // Writes a configuration's ID to pipe
    {
        std::cout<<"writeConfigId("<<configId<<")"<<std::endl;
        write(configId);
    };
    int readConfigId(void) // Reads a configuration's ID from pipe
    {
        int configId=read<unsigned int>();
        std::cout<<"readConfigId() = "<<configId<<std::endl;
        return configId;
    };
};

#endif
