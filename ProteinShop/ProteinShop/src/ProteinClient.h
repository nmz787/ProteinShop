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
ProteinClient - Class representing the client side of a computational
steering environment for protein optimization.
***********************************************************************/

#ifndef PROTEINCLIENT_INCLUDED
#define PROTEINCLIENT_INCLUDED

#include <stdexcept>

/* Forward declarations: */
class ClientServerPipe;
namespace MD {
class Protein;
}

class ProteinClient
{
    /* Embedded classes: */
    public:
    class ProtocolWarning:public std::runtime_error // Exception for recoverable protocol errors
    {
        /* Constructors and destructors: */
        public:
        ProtocolWarning(const std::string& what_arg)
            :std::runtime_error(what_arg)
        {
        };
    };
    
    class ProtocolError:public std::runtime_error // Exception when unexpected protocol messages are received
    {
        /* Constructors and destructors: */
        public:
        ProtocolError(const std::string& what_arg)
            :std::runtime_error(what_arg)
        {
        };
    };
    
    /* Elements: */
    private:
    ClientServerPipe* pipe; // Pipe connected to the protein server
    
    /* Constructors and destructors: */
    public:
    ProteinClient(const char* serverName,int serverPort); // Connects client to server
    ~ProteinClient(void); // Disconnects client from server
    
    void lockServerQueue(void); // Locks the server queue for further access
    void unlockServerQueue(void); // Unlocks server queue to allow optimization to continue
    
    /* Methods that require locked server queue: */
    int getNumConfigurations(void); // Returns number of configurations in server's queue
    unsigned int getConfiguration(int configIndex,MD::Protein* protein); // Retrieves a configuration from server's queue and updates the given protein
    void removeConfiguration(int configIndex); // Removes a configuration from server's queue
    
    /* Methods that work on unlocked server queue: */
    void getOptimizationTree(const char* fileName); // Retrieves entire optimization tree from server and writes it as a graph drawer file
    unsigned int getBestConfiguration(MD::Protein* protein); // Retrieves best configuration from server's queue and updates the given protein
    unsigned int addConfiguration(unsigned int parentId,const MD::Protein* protein); // Adds a new configuration to server's queue
    void getConfigurationById(unsigned int configId,MD::Protein* protein); // Retrieves a configuration by its ID
    void removeConfigurationById(unsigned int configId); // Removes a configuration by its ID
};

#endif
