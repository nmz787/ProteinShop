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

#include <stdio.h>

#include "ClientServerPipe.h"

#include "ProteinClient.h"

/******************************
Methods of class ProteinClient:
******************************/

ProteinClient::ProteinClient(const char* serverName,int serverPort)
    :pipe(new ClientServerPipe(TCPSocket(serverName,serverPort)))
{
    /* Send connect request message: */
    pipe->writeMessage(ClientServerPipe::CONNECT_REQUEST);
    
    /* Wait for connect reply message: */
    if(!pipe->getSocket().waitForData(10,0,false))
        throw ProtocolError("ProteinClient: Server timed out during connection initialization");
	ClientServerPipe::MessageIdType reply=pipe->readMessage();
    
	if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused connection");
    else if(reply!=ClientServerPipe::CONNECT_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message during connection initialization");
}

ProteinClient::~ProteinClient(void)
{
    /* Send disconnect request message: */
    pipe->writeMessage(ClientServerPipe::DISCONNECT_REQUEST);
    delete pipe;
}

void ProteinClient::lockServerQueue(void)
{
    /* Send lock queue request message: */
    pipe->writeMessage(ClientServerPipe::LOCKQUEUE_REQUEST);
    
    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused request to lock queue");
    else if(reply!=ClientServerPipe::LOCKQUEUE_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while locking server queue");
}

void ProteinClient::unlockServerQueue(void)
{
    /* Send unlock queue request message: */
    pipe->writeMessage(ClientServerPipe::UNLOCKQUEUE_REQUEST);
    
    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused request to unlock queue");
    else if(reply!=ClientServerPipe::UNLOCKQUEUE_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while unlocking server queue");
}

int ProteinClient::getNumConfigurations(void)
{
    /* Send get number of configurations request message: */
    pipe->writeMessage(ClientServerPipe::GETNUMCONFIGS_REQUEST);

    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused to send number of configurations");
    else if(reply!=ClientServerPipe::GETNUMCONFIGS_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while querying number of configurations");

    /* Return number of configurations: */
    int numConfigs=pipe->readNumConfigs();
    return numConfigs;
}

unsigned int ProteinClient::getConfiguration(int configIndex,MD::Protein* protein)
{
    /* Send get configuration request message: */
    pipe->writeMessage(ClientServerPipe::GETCONFIG_REQUEST);
    pipe->writeConfigIndex(configIndex);

    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused to send configuration");
    else if(reply!=ClientServerPipe::GETCONFIG_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while retrieving configuration");

    /* Read configuration: */
    if(!pipe->readConfiguration(protein))
        throw ProtocolError("ProteinClient: Mismatching protein size while retrieving configuration");
    
    return pipe->readConfigId();
}

void ProteinClient::removeConfiguration(int configIndex)
{
    /* Send remove configuration request message: */
    pipe->writeMessage(ClientServerPipe::REMOVECONFIG_REQUEST);
    pipe->writeConfigIndex(configIndex);

    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused to remove configuration");
    else if(reply!=ClientServerPipe::REMOVECONFIG_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while removing configuration");
}

void ProteinClient::getOptimizationTree(const char* fileName)
{
    /* Send get optimization tree request message: */
    pipe->writeMessage(ClientServerPipe::GETOPTIMIZATIONTREE_REQUEST);

    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused to send optimization tree");
    else if(reply!=ClientServerPipe::GETOPTIMIZATIONTREE_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while retrieving optimization tree");
    
    /* Retrieve optimization tree: */
    int numTreeNodes=pipe->readNumOptimizationTreeNodes();
    ClientServerPipe::OptimizationTreeNode* treeNodes=new ClientServerPipe::OptimizationTreeNode[numTreeNodes];
    pipe->readOptimizationTreeNodes(numTreeNodes,treeNodes);
    
    /* Open tree file: */
    FILE* file=fopen(fileName,"w");
    if(file==0)
    {
        delete[] treeNodes;
        throw std::runtime_error(std::string("ProteinClient: Unable to open optimization tree file ")+std::string(fileName));
    }

    /* Write tree to file: */
    #if 1
    for(int i=0;i<numTreeNodes;++i)
        fprintf(file,"%d %d %d %g %d\n",treeNodes[i].id,treeNodes[i].parentId,treeNodes[i].manipulated,treeNodes[i].energy,treeNodes[i].status);
    #else
    fprintf(file,"graph G {\n");
    static char* nodeColors[4]={"1.0 0.0 0.0","0.0 1.0 0.0","0.0 0.0 1.0","0.5 0.5 0.5"};
    for(int i=0;i<numTreeNodes;++i)
    {
        char label[80];
        sprintf(label,"%u | %12.4g",treeNodes[i].id,treeNodes[i].energy);
        fprintf(file,"\tcfg%u [label=\"%s\",color=\"%s\"];\n",treeNodes[i].id,label,nodeColors[treeNodes[i].status]);
        if(treeNodes[i].parentId!=0)
        {
            if(treeNodes[i].manipulated)
                fprintf(file,"\tcfg%u -> cfg%u [color=\"1.0 0.0 0.0\"];\n",treeNodes[i].parentId,treeNodes[i].id);
            else
                fprintf(file,"\tcfg%u -> cfg%u;\n",treeNodes[i].parentId,treeNodes[i].id);
        }
    }
    fprintf(file,"}\n");
    #endif
    
    /* Clean up: */
    fclose(file);
    delete[] treeNodes;
}

unsigned int ProteinClient::getBestConfiguration(MD::Protein* protein)
{
    /* Send get configuration request message: */
    pipe->writeMessage(ClientServerPipe::GETBESTCONFIG_REQUEST);

    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused to send configuration");
    else if(reply!=ClientServerPipe::GETBESTCONFIG_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while retrieving configuration");

    /* Read configuration: */
    if(!pipe->readConfiguration(protein))
        throw ProtocolError("ProteinClient: Mismatching protein size while retrieving configuration");
    
    return pipe->readConfigId();
}

unsigned int ProteinClient::addConfiguration(unsigned int parentId,const MD::Protein* protein)
{
    /* Send add configuration request message: */
    pipe->writeMessage(ClientServerPipe::ADDCONFIG_REQUEST);
    
    /* Write configuration: */
    pipe->writeConfigId(parentId);
    pipe->writeConfiguration(protein);

    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused to add configuration");
    else if(reply!=ClientServerPipe::ADDCONFIG_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while adding configuration");
    
    return pipe->readConfigId();
}

void ProteinClient::getConfigurationById(unsigned int configId,MD::Protein* protein)
{
    /* Send get configuration request message: */
    pipe->writeMessage(ClientServerPipe::GETCONFIGBYID_REQUEST);
    pipe->writeConfigId(configId);

    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused to send configuration");
    else if(reply!=ClientServerPipe::GETCONFIGBYID_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while retrieving configuration");

    /* Read configuration: */
    if(!pipe->readConfiguration(protein))
        throw ProtocolError("ProteinClient: Mismatching protein size while retrieving configuration");
}

void ProteinClient::removeConfigurationById(unsigned int configId)
{
    /* Send remove configuration request message: */
    pipe->writeMessage(ClientServerPipe::REMOVECONFIGBYID_REQUEST);
    pipe->writeConfigId(configId);

    /* Wait for reply message: */
    ClientServerPipe::MessageIdType reply=pipe->readMessage();
    if(reply==ClientServerPipe::NACK)
        throw ProtocolWarning("ProteinClient: Server refused to remove configuration");
    else if(reply!=ClientServerPipe::REMOVECONFIGBYID_REPLY)
        throw ProtocolError("ProteinClient: Mismatching message while removing configuration");
}
