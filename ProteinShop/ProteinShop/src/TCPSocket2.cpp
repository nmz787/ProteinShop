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
TCPSocket - Wrapper class for TCP sockets ensuring exception safety.
***********************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include "TCPSocket2.h"
#include "Globals.h"

/**************************
Methods of class TCPSocket:
**************************/

TCPSocket::TCPSocket(int portId,int backlog)
{
    /* Create the socket file descriptor: */
    socketFd=socket(PF_INET,SOCK_STREAM,0);
    if(socketFd<0)
        throw std::runtime_error("TCPSocket: Unable to create socket");
    
    /* Bind the socket file descriptor to the port ID: */
    if(portId>=0)
    {
        struct sockaddr_in socketAddress;
        socketAddress.sin_family=AF_INET;
        socketAddress.sin_port=htons(portId);
        socketAddress.sin_addr.s_addr=htonl(INADDR_ANY);
        if(bind(socketFd,(struct sockaddr*)&socketAddress,sizeof(struct sockaddr_in))==-1)
        {
            close(socketFd);
            char msg[80];
            sprintf(msg,"TCPCSocket: Unable to bind socket to port %d",portId);
            throw std::runtime_error(msg);
        }
    }
    
    /* Start listening on the socket: */
    if(listen(socketFd,backlog)==-1)
    {
        close(socketFd);
        throw std::runtime_error("TCPSocket: Unable to start listening on socket");
    }
}

TCPSocket::TCPSocket(const char* hostname,int portId)
{
    /* Create the socket file descriptor: */
    socketFd=socket(PF_INET,SOCK_STREAM,0);
    if(socketFd<0)
        throw std::runtime_error("TCPSocket: Unable to create socket");
    
    /* Bind the socket file descriptor: */
    struct sockaddr_in mySocketAddress;
    mySocketAddress.sin_family=AF_INET;
    mySocketAddress.sin_port=0;
    mySocketAddress.sin_addr.s_addr=htonl(INADDR_ANY);
    if(bind(socketFd,(struct sockaddr*)&mySocketAddress,sizeof(struct sockaddr_in))==-1)
    {
        close(socketFd);
        throw std::runtime_error("TCPSocket: Unable to bind socket to port");
    }
	else
		msg.Debug(NETID, "TCPSocket: bind socket to port");
    
    /* Lookup host's IP address: */
    struct hostent* hostEntry=gethostbyname(hostname);
    if(hostEntry==0)
    {
        close(socketFd);
        char msg[256];
        sprintf(msg,"TCPSocket: Unable to resolve host name %s",hostname);
        throw std::runtime_error(msg);
    }
	else
		{
        sprintf(msg.buf,"TCPSocket: resolve host name %s",hostname);
		msg.Debug(NETID, msg.buf);
    	}
	struct in_addr hostNetAddress;
    hostNetAddress.s_addr=ntohl(((struct in_addr*)hostEntry->h_addr_list[0])->s_addr);

    /* Connect to the remote host: */
    struct sockaddr_in hostAddress;
    hostAddress.sin_family=AF_INET;
    hostAddress.sin_port=htons(portId);
    hostAddress.sin_addr.s_addr=htonl(hostNetAddress.s_addr);
    
	if(connect(socketFd,(const struct sockaddr*)&hostAddress,sizeof(struct sockaddr_in))==-1)
    {
        close(socketFd);
        char msg[256];
        sprintf(msg,"TCPSocket: Unable to connect to host %s on port %d",hostname,portId);
        throw std::runtime_error(msg);
    }
	else
		{
        sprintf(msg.buf,"TCPSocket: connect to host %s on port %d\n",hostname,portId);
		msg.Debug(NETID, msg.buf);
    	}
}

TCPSocket::TCPSocket(const TCPSocket& source)
    :socketFd(dup(source.socketFd))
{
}

TCPSocket::~TCPSocket(void)
{
    close(socketFd);
}

TCPSocket& TCPSocket::operator=(const TCPSocket& source)
{
    if(this!=&source)
    {
        close(socketFd);
        socketFd=dup(source.socketFd);
    }
    return *this;
}

int TCPSocket::getPortId(void) const
{
    struct sockaddr_in socketAddress;
    #ifdef __SGI_IRIX__
    int socketAddressLen=sizeof(struct sockaddr_in);
    #else
    socklen_t socketAddressLen=sizeof(struct sockaddr_in);
    #endif
    getsockname(socketFd,(struct sockaddr*)&socketAddress,&socketAddressLen);
    return ntohs(socketAddress.sin_port);
}

TCPSocket TCPSocket::accept(void) const
{
    /* Wait for connection attempts: */
    int newSocketFd=::accept(socketFd,0,0);
    if(newSocketFd==-1)
        throw std::runtime_error("TCPSocket: Unable to accept connection");
    return TCPSocket(newSocketFd);
}

bool TCPSocket::waitForData(long timeoutSeconds,long timeoutMicroseconds,bool throwException) const
{
    fd_set readFdSet;
    FD_ZERO(&readFdSet);
    FD_SET(socketFd,&readFdSet);
    struct timeval timeout;
    timeout.tv_sec=timeoutSeconds;
    timeout.tv_usec=timeoutMicroseconds;
    bool dataWaiting=select(socketFd+1,&readFdSet,0,0,&timeout)>0;
    if(throwException&&!dataWaiting)
        throw TimeOut("TCPSocket: Time-out while waiting for data");
    return dataWaiting;
}

size_t TCPSocket::read(void* buffer,size_t count)
{
    char* byteBuffer=reinterpret_cast<char*>(buffer);
    ssize_t numBytesRead=::read(socketFd,byteBuffer,count);
    if(errno==EAGAIN||numBytesRead>0)
		#if __BIG_ENDIAN != 1234 // Big endian is not native format!
        return ntohs(size_t(numBytesRead));
		#else
        return size_t(numBytesRead);
		#endif
    else if(numBytesRead==0)
    {
        /* Other end terminated connection: */
        throw PipeError("TCPSocket: Connection terminated by peer during read");
    }
    else
    {
        /* Consider this a fatal error: */
        throw std::runtime_error("TCPSocket: Fatal error during read");
    }
}

void TCPSocket::blockingRead(void* buffer,size_t count)
{
    char* byteBuffer=reinterpret_cast<char*>(buffer);
    while(count>0)
    {
        ssize_t numBytesRead=::read(socketFd,byteBuffer,count);
        if(numBytesRead!=count)
        {
            if(errno==EAGAIN||numBytesRead>0)
            {
                /* Advance result pointer and retry: */
                count-=numBytesRead;
                byteBuffer+=numBytesRead;
            }
            else if(numBytesRead==0)
            {
                /* Other end terminated connection: */
                throw PipeError("TCPSocket: Connection terminated by peer during read");
            }
            else
            {
                /* Consider this a fatal error: */
                throw std::runtime_error("TCPSocket: Fatal error during read");
            }
        }
        else
            count=0;
    }
}

void TCPSocket::blockingWrite(const void* buffer,size_t count)
{
    const char* byteBuffer=reinterpret_cast<const char*>(buffer);
    while(count>0)
    {
        ssize_t numBytesWritten=::write(socketFd,byteBuffer,count);
        if(numBytesWritten!=count)
        {
            /* Check error condition: */
            if(errno==EAGAIN)
            {
                /* Advance data pointer and try again: */
                count-=numBytesWritten;
                byteBuffer+=numBytesWritten;
            }
            else if(errno==EPIPE)
            {
                /* Other end terminated connection: */
                throw PipeError("TCPSocket: Connection terminated by peer during write");
            }
            else
            {
                /* Consider this a fatal error: */
                throw std::runtime_error("TCPSocket: Fatal error during write");
            }
        }
        else
            count=0;
    }
}
