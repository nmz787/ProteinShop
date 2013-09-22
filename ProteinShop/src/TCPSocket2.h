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

#ifndef TCPSOCKET_INCLUDED
#define TCPSOCKET_INCLUDED

#include <string>
#include <stdexcept>

class TCPSocket
{
    /* Embedded classes: */
    public:
    class PipeError:public std::runtime_error // Exception for unexpected connection termination
    {
        /* Constructors and destructors: */
        public:
        PipeError(const std::string& what_arg)
            :std::runtime_error(what_arg)
        {
        };
    };
    class TimeOut:public std::runtime_error // Exception for time-outs when waiting for data
    {
        /* Constructors and destructors: */
        public:
        TimeOut(const std::string& what_arg)
            :std::runtime_error(what_arg)
        {
        };
    };
    
    /* Elements: */
    private:
    int socketFd; // Internal socket file descriptor
    
    /* Constructors and destructors: */
    private:
    TCPSocket(int sSocketFd) // Creates a TCPSocket wrapper around an existing socket file descriptor (without copying)
        :socketFd(sSocketFd)
    {
    };
    public:
    TCPSocket(int portId,int backlog); // Creates a socket on the local host and starts listening; if portId is negative, random free port is assigned
    TCPSocket(const char* hostname,int portId); // Creates a socket connected to a remote host
    TCPSocket(const TCPSocket& source); // Copy constructor
    ~TCPSocket(void); // Closes a socket
    
    /* Methods: */
    int getFd(void) // Returns low-level socket file descriptor
    {
        return socketFd;
    };
    TCPSocket& operator=(const TCPSocket& source); // Assignment operator
    int getPortId(void) const; // Returns port ID assigned to a socket
    TCPSocket accept(void) const; // Waits for an incoming connection on a listening socket and returns a new socket connected to the initiator
    
    /* I/O methods: */
    bool waitForData(long timeoutSeconds,long timeoutMicroseconds,bool throwException =true) const; // Waits for incoming data on TCP socket; returns true if data is ready; (optionally) throws exception if wait times out
    size_t read(void* buffer,size_t count); // Reads raw buffer from TCP socket; returns number of bytes read
    void blockingRead(void* buffer,size_t count); // Reads raw buffer from TCP socket; blocks until data completely read
    void blockingWrite(const void* buffer,size_t count); // Writes raw buffer to TCP socket; blocks until data completely written
};

#endif
