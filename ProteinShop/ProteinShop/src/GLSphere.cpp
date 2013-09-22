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
GLSphere - Function to render spheres in different modes. Hack job
directly included by client code until I come up with a better way.
***********************************************************************/

#include <Math/Math.h>
#include <GLTemplates.h>
#include <GLTypes.h>

inline void combine(const GLVector<GLdouble,3>& p100,const GLVector<GLdouble,3>& p010,const GLVector<GLdouble,3>& p001,double w0,double w1,double radius)
{
    double w2=1.0-w0-w1;
    GLVector<GLdouble,3> result;
    double resultLen=0.0;
    for(int i=0;i<3;++i)
    {
        result[i]=p100[i]*w0+p010[i]*w1+p001[i]*w2;
        resultLen+=Math::sqr(result[i]);
    }
    resultLen=Math::sqrt(resultLen);
    for(int i=0;i<3;++i)
        result[i]/=resultLen;
    glNormal(result);
    for(int i=0;i<3;++i)
        result[i]*=radius;
    glVertex(result);
}

inline void combine(const GLVector<GLdouble,3>& p00,const GLVector<GLdouble,3>& p10,const GLVector<GLdouble,3>& p01,const GLVector<GLdouble,3>& p11,double wx,double wy,double radius)
{
    GLVector<GLdouble,3> result;
    double resultLen=0.0;
    for(int i=0;i<3;++i)
    {
        double bot=p00[i]*(1.0-wx)+p10[i]*wx;
        double top=p01[i]*(1.0-wx)+p11[i]*wx;
        result[i]=bot*(1.0-wy)+top*wy;
        resultLen+=Math::sqr(result[i]);
    }
    resultLen=Math::sqrt(resultLen);
    for(int i=0;i<3;++i)
        result[i]/=resultLen;
    glNormal(result);
    for(int i=0;i<3;++i)
        result[i]*=radius;
    glVertex(result);
}

void drawSphere(double radius,int subdivision)
{
    typedef GLVector<GLdouble,3> Vector;
    
    const double bX=0.525731112119133606;
    const double bZ=0.850650808352039932;
    static Vector vUnit[12]={Vector(-bX,0.0, bZ),Vector( bX,0.0, bZ),Vector(-bX,0.0,-bZ),Vector( bX,0.0,-bZ),
                             Vector(0.0, bZ, bX),Vector(0.0, bZ,-bX),Vector(0.0,-bZ, bX),Vector(0.0,-bZ,-bX),
                             Vector( bZ, bX,0.0),Vector(-bZ, bX,0.0),Vector( bZ,-bX,0.0),Vector(-bZ,-bX,0.0)};
    static const int stripIndices[12]={0,1,4,8,5,3,2,7,11,6,0,1};
    static const int fanIndices[2][7]={{9,0,4,5,2,11,0},{10,1,6,7,3,8,1}};
    
    /* Render the central triangle strips: */
    for(int strip=0;strip<subdivision;++strip)
    {
        double botW=double(strip)/double(subdivision);
        double topW=double(strip+1)/double(subdivision);
        glBegin(GL_TRIANGLE_STRIP);
        for(int i=0;i<10;i+=2)
        {
            Vector p00=vUnit[stripIndices[i+1]];
            Vector p10=vUnit[stripIndices[i+3]];
            Vector p01=vUnit[stripIndices[i+0]];
            Vector p11=vUnit[stripIndices[i+2]];
            for(int j=0;j<subdivision;++j)
            {
                double leftW=double(j)/double(subdivision);
                // double rightW=double(j+1)/double(subdivision);
                combine(p00,p10,p01,p11,leftW,topW,radius);
                combine(p00,p10,p01,p11,leftW,botW,radius);
            }
            combine(p00,p10,p01,p11,1.0,topW,radius);
            combine(p00,p10,p01,p11,1.0,botW,radius);
        }
        glEnd();
    }
    
    for(int cap=0;cap<2;++cap)
    {
        /* Render the cap triangle strips: */
        for(int strip=0;strip<subdivision-1;++strip)
        {
            double botW=double(strip)/double(subdivision);
            double topW=double(strip+1)/double(subdivision);
            glBegin(GL_TRIANGLE_STRIP);
            combine(vUnit[fanIndices[cap][0]],vUnit[fanIndices[cap][2]],vUnit[fanIndices[cap][1]],topW,0.0,radius);
            for(int i=1;i<6;++i)
            {
                Vector p100=vUnit[fanIndices[cap][0]];
                Vector p010=vUnit[fanIndices[cap][i]];
                Vector p001=vUnit[fanIndices[cap][i+1]];
                for(int j=0;j<subdivision-strip;++j)
                {
                    double leftW=double(j)/double(subdivision);
                    combine(p100,p001,p010,botW,leftW,radius);
                    combine(p100,p001,p010,topW,leftW,radius);
                }
            }
            combine(vUnit[fanIndices[cap][0]],vUnit[fanIndices[cap][2]],vUnit[fanIndices[cap][1]],botW,0.0,radius);
            glEnd();
        }
        
        /* Render the cap triangle fan: */
        glBegin(GL_TRIANGLE_FAN);
        combine(vUnit[fanIndices[cap][0]],vUnit[fanIndices[cap][2]],vUnit[fanIndices[cap][1]],1.0,0.0,radius);
        double botW=double(subdivision-1)/double(subdivision);
        for(int i=1;i<6;++i)
            combine(vUnit[fanIndices[cap][0]],vUnit[fanIndices[cap][i+1]],vUnit[fanIndices[cap][i]],botW,0.0,radius);
        combine(vUnit[fanIndices[cap][0]],vUnit[fanIndices[cap][2]],vUnit[fanIndices[cap][1]],botW,0.0,radius);
        glEnd();
    }
}
