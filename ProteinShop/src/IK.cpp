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
IK - Class to perform inverse kinematics on non-branching chains of
single-axis rotational joints.
***********************************************************************/

#include <Math/Math.h>
#include <Math/Constants.h>

#include "IK.h"

#define IK_FIRSTORDER 1
#define IK_FIRSTORDER_ADAPTIVE 2
#define IK_SECONDORDER 3
#define IK_SECONDORDER_ADAPTIVE 4

#define IK_METHOD IK_FIRSTORDER_ADAPTIVE

/**************************
Methods of class IK::Joint:
**************************/

IK::Transformation IK::Joint::calcJointTransformation(void) const
	{
	return Transformation(offset-Point::origin,Transformation::Rotation::rotateAxis(axis,angle));
	}

/*******************
Methods of class IK:
*******************/

void IK::calcJointTransformations(void)
	{
	/* Accumulate local transformations along the chain: */
	for(int i=0;i<numJoints;++i)
		{
		/* Calculate the current joint's incremental transformation: */
		jointTransformations[i+1]=joints[i].calcJointTransformation();

		/* Accumulate global joint transformations: */
		jointTransformations[i+1].leftMultiply(jointTransformations[i]);
		}
	}

void IK::calcJacobian(void)
	{
	/* Convert the end effector's transformation into an affine matrix: */
	AffineTransformation effectorTransformation=jointTransformations[numJoints];
	
	/* Extract origin and axes from the effector transformation: */
	Point effectorOrigin=effectorTransformation.transform(leafOffset);
	Vector effectorX=effectorTransformation.getDirection(0);
	Vector effectorY=effectorTransformation.getDirection(1);
	Vector effectorZ=effectorTransformation.getDirection(2);
	
	/*********************************************************************
	Calculate the Jacobian matrix.
	*********************************************************************/
	
	/* Calculate the Jacobian matrix column for each joint: */
	for(int columnIndex=0;columnIndex<numJoints;++columnIndex)
		{
		/* Calculate joint offset and axis vector in world coordinates: */
		Point gOffset=jointTransformations[columnIndex].transform(joints[columnIndex].offset);
		Vector gAxis=jointTransformations[columnIndex].transform(joints[columnIndex].axis);
		
		/*******************************************************************
		Calculate translational Jacobian entries.
		*******************************************************************/
		
		/* Calculate dOrigin/dAngle: */
		Vector dOrigin=Geometry::cross(gAxis,effectorOrigin-gOffset);
		for(int i=0;i<3;++i)
			jacobian(i,columnIndex)=dOrigin[i];
		
		/*******************************************************************
		Calculate rotational Jacobian entries.
		*******************************************************************/
		
		/* Calculate d(X*goalY)/dAngle: */
		jacobian(3,columnIndex)=(gAxis[2]*effectorX[0]-gAxis[0]*effectorX[2]); // *orientationWeights[0];
		
		/* Calculate d(Y*goalZ)/dAngle: */
		jacobian(4,columnIndex)=(gAxis[0]*effectorY[1]-gAxis[1]*effectorY[0]); // *orientationWeights[1];
		
		/* Calculate d(Z*goalX)/dAngle: */
		jacobian(5,columnIndex)=(gAxis[1]*effectorZ[2]-gAxis[2]*effectorZ[1]); // *orientationWeights[2];
		}
	}

IK::Scalar IK::calcDeltaP(void)
	{
	/* Convert the end effector's transformation into an affine matrix: */
	AffineTransformation effectorTransformation=jointTransformations[effectorJointIndex];
	
	/* Extract origin and axes from the effector transformation: */
	Point effectorOrigin=effectorJointIndex==numJoints?effectorTransformation.transform(leafOffset):effectorTransformation.getOrigin();
	Vector effectorX=effectorTransformation.getDirection(0);
	Vector effectorY=effectorTransformation.getDirection(1);
	Vector effectorZ=effectorTransformation.getDirection(2);
	
	/*********************************************************************
	Calculate the goal difference vector and the current residual.
	*********************************************************************/
	
	/* Calculate the goal difference vector: */
	Scalar deltaTv[6];
	for(int i=0;i<3;++i)
		deltaT(i,0)=deltaTv[i]=-effectorOrigin[i];
	deltaT(3,0)=deltaTv[3]=-effectorX[1]*orientationWeights[0];
	deltaT(4,0)=deltaTv[4]=-effectorY[2]*orientationWeights[1];
	deltaT(5,0)=deltaTv[5]=-effectorZ[0]*orientationWeights[2];
	
	/* Calculate the residual: */
	Scalar residual=Scalar(0);
	for(int i=0;i<6;++i)
		residual+=Math::sqr(deltaTv[i]);
	residual=Math::sqrt(residual);
	
	/* Multiply the private copy of the goal difference vector with the orientation weights again: */
	for(int i=0;i<3;++i)
		deltaTv[3+i]*=orientationWeights[i];
	
	/*********************************************************************
	Calculate parameter increments for each joint by calculating each
	joint's Jacobian entries, and multiplying them with the parameter
	difference vector.
	*********************************************************************/
	
	/* Calculate the Jacobian matrix column for each joint: */
	for(int columnIndex=0;columnIndex<effectorJointIndex;++columnIndex)
		{
		/* We calculate the Jacobian matrix column implicitly as well: */
		Scalar dAngle=Scalar(0);
		
		/* Calculate joint offset and axis vector in world coordinates: */
		Point gOffset=jointTransformations[columnIndex].transform(joints[columnIndex].offset);
		Vector gAxis=jointTransformations[columnIndex].transform(joints[columnIndex].axis);
		
		/*******************************************************************
		Calculate translational Jacobian entries.
		*******************************************************************/
		
		/* Calculate dOrigin/dAngle: */
		Vector dOrigin=Geometry::cross(gAxis,effectorOrigin-gOffset);
		for(int i=0;i<3;++i)
			dAngle+=dOrigin[i]*deltaTv[i];
		
		/*******************************************************************
		Calculate rotational Jacobian entries.
		
		Conceptually, the rotational Jacobian entries have to be multiplied
		with the orientation weights here, since the goal difference vector
		is multiplied with them. As a simple optimization, we multiply the
		goal difference vector with the weights twice.
		*******************************************************************/
		
		/* Calculate d(X*goalY)/dAngle: */
		dAngle+=(gAxis[2]*effectorX[0]-gAxis[0]*effectorX[2])*deltaTv[3];
		
		/* Calculate d(Y*goalZ)/dAngle: */
		dAngle+=(gAxis[0]*effectorY[1]-gAxis[1]*effectorY[0])*deltaTv[4];
		
		/* Calculate d(Z*goalX)/dAngle: */
		dAngle+=(gAxis[1]*effectorZ[2]-gAxis[2]*effectorZ[1])*deltaTv[5];
		
		deltaP(columnIndex,0)=dAngle;
		}
	
	return residual;
	}

void IK::updateJointAngles(IK::Scalar currentStepSize)
	{
	const Scalar twoPi=Scalar(2)*Math::Constants<Scalar>::pi;
	
	/* Update all joint angles by parameter update vector times step size: */
	for(int i=0;i<numJoints;++i)
		{
		/* Change parameter by update vector multiplied by step size: */
		Scalar jointWeight=Scalar(i+1)/Scalar(numJoints);
		Scalar newAngle=joints[i].angle+deltaP(i,0)*currentStepSize;

		/* Wrap angle to 0...2*pi: */
		newAngle-=Math::floor(newAngle/twoPi)*twoPi;

		/* Set joint angle: */
		joints[i].angle=newAngle;
		}
	}

IK::Scalar IK::performFirstOrderStep(const Transformation& goalTransformation,IK::Scalar currentStepSize)
	{
	/*******************************************************************
	Calculate chain of current global transformations from the root to
	the leaf. The transformations are set up such that the leaf's
	transformation when reaching the goal is the identity
	transformation.
	*******************************************************************/
	
	/* Calculate joint transformations: */
	calcJointTransformations();
	
	/*******************************************************************
	Calculate the joint parameter update vector and the current
	residual.
	*******************************************************************/
	
	Scalar residual=calcDeltaP();
	
	/*******************************************************************
	If the effector joint is not the end effector (bi-directional IK),
	the calculated parameter update vector must be modified to not move
	the end effector. More precisely, the parameter update vector must
	be projected into the full chain's Jacobian matrix' null space.
	*******************************************************************/
	
	if(effectorJointIndex<numJoints)
		{
		/*******************************************************************
		Initialize the parameter update values for joints after the effector
		to zero.
		*******************************************************************/
		
		for(int i=effectorJointIndex;i<numJoints;++i)
			deltaP(i,0)=Scalar(0);
		
		/*******************************************************************
		Calculate the full chain's Jacobian matrix.
		*******************************************************************/
		
		calcJacobian();
		
		/*******************************************************************
		Project the calculated parameter update vector into the Jacobian's
		null space. If IK were linear, this would avoid any motion of the
		end effector.
		*******************************************************************/
		
		/* Calculate J * J^T: */
		DenseMatrix jjt(6,6);
		inplaceTransposed2Multiplication(jjt,jacobian,jacobian);
		
		/* Calculate the parameter update vector's effect on the end effector: */
		DenseMatrix jpv(6,1);
		inplaceMultiplication(jpv,jacobian,deltaP);
		
		/* First projection step: */
		DenseMatrix p1=jjt.solveLinearEquations(jpv);
		
		/* Second projection step: */
		DenseMatrix pvn(numJoints,1);
		inplaceTransposed1Multiplication(pvn,jacobian,p1);
		
		/* Subtract pvn from deltaP to project deltaP into Jacobian's null space: */
		for(int i=0;i<numJoints;++i)
			deltaP(i,0)-=pvn(i,0);
		}
	
	/*******************************************************************
	Perform a first-order integration step towards the goal using the
	parameter update vector.
	*******************************************************************/
	
	updateJointAngles(currentStepSize);
	
	return residual;
	}

IK::Scalar IK::performSecondOrderStep(const IK::Transformation& goalTransformation,IK::Scalar currentStepSize)
	{
	/*******************************************************************
	Save the initial state.
	*******************************************************************/
	
	Scalar* initialJointParameters=new Scalar[numJoints];
	for(int i=0;i<numJoints;++i)
		initialJointParameters[i]=joints[i].angle;
	
	/*******************************************************************
	Perform a first-order step from the initial state to the midpoint of
	the current step size.
	*******************************************************************/
	
	/* Calculate joint transformations: */
	calcJointTransformations();
	
	/* Calculate the parameter update vector and the residual: */
	Scalar residual=calcDeltaP();
	
	/* Move by half the step size: */
	updateJointAngles(currentStepSize*Scalar(0.5));
	
	/*******************************************************************
	Calculate the parameter update vector at the midpoint state.
	*******************************************************************/
	
	/* Calculate joint transformations: */
	calcJointTransformations();
	
	/* Calculate the parameter update vector: */
	calcDeltaP();
	
	/*******************************************************************
	Apply the midpoint parameter update vector to the initial state.
	*******************************************************************/
	
	/* Go back to the initial state: */
	for(int i=0;i<numJoints;++i)
		joints[i].angle=initialJointParameters[i];
	
	/* Move by the full step size: */
	updateJointAngles(currentStepSize);
	
	/* Clean up: */
	delete[] initialJointParameters;
	
	return residual;
	}

void IK::initBookkeeping(int maxNumSteps)
	{
	/* Update bookkeeping arrays: */
	if(maxNumSteps>iterationArraySize)
		{
		delete[] iterationResiduals;
		delete[] iterationStepSizes;
		delete[] iterationErrors;
		
		iterationArraySize=(maxNumSteps*3)/2;
		
		iterationResiduals=new Scalar[iterationArraySize];
		iterationStepSizes=new Scalar[iterationArraySize];
		iterationErrors=new Scalar[iterationArraySize];
		}
	}

IK::IK(int sNumJoints)
	:numJoints(sNumJoints),
	 joints(new Joint[numJoints]),
	 effectorJointIndex(numJoints),
	 jointTransformations(new Transformation[numJoints+1]),
	 jacobian(6,numJoints),
	 deltaT(6,1),
	 deltaP(numJoints,1)
	{
	/* Initialize bookkeeping arrays: */
	iterationArraySize=0;
	iterationResiduals=0;
	iterationStepSizes=0;
	iterationErrors=0;
	}

IK::~IK(void)
	{
	delete[] joints;
	delete[] jointTransformations;
	
	/* Delete bookkeeping arrays: */
	delete[] iterationResiduals;
	delete[] iterationStepSizes;
	delete[] iterationErrors;
	}

void IK::setJoint(int jointIndex,const IK::Point& newOffset,const IK::Vector& newAxis)
	{
	/* Set joint offset: */
	joints[jointIndex].offset=newOffset;
	
	/* Set and normalize joint axis: */
	joints[jointIndex].axis=newAxis;
	joints[jointIndex].axis.normalize();
	}

void IK::setJointAngle(int jointIndex,IK::Scalar newAngle)
	{
	/* Wrap new angle to 0...2*pi: */
	const Scalar twoPi=Scalar(2)*Math::Constants<Scalar>::pi;
	newAngle-=Math::floor(newAngle/twoPi)*twoPi;
	
	/* Set joint angle: */
	joints[jointIndex].angle=newAngle;
	}

void IK::setLeafOffset(const IK::Point& newLeafOffset)
	{
	/* Set leaf offset: */
	leafOffset=newLeafOffset;
	}

void IK::setOrientationWeights(const IK::Scalar newOrientationWeights[3])
	{
	/* Set orientation weights: */
	for(int i=0;i<3;++i)
		orientationWeights[i]=newOrientationWeights[i];
	}

void IK::setMaxResidual(IK::Scalar newMaxResidual)
	{
	/* Set max residual: */
	maxResidual=newMaxResidual;
	}

void IK::setMaxStepError(IK::Scalar newMaxStepError)
	{
	/* Set max step error: */
	maxStepError=newMaxStepError;
	}

void IK::setEffectorJointIndex(int newEffectorJointIndex)
	{
	/* Set effector joint index: */
	effectorJointIndex=newEffectorJointIndex;
	}

void IK::setStepSize(IK::Scalar newStepSize)
	{
	/* Set initial step size: */
	stepSize=newStepSize;
	}

#if IK_METHOD == IK_FIRSTORDER

/***********************************************************************
IK iteration method using first-order steps with fixed step size.
***********************************************************************/

IK::Scalar IK::performSteps(const IK::Transformation& goalTransformation,int maxNumSteps)
	{
	/*********************************************************************
	Initialize the IK object for iteration on the given goal
	transformation.
	*********************************************************************/
	
	/* Initialize bookkeeping arrays: */
	initBookkeeping(maxNumSteps);
	
	/* Initialize transformations with the inverse of the goal transformation: */
	jointTransformations[0]=Geometry::invert(goalTransformation);
	
	/*********************************************************************
	Perform IK steps until the maximum number of steps is reached, or the
	current residual becomes smaller than the preset threshold.
	*********************************************************************/
	
	Scalar residual;
	int step;
	for(step=0;step<maxNumSteps;++step)
		{
		/* Save step size: */
		iterationStepSizes[step]=stepSize;
		
		/*******************************************************************
		Perform a first-order iteration step.
		*******************************************************************/
		
		/* Perform the step and save the computed residual: */
		residual=performFirstOrderStep(goalTransformation,stepSize);
		
		/* Save residual: */
		iterationResiduals[step]=residual;
		
		/* Set error value to default (not computed here): */
		iterationErrors[step]=Scalar(0);
		
		/* Stop iteration if residual smaller than threshold: */
		if(residual<maxResidual)
			break;
		}
	
	/* Update bookkeeping arrays: */
	numIterationSteps=step;
	
	/* Return the final residual: */
	return residual;
	}

#endif

#if IK_METHOD == IK_FIRSTORDER_ADAPTIVE

/***********************************************************************
IK iteration method using first-order steps with adaptive step size.
***********************************************************************/

IK::Scalar IK::performSteps(const IK::Transformation& goalTransformation,int maxNumSteps)
	{
	/*********************************************************************
	Initialize the IK object for iteration on the given goal
	transformation.
	*********************************************************************/
	
	/* Initialize bookkeeping arrays: */
	initBookkeeping(maxNumSteps);
	
	/* Initialize transformations with the inverse of the goal transformation: */
	jointTransformations[0]=Geometry::invert(goalTransformation);
	
	/* Create step undo array: */
	Scalar* oldJointAngles=new Scalar[numJoints];
	
	/* Create parameter increment difference vector: */
	Scalar* deltaDeltaP=new Scalar[numJoints];
	
	/*********************************************************************
	Perform IK steps until the maximum number of steps is reached, or the
	current residual becomes smaller than the preset threshold.
	*********************************************************************/
	
	Scalar residual;
	int step;
	for(step=0;step<maxNumSteps;++step)
		{
		/* Save step size: */
		iterationStepSizes[step]=stepSize;
		
		/*******************************************************************
		Save joint rotation angles to be able to undo a step with too large
		an error.
		*******************************************************************/
		
		for(int i=0;i<numJoints;++i)
			oldJointAngles[i]=joints[i].angle;
		
		/*******************************************************************
		Perform the first first-order iteration step with half the current
		step size.
		*******************************************************************/
		
		residual=performFirstOrderStep(goalTransformation,stepSize*Scalar(0.5f));
		
		/* Save residual: */
		iterationResiduals[step]=residual;
		
		/* Save the parameter update vector for comparison with the second: */
		for(int i=0;i<numJoints;++i)
			deltaDeltaP[i]=deltaP(i,0);
		
		/* Stop iteration if residual is smaller than threshold: */
		if(residual<maxResidual)
			break;
		
		/*******************************************************************
		Perform the second first-order iteration step with half the current
		step size.
		*******************************************************************/
		
		performFirstOrderStep(goalTransformation,stepSize*Scalar(0.5f));
		
		/* Calculate the difference between the two deltaP vectors: */
		for(int i=0;i<numJoints;++i)
			deltaDeltaP[i]-=deltaP(i,0);
		
		/*******************************************************************
		Estimate the current step's error to adjust the step size
		accordingly.
		*******************************************************************/
		
		/* Calculate this step's error: */
		Scalar stepError=Scalar(0);
		for(int i=0;i<numJoints;++i)
			{
			Scalar e=Math::abs(deltaDeltaP[i]);
			if(stepError<e)
				stepError=e;
			}
		stepError*=(stepSize*Scalar(0.5))/residual;
		
		/* Save error value: */
		iterationErrors[step]=stepError;
		
		#if 1
		/* Adapt step size based on step error: */
		if(stepError>maxStepError)
			{
			/*****************************************************************
			Undo this iteration step and decrease the step size for the next
			step.
			*****************************************************************/
			
			#if 1
			/* Reverse the step: */
			for(int i=0;i<numJoints;++i)
				joints[i].angle=oldJointAngles[i];
			#endif
			
			/* Decrease the step size: */
			stepSize*=Scalar(0.98);
			}
		else
			{
			/*****************************************************************
			Increase the step size slightly for the next step.
			*****************************************************************/
			
			stepSize*=Scalar(1.001);
			}
		#endif
		}
	
	/* Update bookkeeping arrays: */
	numIterationSteps=step;
	
	/* Clean up: */
	delete[] oldJointAngles;
	delete[] deltaDeltaP;
	
	/* Return final residual: */
	return residual;
	}

#endif

IK::Transformation IK::getJointTransformation(int jointIndex) const
	{
	/* Accumulate transformations starting at the root of the chain: */
	Transformation result=Transformation::identity;
	for(int i=0;i<jointIndex;++i)
		result*=joints[i].calcJointTransformation();
	
	return result;
	}

IK::Transformation IK::getEffectorTransformation(void) const
	{
	/* Accumulate transformations starting at the root of the chain: */
	Transformation result=Transformation::identity;
	for(int i=0;i<effectorJointIndex;++i)
		result*=joints[i].calcJointTransformation();
	
	if(effectorJointIndex==numJoints)
		{
		/* Tag on the transformation from the leaf joint to the end effector: */
		result*=Transformation::translate(leafOffset-Point::origin);
		}
	
	return result;
	}
