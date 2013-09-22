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

#ifndef IK_INCLUDED
#define IK_INCLUDED

#include <Geometry/Vector.h>
#include <Geometry/Point.h>
#include <Geometry/OrthonormalTransformation.h>
#include <Geometry/AffineTransformation.h>
#include <DenseMatrix.h>

class IK
	{
	/* Embedded classes: */
	public:
	typedef double Scalar;
	typedef Geometry::Vector<Scalar,3> Vector;
	typedef Geometry::Point<Scalar,3> Point;
	typedef Geometry::OrthonormalTransformation<Scalar,3> Transformation;
	typedef Geometry::AffineTransformation<Scalar,3> AffineTransformation;
	
	private:
	struct Joint // Structure for joints with a single axis of unconstrained rotation
		{
		/* Elements: */
		public:
		Point offset; // Origin of joint's axis of rotation
		Vector axis; // Joint's axis of rotation (normalized)
		Scalar angle; // Joint's current rotation angle (in radians)
		
		/* Methods: */
		Transformation calcJointTransformation(void) const; // Calculates transformation between this and next joint in chain
		};
	
	/* Elements: */
	int numJoints; // Number of joints in chain
	Joint* joints; // Array of joints in chain
	Point leafOffset; // Origin of leaf joint's end effector coordinate system
	
	/* IK algorithm parameters: */
	Scalar orientationWeights[3]; // Weights balancing rotational with translational approximation
	Scalar maxResidual; // Residue threshold at which iteration is terminated prematurely
	Scalar maxStepError; // Error threshold at which a step is undone
	int effectorJointIndex; // Index of effector joint (equal to numJoints for IK on entire chain)
	
	/* Transient IK state: */
	Transformation* jointTransformations; // Array of accumulated transformations from the root to the leaf joint
	DenseMatrix jacobian; // Jacobian matrix of the joint chain with respect to the goal transformation
	DenseMatrix deltaT; // Difference vector between current end effector transformation and goal transformation
	DenseMatrix deltaP; // Vector of parameter increments for each joint in the chain
	Scalar stepSize; // Step size of last IK iteration
	
	/* Extra bookkeeping for algorithm analysis purposes: */
	int iterationArraySize; // Size of allocated iteration bookkeeping arrays
	Scalar* iterationResiduals; // Array of residual values from each iteration step
	Scalar* iterationStepSizes; // Array of step sizes for each iteration step
	Scalar* iterationErrors; // Array of step error estimates computed in each step by adaptive methods
	int numIterationSteps; // Number of iteration steps performed before termination (early or scheduled)
	
	/* Private methods: */
	void calcJointTransformations(void); // Accumulates global joint transformations along the chain
	void calcJacobian(void); // Calculates the Jacobian matrix of the entire chain
	Scalar calcDeltaP(void); // Calculates parameter update vector based on the current chain state and the current end effector transformation; returns residual
	void updateJointAngles(Scalar currentStepSize); // Applies the computed parameter update vector to the joints' angles
	Scalar performFirstOrderStep(const Transformation& goalTransformation,Scalar currentStepSize); // Performs a first-order iteration step towards the given goal transformation with the given step size; returns residual
	Scalar performSecondOrderStep(const Transformation& goalTransformation,Scalar currentStepSize); // Performs a second-order iteration step towards the given goal transformation with the given step size; returns residual
	
	/* Bookkeeping methods: */
	void initBookkeeping(int maxNumSteps);
	
	/* Constructors and destructors: */
	public:
	IK(int sNumJoints); // Creates uninitialized IK chain of given number of sensors
	~IK(void);
	
	/* Methods: */
	void setJoint(int jointIndex,const Point& newOffset,const Vector& newAxis); // Sets joint's offset and axis
	void setJointAngle(int jointIndex,Scalar newAngle); // Sets joint's current angle of rotation (wraps to 0...2*pi)
	void setLeafOffset(const Point& newLeafOffset); // Sets the leaf offset in leaf joint's coordinate system
	int getNumJoints(void) const // Returns number of joints in the IK chain
		{
		return numJoints;
		};
	const Point& getJointOffset(int jointIndex) const // Returns joint's offset in local joint coordinates
		{
		return joints[jointIndex].offset;
		};
	const Vector& getJointAxis(int jointIndex) const // Returns joint's rotation axis in local joint coordinates
		{
		return joints[jointIndex].axis;
		};
	const Point& getLeafOffset(void) const // Returns end effector's offset in last joint's local coordinates
		{
		return leafOffset;
		};
	int getEffectorJointIndex(void) const // Returns index of current effector joint
		{
		return effectorJointIndex;
		};
	void setOrientationWeights(const Scalar newOrientationWeights[3]); // Sets orientation weight parameters for IK
	void setMaxResidual(Scalar newMaxResidual); // Sets maximum residual threshold for IK
	void setMaxStepError(Scalar newMaxStepError); // Sets maximum error threshold for adaptive IK
	void setEffectorJointIndex(int newEffectorJointIndex); // Sets the index of the effector joint (for constrained IK)
	void setStepSize(Scalar newStepSize); // Sets initial step size for IK iteration
	Scalar performSteps(const Transformation& goalTransformation,int maxNumSteps); // Performs a set of IK iterations; returns last residual value
	Scalar getStepSize(void) const // Returns current iteration step size
		{
		return stepSize;
		};
	Transformation getJointTransformation(int jointIndex) const; // Returns global transformation of given joint in root joint's coordinate system
	Scalar getJointAngle(int jointIndex) const // Returns the current rotation angle of a joint
		{
		return joints[jointIndex].angle;
		};
	Transformation getEffectorTransformation(void) const; // Returns global transformation of end effector in root joint's coordinate system
	
	/* Bookkeeping access methods: */
	int getNumIterationSteps(void) const	
		{
		return numIterationSteps;
		};
	Scalar getIterationResidual(int stepIndex) const
		{
		return iterationResiduals[stepIndex];
		};
	Scalar getIterationStepSize(int stepIndex) const
		{
		return iterationStepSizes[stepIndex];
		};
	Scalar getIterationError(int stepIndex) const
		{
		return iterationErrors[stepIndex];
		};
	};

#endif
