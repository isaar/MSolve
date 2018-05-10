using ISAAR.MSolve.Numerical.Interfaces;
using System;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    /**
     * Class for quaternions to represent finite rotations.
     *
     * @author Theofilos Manitaras
     */
    public class Quaternion
    {

        private static int VECTOR_COMPONENT_COUNT = 3;

        /**
	     * Creates a quaternion from individual scalar and vector parts.
	     *
	     * @param scalarPart
	     *            The scalar part
	     * @param vectorPart
	     *            The vector part
	     * @return The new quaternion
	     */
        public static Quaternion CreateFromIndividualParts(double scalarPart, Vector vectorPart)
        {
            return new Quaternion(scalarPart, vectorPart);
        }

        /**
	     * Creates a new quaternion of zero angle.
	     *
	     * @return The new quaternion
	     */
        public static Quaternion OfZeroAngle()
        {
            return new Quaternion(1.0, new Vector(VECTOR_COMPONENT_COUNT));
        }

        private double scalarPart;
        private readonly Vector vectorPart;

        public Quaternion(double scalarPart, Vector vectorPart)
        {
            this.scalarPart = scalarPart;
            this.vectorPart = vectorPart;
        }

        /**
	     * Rotates the quaternion by an incremental rotation given in vector form.
	     *
	     * @param incrementalRotation
	     *            The incremental rotation
	     */
        public void ApplyIncrementalRotation(Vector incrementalRotation)
        {
            var vectorPartIncrement = new Vector(VECTOR_COMPONENT_COUNT);
            incrementalRotation.CopyTo(vectorPartIncrement.Data, 0);
            vectorPartIncrement.Scale(0.5);
            double squareRootArgument = 1.0 - vectorPartIncrement * vectorPartIncrement;
            //Preconditions.checkArgument(squareRootArgument > 0.0, "Very large rotation increment applied");
            double scalarPartIncrement = Math.Sqrt(squareRootArgument);
            double updatedScalarPart = (scalarPartIncrement * this.scalarPart) - vectorPartIncrement * (this.vectorPart);
            Vector updatedVectorPart = this.vectorPart;
            updatedVectorPart.Scale(scalarPartIncrement);
            updatedVectorPart.Add(scalarPart * vectorPartIncrement);
            Vector incrementCrossVectorPart = vectorPartIncrement ^ this.vectorPart;
            updatedVectorPart.Add(incrementCrossVectorPart);
            this.scalarPart = updatedScalarPart;
            updatedVectorPart.CopyTo(this.vectorPart.Data, 0);
        }

        /**
	     * Rotates the quaternion to a new quaternion by an incremental rotation given in vector form.
	     *
	     * @param incrementalRotation
	     *            The incremental rotation
	     * @return The rotated quaternion
	     */
        public Quaternion ApplyIncrementalRotationToNew(Vector incrementalRotation)
        {
            var vectorPartIncrement = new Vector(VECTOR_COMPONENT_COUNT);
            incrementalRotation.CopyTo(vectorPartIncrement.Data, 0);
            vectorPartIncrement.Scale(0.5);
            double squareRootArgument = 1.0 - vectorPartIncrement * vectorPartIncrement;
            //Preconditions.checkArgument(squareRootArgument > 0.0, "Very large rotation increment applied");
            double scalarPartIncrement = Math.Sqrt(squareRootArgument);
            double updatedScalarPart = (scalarPartIncrement * this.scalarPart) - vectorPartIncrement * (this.vectorPart);
            Vector updatedVectorPart = this.vectorPart;
            updatedVectorPart.Scale(scalarPartIncrement);
            updatedVectorPart.Add(scalarPart * vectorPartIncrement);
            Vector incrementCrossVectorPart = vectorPartIncrement ^ this.vectorPart;
            updatedVectorPart.Add(incrementCrossVectorPart);
            return new Quaternion(updatedScalarPart, updatedVectorPart);
        }

        public Matrix2D GetRotationMatrix()
        {
            Matrix2D rotationMatrix = new Matrix2D(VECTOR_COMPONENT_COUNT, VECTOR_COMPONENT_COUNT);
            double r0 = this.scalarPart;
            double r1 = this.vectorPart[0];
            double r2 = this.vectorPart[1];
            double r3 = this.vectorPart[2];
            double r0Squared = r0 * r0;
            double r1Squared = r1 * r1;
            double r2Squared = r2 * r2;
            double r3Squared = r3 * r3;

            rotationMatrix[0, 0] = (r0Squared + r1Squared) - r2Squared - r3Squared;
            rotationMatrix[0, 1] = 2.0 * ((r1 * r2) - (r0 * r3));
            rotationMatrix[0, 2] = 2.0 * ((r1 * r3) + (r0 * r2));
            rotationMatrix[1, 0] = 2.0 * ((r2 * r1) + (r0 * r3));
            rotationMatrix[1, 1] = (r0Squared + r2Squared) - r1Squared - r3Squared;
            rotationMatrix[1, 2] = 2.0 * ((r2 * r3) - (r0 * r1));
            rotationMatrix[2, 0] = 2.0 * ((r3 * r1) - (r0 * r2));
            rotationMatrix[2, 1] = 2.0 * ((r3 * r2) + (r0 * r1));
            rotationMatrix[2, 2] = (r0Squared + r3Squared) - r1Squared - r2Squared;

            return rotationMatrix;
        }

        public double ScalarPart { get { return this.scalarPart; } }
        public Vector VectorPart { get { return this.vectorPart; } }
    }
}
