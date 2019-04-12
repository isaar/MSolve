using System;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Geometry
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
        public static Quaternion CreateFromIndividualParts(double scalarPart, double[] vectorPart)
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
            return new Quaternion(1.0, new double[VECTOR_COMPONENT_COUNT]);
        }

        private double scalarPart;
        private readonly double[] vectorPart;

        public Quaternion(double scalarPart, double[] vectorPart)
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
        public void ApplyIncrementalRotation(double[] incrementalRotation)
        {
            Debug.Assert(incrementalRotation.Length == VECTOR_COMPONENT_COUNT);
            double[] vectorPartIncrement = incrementalRotation.Copy();
            //var vectorPartIncrement = new Vector(VECTOR_COMPONENT_COUNT);
            //incrementalRotation.CopyTo(vectorPartIncrement.Data, 0);
            vectorPartIncrement.ScaleIntoThis(0.5);
            double squareRootArgument = 1.0 - vectorPartIncrement.DotProduct(vectorPartIncrement);
            //Preconditions.checkArgument(squareRootArgument > 0.0, "Very large rotation increment applied");
            double scalarPartIncrement = Math.Sqrt(squareRootArgument);
            double updatedScalarPart = (scalarPartIncrement * this.scalarPart) - vectorPartIncrement.DotProduct(this.vectorPart);

            //TODO: Is the following correct? Vector updatedVectorPart = this.vectorPart; does not copy the original vector.
            //TODO: Why not use this.vectorPart for all the next? 
            double[] updatedVectorPart = this.vectorPart;
            updatedVectorPart.ScaleIntoThis(scalarPartIncrement);
            updatedVectorPart.AddIntoThis(vectorPartIncrement.Scale(scalarPart));
            double[] incrementCrossVectorPart = vectorPartIncrement.CrossProduct(this.vectorPart);
            updatedVectorPart.AddIntoThis(incrementCrossVectorPart);
            this.scalarPart = updatedScalarPart;
            this.vectorPart.CopyFrom(updatedVectorPart); //TODO: This does nothing. updatedVectorPart is a reference to the same memory as this.vectorPart
        }

        /**
	     * Rotates the quaternion to a new quaternion by an incremental rotation given in vector form.
	     *
	     * @param incrementalRotation
	     *            The incremental rotation
	     * @return The rotated quaternion
	     */
        public Quaternion ApplyIncrementalRotationToNew(double[] incrementalRotation)
        {
            Debug.Assert(incrementalRotation.Length == VECTOR_COMPONENT_COUNT);
            double[] vectorPartIncrement = incrementalRotation.Copy();
            vectorPartIncrement.ScaleIntoThis(0.5);
            double squareRootArgument = 1.0 - vectorPartIncrement.DotProduct(vectorPartIncrement);
            //Preconditions.checkArgument(squareRootArgument > 0.0, "Very large rotation increment applied");
            double scalarPartIncrement = Math.Sqrt(squareRootArgument);
            double updatedScalarPart = (scalarPartIncrement * this.scalarPart) - vectorPartIncrement.DotProduct(this.vectorPart);

            //TODO: Is the following correct? Vector updatedVectorPart = this.vectorPart; does not copy the original vector.
            //TODO: The next part creates a new Quaternion (using legacy linear algebra classes), but the existing ones is not copied as per the author's intention
            double[] updatedVectorPart = this.vectorPart;
            updatedVectorPart.ScaleIntoThis(scalarPartIncrement);
            updatedVectorPart.AddIntoThis(vectorPartIncrement.Scale(scalarPart));
            double[] incrementCrossVectorPart = vectorPartIncrement.CrossProduct(this.vectorPart);
            updatedVectorPart.AddIntoThis(incrementCrossVectorPart);
            return new Quaternion(updatedScalarPart, updatedVectorPart);
        }

        public Matrix GetRotationMatrix()
        {
            var rotationMatrix = Matrix.CreateZero(VECTOR_COMPONENT_COUNT, VECTOR_COMPONENT_COUNT);
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

        public double ScalarPart => this.scalarPart;
        public double[] VectorPart => this.vectorPart;
    }
}
