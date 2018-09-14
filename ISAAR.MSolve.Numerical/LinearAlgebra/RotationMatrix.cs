using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Numerical.Matrices
{
    public class RotationMatrix
    {
        /**
         * Calculates a rotation matrix from three orthonormal vectors.
         *
         * @param vectorX
         *            The vector along the X axis
         * @param vectorY
         *            The vector along the Y axis
         * @param vectorZ
         *            The vector along the Z axis
         * @return The rotation matrix
         */
        public static Matrix2D CalculateFromOrthonormalVectors(Vector vectorX, Vector vectorY, Vector vectorZ)
        {
            var rotationMatrix = new Matrix2D(3, 3);
            rotationMatrix[0, 0] = vectorX[0];
            rotationMatrix[0, 1] = vectorY[0];
            rotationMatrix[0, 2] = vectorZ[0];

            rotationMatrix[1, 0] = vectorX[1];
            rotationMatrix[1, 1] = vectorY[1];
            rotationMatrix[1, 2] = vectorZ[1];

            rotationMatrix[2, 0] = vectorX[2];
            rotationMatrix[2, 1] = vectorY[2];
            rotationMatrix[2, 2] = vectorZ[2];

            return rotationMatrix;
        }

        /**
         * Calculates the rotations matrix given initial and rotated vectors.
         *
         * @param initialVector
         *            The initial vector
         * @param rotatedVector
         *            The rotated vector
         * @return The rotation matrix
         */
        public static Matrix2D CalculateRotationMatrix(Vector initialVector, Vector rotatedVector)
        {
            Vector initialVectorLocal = new Vector(3);
            Vector rotatedVectorLocal = new Vector(3);

            double normInitial = initialVector.Norm;
            double normRotated = rotatedVector.Norm;
            initialVector.CopyTo(initialVectorLocal.Data, 0);
            rotatedVector.CopyTo(rotatedVectorLocal.Data, 0);
            initialVectorLocal.Scale(1d / normInitial);
            rotatedVectorLocal.Scale(1d / normRotated);

            var vectorSum = new Vector(initialVectorLocal + rotatedVectorLocal);
            double vectorSumNorm = vectorSum.Norm;

            var rotationMatrix = new Matrix2D(3, 3);
            rotationMatrix[0, 0] = 1d;
            rotationMatrix[1, 1] = 1d;
            rotationMatrix[2, 2] = 1d;

            // @Theo
            rotationMatrix.LinearCombinationGOAT(new[] { 2d }, new[] { Matrix2D.FromVector(rotatedVectorLocal.Data) * Matrix2D.FromVectorTranspose(initialVectorLocal.Data) });
            rotationMatrix.LinearCombinationGOAT(new[] { -2d / (vectorSumNorm * vectorSumNorm) }, new[] { Matrix2D.FromVector(vectorSum.Data) * Matrix2D.FromVectorTranspose(vectorSum.Data) });
            return rotationMatrix;
        }

        public static Matrix2D CalculateRotationMatrixBeam2D(double phi)
        {
            var rotationMatrix = new Matrix2D(3, 3);
            rotationMatrix[0, 0] = Math.Cos(phi);
            rotationMatrix[0, 1] = -Math.Sin(phi);
            rotationMatrix[1, 0] = Math.Sin(phi);
            rotationMatrix[1, 1] = Math.Cos(phi);
            rotationMatrix[2, 2] = 1d;

            return rotationMatrix;
        }

        private RotationMatrix()
        {
            // DO Nothing
        }

    }
}
