using System;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Geometry
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
        public static Matrix CalculateFromOrthonormalVectors(double[] vectorX, double[] vectorY, double[] vectorZ)
        {
            var rotationMatrix = Matrix.CreateZero(3, 3);
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
        public static Matrix CalculateRotationMatrix(double[] initialVector, double[] rotatedVector)
        {
            double[] initialVectorLocal = initialVector.Copy();
            double[] rotatedVectorLocal = rotatedVector.Copy();

            double normInitial = initialVector.Norm2();
            double normRotated = rotatedVector.Norm2();
            initialVectorLocal.ScaleIntoThis(1d / normInitial);
            rotatedVectorLocal.ScaleIntoThis(1d / normRotated);

            double[] vectorSum = initialVectorLocal.Add(rotatedVectorLocal);
            double vectorSumNorm = vectorSum.Norm2();

            var rotationMatrix = Matrix.CreateZero(3, 3);
            rotationMatrix[0, 0] = 1d;
            rotationMatrix[1, 1] = 1d;
            rotationMatrix[2, 2] = 1d;

            // @Theo
            rotationMatrix.AxpyIntoThis(rotatedVectorLocal.TensorProduct(initialVectorLocal), 2.0);
            rotationMatrix.AxpyIntoThis(vectorSum.TensorProduct(vectorSum), -2.0 / (vectorSumNorm * vectorSumNorm));
            //rotationMatrix.LinearCombinationGOAT(new[] { 2d }, new[] { Matrix2D.FromVector(rotatedVectorLocal.Data) * Matrix2D.FromVectorTranspose(initialVectorLocal.Data) });
            //rotationMatrix.LinearCombinationGOAT(new[] { -2d / (vectorSumNorm * vectorSumNorm) }, new[] { Matrix2D.FromVector(vectorSum.Data) * Matrix2D.FromVectorTranspose(vectorSum.Data) });

            return rotationMatrix;
        }

        public static Matrix CalculateRotationMatrixBeam2D(double phi)
        {
            var rotationMatrix = Matrix.CreateZero(3, 3);
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
