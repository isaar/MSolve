using System;
using ISAAR.MSolve.LinearAlgebra.Exceptions;

namespace ISAAR.MSolve.LinearAlgebra.Commons
{
    /// <summary>
    /// Implementations of analytic linear algebra formulas for small matrices and vectors.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal static class AnalyticFormulas
    {
        internal const double determinantTolerance = 1e-10;

        internal static double Matrix2x2ColMajorDeterminant(double[] matrix)
        {
            // a00 = matrix[0], a01 = matrix[2]
            // a10 = matrix[1], a11 = matrix[3]

            // Leibniz formula:
            return matrix[0] * matrix[3] - matrix[1] * matrix[2];
        }

        internal static (double[] inverse, double determinant) Matrix2x2ColMajorInvert(double[] matrix)
        {
            // a00 = matrix[0], a01 = matrix[2]
            // a10 = matrix[1], a11 = matrix[3]

            // Leibniz formula:
            double det = matrix[0] * matrix[3] - matrix[1] * matrix[2];
            if (Math.Abs(det) < determinantTolerance)
            {
                throw new SingularMatrixException($"Determinant == {det}");
            }

            // Cramer's rule: inverse = 1/det * [a11 -a01; -a10 a00]
            double[] inverse = new double[4];
            inverse[0] = matrix[3] / det;
            inverse[1] = -matrix[1] / det;
            inverse[2] = -matrix[2] / det;
            inverse[3] = matrix[0] / det;

            return (inverse, det);
        }

        internal static double Matrix3x3ColMajorDeterminant(double[] matrix)
        {
            // a00 = matrix[0], a01 = matrix[3], a02 = matrix[6]
            // a10 = matrix[1], a11 = matrix[4], a12 = matrix[7]
            // a20 = matrix[2], a21 = matrix[5], a22 = matrix[8]

            // Laplace formula: 
            // det = a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20)  + a02 * (a10 * a21 - a11 * a20)
            return matrix[0] * (matrix[4] * matrix[8] - matrix[7] * matrix[5])
                - matrix[3] * (matrix[1] * matrix[8] - matrix[7] * matrix[2])
                + matrix[6] * (matrix[1] * matrix[5] - matrix[4] * matrix[2]);
        }

        internal static (double[] inverse, double determinant) Matrix3x3ColMajorInvert(double[] matrix)
        {
            // a00 = matrix[0], a01 = matrix[3], a02 = matrix[6]
            // a10 = matrix[1], a11 = matrix[4], a12 = matrix[7]
            // a20 = matrix[2], a21 = matrix[5], a22 = matrix[8]

            // Minors of the first row (with signs)
            double c00 = matrix[4] * matrix[8] - matrix[7] * matrix[5];     // c00 = + a11*a22 - a12*a21
            double c01 = -matrix[1] * matrix[8] + matrix[7] * matrix[2];    // c01 = - a10*a22 + a12*a20
            double c02 = matrix[1] * matrix[5] - matrix[4] * matrix[2];     // c02 = + a10*a21 - a11*a20
            double det = matrix[0] * c00 + matrix[3] * c01 + matrix[6] * c02; // Laplace: det = a00*c00 + a01*c01 + a02*c02 
            if (Math.Abs(det) < determinantTolerance)
            {
                throw new SingularMatrixException($"Determinant == {det}");
            }

            // Cramer's rule: inverse = 1/det * transpose(C), C = matrix of minors
            double[] inverse = new double[9];
            inverse[0] = c00 / det; // inv[0,0]: c10
            inverse[1] = c01 / det; // inv[1,0]: c01
            inverse[2] = c02 / det; // inv[2,0]: c02
            inverse[3] = (-matrix[3] * matrix[8] + matrix[6] * matrix[5]) / det;    // inv[0,1]: c10 = - a01*a22 + a02*a21
            inverse[4] = (matrix[0] * matrix[8] - matrix[6] * matrix[2]) / det;     // inv[1,1]: c11 = + a00*a22 - a02*a20
            inverse[5] = (-matrix[0] * matrix[5] + matrix[3] * matrix[2]) / det;    // inv[2,1]: c12 = - a00*a21 + a01*a20
            inverse[6] = (matrix[3] * matrix[7] - matrix[6] * matrix[4]) / det;     // inv[0,2]: c20 = + a01*a12 - a02*a11
            inverse[7] = (-matrix[0] * matrix[7] + matrix[6] * matrix[7]) / det;    // inv[1,2]: c21 = - a00*a12 + a02*a12
            inverse[8] = (matrix[0] * matrix[4] - matrix[3] * matrix[1]) / det;     // inv[2,2]: c22 = + a00*a11 - a01*a10

            return (inverse, det);
        }
    }
}
