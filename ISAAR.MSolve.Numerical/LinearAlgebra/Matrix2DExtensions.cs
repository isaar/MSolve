using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public static class Matrix2DExtensions
    {
        public static (Matrix2D inverse, double determinant) Invert2x2AndDeterminant(this Matrix2D matrix, 
            double determinantTolerance = 1E-10)
        {
            if ((matrix.Rows != 2) || (matrix.Rows != 2)) throw new ArgumentException("This only works for 2-by-2 matrices");

            // Leibniz formula:
            double det = matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
            if (Math.Abs(det) < determinantTolerance)
            {
                throw new ArgumentException($"This matrix is considered singular, because its determinant = {det}" 
                    + " is under the tolerance {determinantTolerance}.");
            }

            // Cramer's rule: inverse = 1/det * [a11 -a01; -a10 a00]
            double[,] inverse = new double[,] 
            { 
                { matrix[1, 1] / det, -matrix[0, 1] / det }, 
                { -matrix[1, 0] / det, matrix[0, 0] / det }
            };
            return (new Matrix2D(inverse), det);
        }
    }
}
