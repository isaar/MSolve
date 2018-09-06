using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;

//TODO: need special 2x2 Matrix class
//TODO: once we know that an exception will be thrown, try to pinpoint the error: wrong node order, clockwise node order, the  
//      element's shape is too distorted, midpoints are too close to corners in quadratic elements, etc...
namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// This class encapsulates the determinant and inverse of the Jacobian matrix for a 2D mapping.
    /// Let f be a mapping: x \in R^2 -> f(x) \in R^2. The Jacobian matrix of the mapping is 
    /// J = [df_1/dx_1 df_1/dx_2; df_2/dx_1 df_2/dx_2]. Note that some sources call the transpose of this matrix as J. In FEM we 
    /// are usually interested in the determinant and inverse of the Jacobian matrix. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Jacobian2D 
    {
        private const double determinantTolerance = 1E-8; // This needs to be in a static settings class.

        private readonly double[,] inverseJ;

        /// <summary>
        /// The caller (usually the interpolation class) assumes responsibility for matching the nodes to the shape function 
        /// derivatives.
        /// </summary>
        /// <param name="nodes">The nodes used for the interpolation.</param>
        /// <param name="naturalDerivatives">The shape function derivatives at a specific integration point.</param>
        public Jacobian2D(IReadOnlyList<Node2D> nodes, double[,] naturalDerivatives)
        {
            // The original matrix is not stored. Only the inverse and the determinant
            double[,] jacobianMatrix = CalculateJacobianMatrix(nodes, naturalDerivatives);
            (inverseJ, Determinant) = InvertAndDeterminant(jacobianMatrix);
            if (Determinant < determinantTolerance)
            {
                throw new ArgumentException("Jacobian determinant is negative or under the allowed tolerance"
                    + $" ({Determinant} < {determinantTolerance}). Check the order of nodes or the element geometry.");
            }
        }

        /// <summary>
        /// The determinant of the original Jacobian matrix, not its inverse.
        /// </summary>
        public double Determinant { get; }

        public double[] TransformNaturalDerivativesToCartesian(double[] naturalGradient)
        {
            // naturalGradient * inverseJ = 1-by-2 * 2-by-2 = 1-by-2
            var result = new double[2];
            result[0] = naturalGradient[0] * inverseJ[0, 0] + naturalGradient[1] * inverseJ[1, 0];
            result[1] = naturalGradient[0] * inverseJ[0, 1] + naturalGradient[1] * inverseJ[1, 1];
            return result;
        }

        public double[] TransformNaturalDerivativesToCartesian(double derivativeXi, double derivativeEta)
        {
            // naturalGradient * inverseJ = 1-by-2 * 2-by-2 = 1-by-2
            var result = new double[2];
            result[0] = derivativeXi * inverseJ[0, 0] + derivativeEta * inverseJ[1, 0];
            result[1] = derivativeXi * inverseJ[0, 1] + derivativeEta * inverseJ[1, 1];
            return result;
        }

        private static double[,] CalculateJacobianMatrix(IReadOnlyList<Node2D> nodes, double[,] naturalDerivatives)
        {
            var J = new double[2, 2];
            for (int nodeIndex = 0; nodeIndex < nodes.Count; ++nodeIndex)
            {
                double x = nodes[nodeIndex].X;
                double y = nodes[nodeIndex].Y;
                double N_xi = naturalDerivatives[nodeIndex, 0];
                double N_eta = naturalDerivatives[nodeIndex, 1];

                J[0, 0] += N_xi * x;
                J[0, 1] += N_eta * x;
                J[1, 0] += N_xi * y;
                J[1, 1] += N_eta * y;
            }
            return J;
        }

        private static (double[,] inverse, double det) InvertAndDeterminant(double[,] matrix)
        {
            // Leibniz formula:
            double det = matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
            if (Math.Abs(det) < determinantTolerance) throw new Exception(
                $"|Determinant| = {Math.Abs(det)} < tolerance = {determinantTolerance}. The matrix is singular");

            // Cramer's rule: inverse = 1/det * [a11 -a01; -a10 a00]
            double[,] inverse = new double[,] 
            { 
                { matrix[1, 1] / det, -matrix[0, 1] / det }, 
                { -matrix[1, 0] / det, matrix[0, 0] / det }
            };

            return (inverse, det);
        }
    }
}
