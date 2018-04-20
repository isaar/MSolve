using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Interpolation
{
    class Jacobian2D //TODO: need special 2x2 Matrix class probably
    {
        private const double DETERMINANT_TOLERANCE = 0.00000001; // This needs to be in a static settings class.
        private const int DIMENSION = 2;

        private readonly Matrix2by2 inverseJ; // I need a Matrix view for this
        public double Determinant { get; }

        /// <summary>
        /// Caller (the interpolation class) assumes responsibility for matching the nodes to the shape function 
        /// derivatives.
        /// </summary>
        /// <param name="nodes">The nodes used for the interpolation.</param>
        /// <param name="naturalDerivatives">The shape function derivatives at a specific integration point. </param>
        public Jacobian2D(IReadOnlyList<ICartesianPoint2D> nodes, double[,] naturalDerivatives)
        {
            // The original matrix is not stored. Only the inverse and the determinant
            Matrix2by2 jacobianMatrix = CalculateJacobianMatrix(nodes, naturalDerivatives);
            (inverseJ, Determinant) = jacobianMatrix.InvertAndDetermninant();
            if (Determinant < DETERMINANT_TOLERANCE)
            {
                throw new ArgumentException("Jacobian determinant is negative or under the allowed tolerance"
                    + $" ({Determinant} < {DETERMINANT_TOLERANCE}). Check the order of nodes or the element geometry.");
            }
        }

        public Vector2 TransformNaturalDerivativesToCartesian(Vector2 naturalGradient)
        {
            return naturalGradient * inverseJ;
        }

        private static Matrix2by2 CalculateJacobianMatrix(IReadOnlyList<ICartesianPoint2D> nodes, 
            double[,] naturalDerivatives)
        {
            var J = new double[DIMENSION, DIMENSION];
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
            return Matrix2by2.CreateFromArray(J);
        }
    }
}
