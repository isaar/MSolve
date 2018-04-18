using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Interpolation
{
    class Jacobian2D //TODO: need special 2x2 Matrix class probably
    {
        private const double DETERMINANT_TOLERANCE = 0.00000001; // This needs to be in a static settings class.
        private const int DIMENSION = 2;

        private readonly Matrix inverseJ; // I need a Matrix view for this
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
            Matrix jacobianMatrix = CalculateJacobianMatrix(nodes, naturalDerivatives);
            (inverseJ, Determinant) = jacobianMatrix.InvertAndDetermninant();
            if (Determinant < DETERMINANT_TOLERANCE)
            {
                throw new ArgumentException(String.Format(
                    "Jacobian determinant is negative or under tolerance ({0} < {1}). Check the order of nodes or the element geometry.",
                    Determinant, DETERMINANT_TOLERANCE));
            }
        }

        public Tuple<double, double> TransformNaturalDerivativesToCartesian(double derivativeXi, double derivativeEta)
        {
            double derivativeX = derivativeXi * inverseJ[0, 0] + derivativeEta * inverseJ[0, 1];
            double derivativeY = derivativeXi * inverseJ[1, 0] + derivativeEta * inverseJ[1, 1];
            return new Tuple<double, double>(derivativeX, derivativeY);
        }

        private static Matrix CalculateJacobianMatrix(IReadOnlyList<ICartesianPoint2D> nodes, 
            double[,] naturalDerivatives)
        {
            var J = Matrix.CreateZero(DIMENSION, DIMENSION);
            for (int nodeIndex = 0; nodeIndex < nodes.Count; ++nodeIndex)
            {
                double x = nodes[nodeIndex].X;
                double y = nodes[nodeIndex].Y;
                double N_xi = naturalDerivatives[nodeIndex, 0];
                double N_eta = naturalDerivatives[nodeIndex, 1];

                J[0, 0] += N_xi * x;
                J[0, 1] += N_xi * y;
                J[1, 0] += N_eta * x;
                J[1, 1] += N_eta * y;
            }
            return J;
        }
    }
}
