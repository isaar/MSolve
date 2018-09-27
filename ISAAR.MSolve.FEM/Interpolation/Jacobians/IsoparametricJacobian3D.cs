using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.FEM.Interpolation.Jacobians
{
    /// <summary>
    /// This class encapsulates the determinant and inverse of the Jacobian matrix for a 3D isoparametric mapping.
    /// Let f be a mapping: x \in R^3 -> f(x) \in R^3. The Jacobian matrix of the mapping is: 
    /// J = [df_1/dx_1 df_1/dx_2 df_1/dx_3; df_2/dx_1 df_2/dx_2 df_2/dx_3; df_3/dx_1 df_3/dx_2 df_3/dx_3]. 
    /// Note that some sources call the transpose of this matrix as J. In FEM we are usually interested in the determinant and 
    /// inverse of the Jacobian matrix.
    /// Authors: Dimitris Tsapetis
    /// </summary>
    public class IsoparametricJacobian3D
	{
		private const double determinantTolerance = 1E-8;

		private readonly Matrix2D inverseJ;

		/// <summary>
		/// The caller (usually the interpolation class) assumes responsibility for matching the nodes to the shape function 
		/// derivatives.
		/// </summary>
		/// <param name="nodes">The nodes used for the interpolation.</param>
		/// <param name="naturalCoordinates">The shape function derivatives at a specific integration point.</param>
		public IsoparametricJacobian3D(IReadOnlyList<Node3D> nodes, Matrix2D naturalDerivatives)
		{
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

        /// <summary>
        /// Transforms the gradient of a vector-valued function from the natural to the global cartesian coordinate system.
        /// </summary>
        /// <param name="naturalGradient">The gradient of a vector-valued function in the natural coordinate system. Each row 
        ///     corresponds to the gradient of a single component of the vector function. Each column corresponds to the 
        ///     derivatives of all components with respect to a single coordinate.</param>
        public Matrix2D TransformNaturalDerivativesToCartesian(Matrix2D naturalGradient) => naturalGradient * inverseJ;

        public double[] TransformNaturalDerivativesToCartesian(double[] naturalGradient)
		{
			var result = new double[3];
			result[0] = naturalGradient[0] * inverseJ[0, 0] + naturalGradient[1] * inverseJ[1, 0] +
			            naturalGradient[2] * inverseJ[2, 0];

			result[1] = naturalGradient[0] * inverseJ[0, 1] + naturalGradient[1] * inverseJ[1, 1] +
			            naturalGradient[2] * inverseJ[2, 1];

			result[2] = naturalGradient[0] * inverseJ[0, 2] + naturalGradient[1] * inverseJ[1, 2] +
			            naturalGradient[2] * inverseJ[2, 2];
			return result;
		}

		public double[] TransformNaturalDerivativesToCartesian(double derivatiXi, double derivativeEta, double derivativeZeta)
		{
			var result = new double[3];
			result[0] = derivatiXi * inverseJ[0, 0] + derivativeEta * inverseJ[1, 0] +
			            derivativeZeta * inverseJ[2, 0];

			result[1] = derivatiXi * inverseJ[0, 1] + derivativeEta * inverseJ[1, 1] +
			            derivativeZeta * inverseJ[2, 1];

			result[2] = derivatiXi * inverseJ[0, 2] + derivativeEta * inverseJ[1, 2] +
			            derivativeZeta * inverseJ[2, 2];
			return result;
		}

		private static double[,] CalculateJacobianMatrix(IReadOnlyList<Node3D> nodes, Matrix2D naturalDerivatives)
		{
			var jacobianMatrix = new double[3,3];

			for (int nodeIndex = 0; nodeIndex < nodes.Count; nodeIndex++)
			{
				jacobianMatrix[0, 0] += naturalDerivatives[nodeIndex, 0] * nodes[nodeIndex].X;
				jacobianMatrix[0, 1] += naturalDerivatives[nodeIndex, 1] * nodes[nodeIndex].X;
				jacobianMatrix[0, 2] += naturalDerivatives[nodeIndex, 2] * nodes[nodeIndex].X;

				jacobianMatrix[1, 0] += naturalDerivatives[nodeIndex, 0] * nodes[nodeIndex].Y;
				jacobianMatrix[1, 1] += naturalDerivatives[nodeIndex, 1] * nodes[nodeIndex].Y;
				jacobianMatrix[1, 2] += naturalDerivatives[nodeIndex, 2] * nodes[nodeIndex].Y;

				jacobianMatrix[2, 0] += naturalDerivatives[nodeIndex, 0] * nodes[nodeIndex].Z;
				jacobianMatrix[2, 1] += naturalDerivatives[nodeIndex, 1] * nodes[nodeIndex].Z;
				jacobianMatrix[2, 2] += naturalDerivatives[nodeIndex, 2] * nodes[nodeIndex].Z;
			}

			return jacobianMatrix;
		}

		private static (Matrix2D inverse, double determinant) InvertAndDeterminant(double[,] jacobianMatrix)
		{
			double determinant = jacobianMatrix[0, 0] *
			                     (jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[2, 1] * jacobianMatrix[1, 2])
			                     - jacobianMatrix[0, 1] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 2] -
			                                               jacobianMatrix[2, 0] * jacobianMatrix[1, 2])
			                     + jacobianMatrix[0, 2] * (jacobianMatrix[1, 0] * jacobianMatrix[2, 1] -
			                                               jacobianMatrix[2, 0] * jacobianMatrix[1, 1]);
			if (Math.Abs(determinant) < determinantTolerance) throw new Exception(
				$"|Determinant| = {Math.Abs(determinant)} < tolerance = {determinantTolerance}. The matrix is singular");

			var inverseJacobian = new double[3,3];

			inverseJacobian[0, 0] = (jacobianMatrix[1, 1] * jacobianMatrix[2, 2] - jacobianMatrix[1, 2] * jacobianMatrix[2, 1]) / determinant;
			inverseJacobian[0, 1] = (jacobianMatrix[0, 2] * jacobianMatrix[2, 1] - jacobianMatrix[0, 1] * jacobianMatrix[2, 2]) / determinant;
			inverseJacobian[0, 2] = (jacobianMatrix[0, 1] * jacobianMatrix[1, 2] - jacobianMatrix[0, 2] * jacobianMatrix[1, 1]) / determinant;

			inverseJacobian[1, 0] = (jacobianMatrix[1, 2] * jacobianMatrix[2, 0] - jacobianMatrix[1, 0] * jacobianMatrix[2, 2]) / determinant;
			inverseJacobian[1, 1] = (jacobianMatrix[0, 0] * jacobianMatrix[2, 2] - jacobianMatrix[0, 2] * jacobianMatrix[2, 0]) / determinant;
			inverseJacobian[1, 2] = (jacobianMatrix[0, 2] * jacobianMatrix[1, 0] - jacobianMatrix[0, 0] * jacobianMatrix[1, 2]) / determinant;

			inverseJacobian[2, 0] = (jacobianMatrix[1, 0] * jacobianMatrix[2, 1] - jacobianMatrix[1, 1] * jacobianMatrix[2, 0]) / determinant;
			inverseJacobian[2, 1] = (jacobianMatrix[0, 1] * jacobianMatrix[2, 0] - jacobianMatrix[0, 0] * jacobianMatrix[2, 1]) / determinant;
			inverseJacobian[2, 2] = (jacobianMatrix[0, 0] * jacobianMatrix[1, 1] - jacobianMatrix[0, 1] * jacobianMatrix[1, 0]) / determinant;

			return (new Matrix2D(inverseJacobian), determinant);
		}
	}
}
