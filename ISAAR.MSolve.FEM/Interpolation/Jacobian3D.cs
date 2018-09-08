using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// This class encapsulates the determinant and inverse of the Jacobian matrix for a 3D mapping.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class Jacobian3D
	{
		private const double determinantTolerance = 1E-8;

		private readonly double[,] inverseJ;

		/// <summary>
		/// The caller (usually the interpolation class) assumes responsibility for matching the nodes to the shape function 
		/// derivatives.
		/// </summary>
		/// <param name="nodes">The nodes used for the interpolation.</param>
		/// <param name="naturalCoordinates">The shape function derivatives at a specific integration point.</param>
		public Jacobian3D(IReadOnlyList<Node3D> nodes, double[,] naturalDerivatives)
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

		private static double[,] CalculateJacobianMatrix(IReadOnlyList<Node3D> nodes, double[,] naturalDerivatives)
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

		private static (double[,] inverse, double determinant) InvertAndDeterminant(double[,] jacobianMatrix)
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

			return (inverseJacobian,determinant);
		}
	}
}
