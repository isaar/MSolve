using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// Stores the shape functions of an interpolation, evaluated at a certain natural point of a finite element.
	/// Authors: Dimitris Tsapetis
	/// </summary>
    public class EvalShapeFunctions3D
	{
		private readonly double[] shapeFunctions;

		public EvalShapeFunctions3D(double[] shapeFunctions)
		{
			this.shapeFunctions = shapeFunctions;
		}

		/// <summary>
		/// The value of the stored shape function that corresponds to the node with local index <paramref name="nodeIdx"/>.
		/// </summary>
		/// <param name="nodeIdx">The local index of the node, namely its order among the nodes of the finite element.</param>
		/// <returns></returns>
		public double this[int nodeIdx] => shapeFunctions[nodeIdx];

		/// <summary>
		/// The shape function matrix is a 3-by3n, where n = is the number of shape functions. Row 0 corresponds to dof X,
		/// row1 corresponds to dof Y and row 2 to dof Z.
		/// </summary>
		/// <returns></returns>
		public Matrix2D BuildShapeFunctionMatrix()
		{
			var array2D = new double[3, 3 * shapeFunctions.Length];
			for (int i = 0; i < shapeFunctions.Length; i++)
			{
				array2D[0, 3 * i] = shapeFunctions[i];
				array2D[1, 3 * i + 1] = shapeFunctions[i];
				array2D[2, 3 * i + 2] = shapeFunctions[i];
			}

			return new Matrix2D(array2D);
		}


		public CartesianPoint3D TransformPointToGlobalCartesian(IReadOnlyList<Node3D> nodes,
			NaturalPoint3D naturalCoordinates)
		{
			if (nodes.Count != shapeFunctions.Length) throw new ArgumentException(
				$"There are {shapeFunctions.Length} evaluated shape functions stored, but {nodes.Count} were passed in.");
			double x = 0, y = 0, z = 0;
			for (int i = 0; i < shapeFunctions.Length; ++i)
			{
				Node3D node = nodes[i];
				double val = shapeFunctions[i];
				x += val * node.X;
				y += val * node.Y;
				z += val * node.Z;
			}

			return new CartesianPoint3D(x, y, z);
		}
	}
}
