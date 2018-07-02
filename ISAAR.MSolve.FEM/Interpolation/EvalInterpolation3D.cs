using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// Stores the shape functions, 1st order derivatives with respect to the global cartesian coordinates and the Jacobian
	/// of an interpolation, evaluated at a certain natural point of a finite element. These quantities are needed in many 
	/// places, thus passing an instance of this class is less verbose and error prone.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class EvalInterpolation3D
    {
		private readonly double[] shapeFunctions;
	    private readonly double[][] shapeGradientCartesian;

	    public EvalInterpolation3D(double[] shapeFunctions, double[,] shapeGradientsNatural, Jacobian3D jacobian)
	    {
		    int numberOfNodes = shapeFunctions.Length;
		    if (shapeGradientsNatural.GetLength(0) != numberOfNodes) throw new ArgumentException($"There are {shapeFunctions.Length}"
		        + $" evaluated shape functions, but {shapeGradientsNatural.GetLength(0)} evaluated natural shape derivatives.");
		    this.shapeFunctions = shapeFunctions;
			this.shapeGradientCartesian= new double[numberOfNodes][];
		    for (int i = 0; i < numberOfNodes; i++)
		    {
			    this.shapeGradientCartesian[i] = jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[i, 0],
				    shapeGradientsNatural[i, 1], shapeGradientsNatural[i, 2]);
		    }
		    this.Jacobian = jacobian;
	    }

	    /// <summary>
	    /// The inverse Jacobian matrix of the interpolation and its determinant.
	    /// </summary>
	    public Jacobian3D Jacobian { get; }

	    public Matrix2D BuildShapeFunctionMatrix()
	    {
			var array2D= new double[3,3*shapeFunctions.Length];
		    for (int i = 0; i < shapeFunctions.Length; i++)
		    {
			    array2D[0, 3 * i] = shapeFunctions[i];
			    array2D[1, 2 * i + 1] = shapeFunctions[i];
			    array2D[2, 3 * i + 2] = shapeFunctions[i];
		    }
			return new Matrix2D(array2D);
	    }

	    /// <summary>
	    /// The value of the stored shape function that corresponds to the node with local index <paramref name="nodeIdx"/>.
	    /// </summary>
	    /// <param name="nodeIdx">The local index of the node, namely its order among the nodes of the finite element.</param>
	    /// <returns></returns>
		public double GetShapeFunction(int nodeIdx) => shapeFunctions[nodeIdx];

	    /// <summary>
	    /// The values of the stored shape function derivatives, with respect to the global cartesian coordinates, that 
	    /// correspond to the node with local index <paramref name="nodeIdx"/>.
	    /// </summary>
	    /// <param name="nodeIdx">The local index of the node, namely its order among the nodes of the finite element.</param>
	    /// <returns></returns>
		public IReadOnlyList<double> GetShapeGradientCartesian(int nodeIdx) => shapeGradientCartesian[nodeIdx];

	    public CartesianPoint3D TransformPointNaturalToGlobalCartesian(IReadOnlyList<Node3D> nodes, NaturalPoint3D naturalCoordinates)
	    {
		    if (nodes.Count != shapeFunctions.Length) throw new ArgumentException(
			    $"There are {shapeFunctions.Length} evaluated shape functions stored, but {nodes.Count} were passed in.");
		    double x = 0, y = 0, z = 0;
		    for (int i = 0; i < shapeFunctions.Length; i++)
		    {
			    Node3D node = nodes[i];
			    x += shapeFunctions[i] * node.X;
			    y += shapeFunctions[i] * node.Y;
			    z += shapeFunctions[i] * node.Z;
			}
			return new CartesianPoint3D(x,y,z);
	    }

	}
}
