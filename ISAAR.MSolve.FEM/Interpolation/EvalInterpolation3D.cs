using System;
using System.Collections.Generic;
using System.Text;
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
			var array2D= new double();
	    }
	}
}
