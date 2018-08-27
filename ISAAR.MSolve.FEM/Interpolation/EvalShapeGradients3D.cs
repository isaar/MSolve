using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// Stores the 1st order shape function derivatives with respect to the global cartesian coordinates and the Jacobian
	/// of an interpolation, evaluated at a certain natural point of a finite element. These quantities are needed in many 
	/// places, thus passing an instance of this class is less verbose and error prone.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class EvalShapeGradients3D
    {
		private readonly double[][] shapeGradientsCartesian;

	    public EvalShapeGradients3D(double[,] shapeGradientNatural, Jacobian3D jacobian)
	    {
		    int numberOfNodes = shapeGradientNatural.GetLength(0);
			this.shapeGradientsCartesian=new double[numberOfNodes][];
		    for (int i = 0; i < numberOfNodes; i++)
		    {
			    this.shapeGradientsCartesian[i] = jacobian.TransformNaturalDerivativesToCartesian(shapeGradientNatural[i, 0],
				    shapeGradientNatural[i, 1], shapeGradientNatural[i, 2]);
		    }

		    this.Jacobian = jacobian;
	    }

		/// <summary>
		/// The inverse Jacobian matrix of the interpolation and its determinant.s
		/// </summary>
		public Jacobian3D Jacobian { get; }

		/// <summary>
		/// The values of the stored shape function derivatives, with respect to the global cartesian coordinates, that 
		/// correspond to the node with local index <paramref name="nodeIdx"/>.
		/// </summary>
		/// <param name="nodeIdx">The local index of the node, namely its order among the nodes of the finite element.</param>
		/// <returns></returns>
		public IReadOnlyList<double> this[int nodeIdx] => shapeGradientsCartesian[nodeIdx];
    }
}
