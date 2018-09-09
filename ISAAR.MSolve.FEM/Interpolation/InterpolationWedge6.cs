using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a wedge finite element with 6 nodes. Linear shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
    public class InterpolationWedge6:IsoparametricInterpolation3DBase
    {
		private static readonly InterpolationWedge6 uniqueInstance= new InterpolationWedge6();

	    private InterpolationWedge6() : base(6)
	    {
		    NodalNaturalCoordinates = new NaturalPoint3D[]
		    {
			    new NaturalPoint3D(-1, 1, 0),
			    new NaturalPoint3D(-1, 0, 1),
			    new NaturalPoint3D(-1, 0, 0),
			    new NaturalPoint3D(1, 1, 0),
			    new NaturalPoint3D(1, 0, 1),
			    new NaturalPoint3D(1, 0, 0),
		    };
	    }

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationWedge6"/> object for the whole program. Thread safe.
		/// </summary>
		public static InterpolationWedge6 UniqueInstance => uniqueInstance;

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> node)=> throw new NotImplementedException("Implementation pending");

	    protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
	    {
		    var values = new double[6];

		    values[0] = 0.5 * eta * (1 - xi);
		    values[1] = 0.5 * zeta * (1 - xi);
		    values[2] = 0.5 * (1 - eta - zeta) * (1 - xi);
		    values[3] = 0.5 * eta * (xi + 1);
		    values[4] = 0.5 * zeta * (xi + 1);
		    values[5] = 0.5*(1-eta-zeta)*(xi+1); 

		    return values;
	    }

		// Evaluation based on: https://www.code-aster.org/V2/doc/v11/en/man_r/r3/r3.01.01.pdf
	    protected sealed override double[,] EvaluateGradientsAt(double xi, double eta, double zeta)
	    {
		    var x = xi;
		    var y = eta;
		    var z = zeta;

		    var derivatives = new double[6, 3];

		    derivatives[0, 0] = -y / 2;
		    derivatives[1, 0] = -z / 2;
		    derivatives[2, 0] = (y  + z  - 1) / 2;
		    derivatives[3, 0] = y / 2;
		    derivatives[4, 0] = z / 2;
		    derivatives[5, 0] =  (1 - z  - y) / 2;

		    derivatives[0, 1] = (1- x) / 2;
		    derivatives[1, 1] = 0.0;
		    derivatives[2, 1] = (x - 1) / 2;
		    derivatives[3, 1] = (x + 1) / 2;
		    derivatives[4, 1] = 0.0;
		    derivatives[5, 1] = (-x - 1) / 2;

		    derivatives[0, 2] = 0.0;
		    derivatives[1, 2] = (1 - x) / 2;
		    derivatives[2, 2] = (x - 1) / 2;
		    derivatives[3, 2] = 0.0;
		    derivatives[4, 2] = (x + 1) / 2;
		    derivatives[5, 2] = (-x - 1) / 2;

		    return derivatives;
	    }
    }
}
