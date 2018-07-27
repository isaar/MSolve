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
			NodalNaturalCoordinates= new NaturalPoint3D[]
			{
				new NaturalPoint3D(0,0,0), 
				new NaturalPoint3D(1,0,0),
				new NaturalPoint3D(0,1,0),
				new NaturalPoint3D(0,0,1),
				new NaturalPoint3D(1,0,1),
				new NaturalPoint3D(0,1,1),
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

	    protected override double[] EvaluateAt(double xi, double eta, double zeta)
	    {
		    var values = new double[6];

		    values[0] = 1 / 6.0 * (1 + 2 * xi) * (1 - zeta);
		    values[1] = 1 / 6.0 * (1 - xi - Math.Sqrt(3) * eta) * (1 - zeta);
		    values[2] = 1 / 6.0 * (1 - xi + Math.Sqrt(3) * eta) * (1 - zeta);
			values[3] = 1 / 6.0 * (1 + 2 * xi) * (1 + zeta); 
		    values[4] = 1 / 6.0 * (1 - xi - Math.Sqrt(3) * eta) * (1 + zeta);
		    values[5] = 1 / 6.0 * (1 - xi + Math.Sqrt(3) * eta) * (1 + zeta); 

		    return values;
	    }

		// Evaluation based on: http://www.softeng.rl.ac.uk/st/projects/felib4/Docs/html/Level-0/wdg6/doc-wdg6.pdf
	    protected override double[,] EvaluateGradientsAt(double xi, double eta, double zeta)
	    {
		    var x = xi;
		    var y = eta;
		    var z = zeta;

		    var derivatives = new double[6, 3];

		    derivatives[0, 0] = 1 / 3.0 - z / 3;
		    derivatives[1, 0] = z / 6 - 1 / 6.0;
		    derivatives[2, 0] = z / 6 - 1 / 6.0;
		    derivatives[3, 0] = z / 3 + 1 / 3.0;
		    derivatives[4, 0] = -z / 6 - 1 / 6.0;
		    derivatives[5, 0] = -z / 6 - 1 / 6.0;

		    derivatives[0, 1] = 0.0;
		    derivatives[1, 1] = Math.Sqrt(3) * (z - 1) / 6;
		    derivatives[2, 1] = -Math.Sqrt(3) * (z - 1) / 6;
		    derivatives[3, 1] = 0.0;
		    derivatives[4, 1] = -Math.Sqrt(3) * (z + 1) / 6;
		    derivatives[5, 1] = Math.Sqrt(3) * (z + 1) / 6;

		    derivatives[0, 2] = -x / 3 - 1 / 6.0;
		    derivatives[1, 2] = x / 6 + Math.Sqrt(3) * y / 6 - 1 / 6.0;
		    derivatives[2, 2] = x / 6 - Math.Sqrt(3) * y / 6 - 1 / 6.0;
		    derivatives[3, 2] = x / 3 + 1 / 6.0;
		    derivatives[4, 2] = 1 / 6.0 - Math.Sqrt(3) * y / 6 - x / 6;
		    derivatives[5, 2] = Math.Sqrt(3) * y / 6 - x / 6 + 1 / 6.0;


		    return derivatives;
	    }
    }
}
