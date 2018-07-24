using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a tetrahedral finite element with 10 nodes. Quadratic shape function.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis.
	/// </summary>
    public class InterpolationTet10:IsoparametricInterpolation3DBase
    {
	    private static readonly InterpolationTet10 uniqueInstance = new InterpolationTet10();

		private InterpolationTet10() : base(10)
	    {
		    NodalNaturalCoordinates = new NaturalPoint3D[]
		    {
			    new NaturalPoint3D(0,0,0),
			    new NaturalPoint3D(1,0,0),
			    new NaturalPoint3D(0,1,0),
			    new NaturalPoint3D(0,0,1),

			    new NaturalPoint3D(0.5,0,0),
			    new NaturalPoint3D(0.5,0.5,0),
			    new NaturalPoint3D(0,0.5,0),

			    new NaturalPoint3D(0,0,0.5),
			    new NaturalPoint3D(0,0.5,0.5),
			    new NaturalPoint3D(0.5,0,0.5),
			};
		}

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
	    public override IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <exception cref="InterpolationTet10"/> object for the whole program. Thread safe.
		/// </summary>
	    public static InterpolationTet10 UniqueInstance => uniqueInstance;

		/// <summary>
		/// The reverse mapping for this interpolation, namely from global cartesian coordinates to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> node)=> throw new NotImplementedException("Requires iterative procedure.");


		/// <summary>
		/// Returns the shape functions a tetrahedral quadratic element evaluated on a single point.
		/// Implementation is based on <see cref="https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch10.d/AFEM.Ch10.pdf">Carlos Felippa - Introduction to Finite Element Methods</see>
		/// </summary>
		/// <param name="xi"></param>
		/// <param name="eta"></param>
		/// <param name="zeta"></param>
		/// <returns></returns>
		protected override double[] EvaluateAt(double xi, double eta, double zeta)
		{
			var z1 = xi;
			var z2 = eta;
			var z3 = zeta;
			var z4 = 1 - z1 - z2 - z3;

			var values = new double[10];
			values[0] = z1 * (2 * z1 - 1);
			values[1] = z2 * (2 * z2 - 1);
			values[2] = z3 * (2 * z3 - 1);
			values[3] = z4 * (2 * z4 - 1);
			values[4] = 4 * z1 * z2;
			values[5] = 4 * z2 * z3;
			values[6] = 4 * z3 * z1;
			values[7] = 4 * z1 * z4;
			values[8] = 4 * z2 * z4;
			values[9] = 4 * z3 * z4;
			return values;
		}


	    protected override double[,] EvaluateGradientsAt(double xi, double eta, double zeta)
	    {
		    var derivatives = new double[10, 3];

		    derivatives[0, 0] = 4 * xi - 1;
		    derivatives[1, 0] = 0.0;
		    derivatives[2, 0] = 0.0;
		    derivatives[3, 0] = -3 + 4 * eta + 4 * zeta + 4 * xi;
		    derivatives[4, 0] = 4 * eta;
		    derivatives[5, 0] = 0.0;
		    derivatives[6, 0] = 4 * zeta;
		    derivatives[7, 0] = 4 - 8 * xi - 4 * eta - 4 * zeta;
		    derivatives[8, 0] = -4 * eta;
		    derivatives[9, 0] = -4 * zeta;

			derivatives[0, 1] = 0.0;
		    derivatives[1, 1] = 4*eta-1;
		    derivatives[2, 1] = 0.0;
		    derivatives[3, 1] = -3 + 4 * xi + 4 * zeta + 4 * eta;
		    derivatives[4, 1] = 4 * xi;
		    derivatives[5, 1] = 4 * zeta;
		    derivatives[6, 1] = 0.0;
		    derivatives[7, 1] = -4 * xi;
		    derivatives[8, 1] = 4 - 4 * xi - 8 * eta - 4 * zeta;
		    derivatives[9, 1] = -4 * zeta;

		    derivatives[0, 2] = 0.0;
		    derivatives[1, 2] = 0.0;
		    derivatives[2, 2] = 4 * zeta - 1;
		    derivatives[3, 2] = -3+4*eta+4*xi+4*zeta;
		    derivatives[4, 2] = 0.0;
		    derivatives[5, 2] = 4 * eta;
		    derivatives[6, 2] = 4 * xi;
		    derivatives[7, 2] = -4 * xi;
		    derivatives[8, 2] = -4 * eta;
		    derivatives[9, 2] = 4 - 4 * xi - 4 * eta - 8 * zeta;
		    return derivatives;
	    }
    }
}
