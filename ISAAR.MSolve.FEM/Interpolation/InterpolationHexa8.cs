using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a hexahedral finite element with 8 nodes. Linear shape functions.
	/// Implements singleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
    public class InterpolationHexa8:IsoparametricInterpolation3DBase
    {
		private static readonly InterpolationHexa8 uniqueInstance= new InterpolationHexa8();

	    private InterpolationHexa8() : base(8)
	    {
		    NodalNaturalCoordinates = new NaturalPoint3D[]
		    {
			    new NaturalPoint3D(-1, -1, -1),
			    new NaturalPoint3D(1, -1, -1),
			    new NaturalPoint3D(1, 1, -1),
			    new NaturalPoint3D(-1, 1, -1),
			    new NaturalPoint3D(-1, -1, 1),
			    new NaturalPoint3D(1, -1, 1),
			    new NaturalPoint3D(1, 1, 1),
			    new NaturalPoint3D(-1, 1, 1),
		    };
	    }

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
		/// nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
		public override IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique <see cref="InterpolationHexa8"/> object for the whole program. Thread safe.
		/// </summary>
	    public static InterpolationHexa8 UniqueInstance => uniqueInstance;


	    /// <summary>
	    /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
	    /// </summary>
	    /// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
	    /// <returns></returns>
	    public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> nodes) =>
		    new InverseInterpolationHexa8(nodes);

	    protected sealed  override double[] EvaluateAt(double xi, double eta, double zeta)
	    {
		    var values = new double[8];
		    values[0] = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
		    values[1] = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
		    values[2] = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
		    values[3] = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
		    values[4] = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
		    values[5] = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
		    values[6] = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
		    values[7] = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);
		    return values;
	    }

	    protected sealed override double[,] EvaluateGradientsAt(double xi, double eta, double zeta)
	    {
		    var derivatives = new double[8, 3];
		    derivatives[0, 0] = -1 / 8.0 * (1 - eta) * (1 - zeta);
		    derivatives[1, 0] = 1 / 8.0 * (1 - eta) * (1 - zeta);
		    derivatives[2, 0] = 1 / 8.0 * (1 + eta) * (1 - zeta);
		    derivatives[3, 0] = -1 / 8.0 * (1 + eta) * (1 - zeta);
		    derivatives[4, 0] = -1 / 8.0 * (1 - eta) * (1 + zeta);
		    derivatives[5, 0] = 1 / 8.0 * (1 - eta) * (1 + zeta);
		    derivatives[6, 0] = 1 / 8.0 * (1 + eta) * (1 + zeta);
		    derivatives[7, 0] = -1 / 8.0 * (1 + eta) * (1 + zeta);

		    derivatives[0, 1] = -1 / 8.0 * (1 - xi) * (1 - zeta);
		    derivatives[1, 1] = -1 / 8.0 * (1 + xi) * (1 - zeta);
		    derivatives[2, 1] = 1 / 8.0 * (1 + xi) * (1 - zeta);
		    derivatives[3, 1] = 1 / 8.0 * (1 - xi) * (1 - zeta);
		    derivatives[4, 1] = -1 / 8.0 * (1 - xi) * (1 + zeta);
		    derivatives[5, 1] = -1 / 8.0 * (1 + xi) * (1 + zeta);
		    derivatives[6, 1] = 1 / 8.0 * (1 + xi) * (1 + zeta);
		    derivatives[6, 1] = 1 / 8.0 * (1 - xi) * (1 + zeta);

		    derivatives[0, 2] = -1 / 8.0 * (1 - xi) * (1 - eta);
		    derivatives[1, 2] = -1 / 8.0 * (1 + xi) * (1 - eta);
		    derivatives[2, 2] = -1 / 8.0 * (1 + xi) * (1 + eta);
		    derivatives[3, 2] = -1 / 8.0 * (1 - xi) * (1 + eta);
		    derivatives[4, 2] = 1 / 8.0 * (1 - xi) * (1 - eta);
		    derivatives[5, 2] = 1 / 8.0 * (1 + xi) * (1 - eta);
		    derivatives[6, 2] = 1 / 8.0 * (1 + xi) * (1 + eta);
		    derivatives[7, 2] = 1 / 8.0 * (1 - xi) * (1 + eta);
		    return derivatives;
	    }
    }
}
