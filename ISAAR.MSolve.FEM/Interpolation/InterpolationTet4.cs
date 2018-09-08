using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.FEM.Interpolation
{
	/// <summary>
	/// Isoparametric interpolation of a tetrahedral finite element with 4 nodes. Linear shape functions.
	/// Implements sigleton pattern.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InterpolationTet4:IsoparametricInterpolation3DBase
    {
		private static readonly  InterpolationTet4 uniqueInstance= new InterpolationTet4();

	    private InterpolationTet4() : base(4)
	    {
			NodalNaturalCoordinates= new NaturalPoint3D[]
			{
				new NaturalPoint3D(0,0,0),
				new NaturalPoint3D(1,0,0),
				new NaturalPoint3D(0,1,0),
				new NaturalPoint3D(0,0,1),
			};
	    }

		/// <summary>
		/// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order
		/// of these nodes matches the order of the shape functions and is always the same for each element.
		/// </summary>
	    public override IReadOnlyList<NaturalPoint3D> NodalNaturalCoordinates { get; }

		/// <summary>
		/// Get the unique instance <see cref="InterpolationTet4"/> object for the whole program. Thread safe.
		/// </summary>
	    public static InterpolationTet4 UniqueInstance => uniqueInstance;

		/// <summary>
		/// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
		/// </summary>
		/// <param name="node">The nodes of the finite element in the global cartesian coordinate system.</param>
		/// <returns></returns>
		// TODO: Find and implement inverse mapping for Tet4.
	    public override IInverseInterpolation3D CreateInverseMappingFor(IReadOnlyList<Node3D> node) => throw new NotImplementedException("Not implemented yet.");

		/// <summary>
		/// Returns the shape functions a tetrahedral linear element evaluated on a single point.
		/// Implementation is based on <see cref="https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf">Carlos Felippa - Introduction to Finite Element Methods</see>
		/// </summary>
		/// <param name="xi"></param>
		/// <param name="eta"></param>
		/// <param name="zeta"></param>
		/// <returns></returns>
		protected sealed override double[] EvaluateAt(double xi, double eta, double zeta)
	    {
		    var values = new double[4];
		    values[0] = 1 - xi - eta - zeta;
			values[1] = xi;
		    values[2] = eta;
		    values[3] = zeta;
		    
		    return values;
	    }

	    protected sealed override double[,] EvaluateGradientsAt(double xi, double eta, double zeta)
	    {
		    var derivatives = new double[4, 3];
		    derivatives[0, 0] = -1.0;
		    derivatives[0, 1] = -1.0;
		    derivatives[0, 2] = -1.0;

			derivatives[1, 0] = 1.0;
		    derivatives[1, 1] = 0.0;
		    derivatives[1, 2] = 0.0;

			derivatives[2, 0] = 0.0;
		    derivatives[2, 1] = 1.0;
		    derivatives[2, 2] = 0.0;

			derivatives[3, 0] = 0.0;
		    derivatives[3, 1] = 0.0;
		    derivatives[3, 2] = 1.0;
		    
		    return derivatives;
	    }
    }
}
