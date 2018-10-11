using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Points;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
    /// <summary>
    /// Enum class with the 3D integration rules for tetrahedra of varying orders. These are not tensor product of
    /// simple <see cref="GaussLegendre1D_old"/> rules. Quadrature rules were provided in https://www.code-aster.org/V2/doc/v11/fr/man_r/r3/r3.01.01.pdf
    /// Authors: Dimitris Tsapetis
    /// </summary>
    public sealed class PyramidQuadrature:IQuadrature3D
	{
		public static readonly PyramidQuadrature Points5 = new PyramidQuadrature(
			new GaussPoint3D(0.5, 0, 0.1531754163448146, 2.0 / 15.0),
			new GaussPoint3D(0, 0.5, 0.1531754163448146, 2.0 / 15.0),
			new GaussPoint3D(-0.5, 0, 0.1531754163448146, 2.0 / 15.0),
			new GaussPoint3D(0, -0.5, 0.1531754163448146, 2.0 / 15.0),
			new GaussPoint3D(0, 0, 0.6372983346207416, 2.0 / 15.0));

		public static readonly PyramidQuadrature Points6 = new PyramidQuadrature(
			new GaussPoint3D(0.5702963741068025, 0, 0.1666666666666666, 0.1024890634400000),
			new GaussPoint3D(0, 0.5702963741068025, 0.1666666666666666, 0.1024890634400000),
			new GaussPoint3D(-0.5702963741068025, 0, 0.1666666666666666, 0.1024890634400000),
			new GaussPoint3D(0, -0.5702963741068025, 0.1666666666666666, 0.1024890634400000),
			new GaussPoint3D(0, 0, 0.08063183038464675, 0.1100000000000000),
			new GaussPoint3D(0, 0, 0.6098484849057127, 0.1467104129066667));

		// Find and implement missing quadratures.

		private PyramidQuadrature(params GaussPoint3D[] points)
	    {
		    this.IntegrationPoints = new List<GaussPoint3D>(points);
	    }

	    public IReadOnlyList<GaussPoint3D> IntegrationPoints { get; }
	}
}
