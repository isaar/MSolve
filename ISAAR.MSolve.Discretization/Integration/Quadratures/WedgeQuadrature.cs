using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
	/// <summary>
	/// Enum class with the 2D integration rules for wedges of varying orders. These are not tensor product of
	/// simple <see cref="GaussLegendre1D_old"/> rules. Quadrature rules were provided in https://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public sealed class WedgeQuadrature:IQuadrature3D
	{
		/// <summary>
		/// 6-integration point quadrature for order 3 of axis x, and order 2 for both axis y and z.
		/// </summary>
		public static readonly WedgeQuadrature Points6 = new WedgeQuadrature(
			new GaussPoint3D(-1 / Math.Sqrt(3), 0.5, 0.5, 1 / 6.0),
			new GaussPoint3D(-1 / Math.Sqrt(3), 0, 0.5, 1 / 6.0),
			new GaussPoint3D(-1 / Math.Sqrt(3), 0.5, 0, 1 / 6.0),
			new GaussPoint3D(1 / Math.Sqrt(3), 0.5, 0.5, 1 / 6.0),
			new GaussPoint3D(1 / Math.Sqrt(3), 0, 0.5, 1 / 6.0),
			new GaussPoint3D(1 / Math.Sqrt(3), 0.5, 0, 1 / 6.0));

		/// <summary>
		/// 8-integration point rule for order 3 for all axes.
		/// </summary>
		public static readonly WedgeQuadrature Points8 = new WedgeQuadrature(
			new GaussPoint3D(-0.577350269189626, 1 / 3.0, 1 / 3.0, -27.0 / 96.0),
			new GaussPoint3D(-0.577350269189626, 0.6, 0.2, 25.0 / 96.0),
			new GaussPoint3D(-0.577350269189626, 0.2, 0.6, 25.0 / 96.0),
			new GaussPoint3D(-0.577350269189626, 0.2, 0.0, 25.0 / 96.0),
			new GaussPoint3D(0.577350269189626, 1 / 3.0, 1 / 3.0, -27.0 / 96.0),
			new GaussPoint3D(0.577350269189626, 0.6, 0.2, 25.0 / 96.0),
			new GaussPoint3D(0.577350269189626, 0.2, 0.6, 25.0 / 96.0),
			new GaussPoint3D(0.577350269189626, 0.2, 0.0, 25.0 / 96.0));

		/// <summary>
		/// 8-integration point rule for order 3 for all axes.
		/// </summary>
		public static readonly WedgeQuadrature Points21 = new WedgeQuadrature(
			new GaussPoint3D(-Math.Sqrt(3.0 / 5.0), 1 / 3.0, 1 / 3.0, 5.0 / 9.0 * 9.0 / 80.0),
			new GaussPoint3D(-Math.Sqrt(3.0 / 5.0), Math.Sqrt(3.0 / 5.0), Math.Sqrt(3.0 / 5.0), 5.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(-Math.Sqrt(3.0 / 5.0), 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, Math.Sqrt(3.0 / 5.0), 5.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(-Math.Sqrt(3.0 / 5.0), Math.Sqrt(3.0 / 5.0), 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, 5.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(-Math.Sqrt(3.0 / 5.0), (6 + Math.Sqrt(15)) / 24.0, (6 + Math.Sqrt(15)) / 24.0, 5.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(-Math.Sqrt(3.0 / 5.0), 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, (6 + Math.Sqrt(15)) / 24.0, 5.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(-Math.Sqrt(3.0 / 5.0), (6 + Math.Sqrt(15)) / 24.0, 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, 5.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0),

			new GaussPoint3D(0, 1 / 3.0, 1 / 3.0, 5.0 / 9.0 * 9.0 / 80.0),
			new GaussPoint3D(0, Math.Sqrt(3.0 / 5.0), Math.Sqrt(3.0 / 5.0), 8.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(0, 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, Math.Sqrt(3.0 / 5.0), 8.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(0, Math.Sqrt(3.0 / 5.0), 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, 8.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(0, (6 + Math.Sqrt(15)) / 24.0, (6 + Math.Sqrt(15)) / 24.0, 8.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(0, 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, (6 + Math.Sqrt(15)) / 24.0, 8.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(0, (6 + Math.Sqrt(15)) / 24.0, 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, 8.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0),

			new GaussPoint3D(Math.Sqrt(3.0 / 5.0), 1 / 3.0, 1 / 3.0, 5.0 / 9.0 * 9.0 / 80.0),
			new GaussPoint3D(Math.Sqrt(3.0 / 5.0), Math.Sqrt(3.0 / 5.0), Math.Sqrt(3.0 / 5.0), 5.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(Math.Sqrt(3.0 / 5.0), 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, Math.Sqrt(3.0 / 5.0), 5.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(Math.Sqrt(3.0 / 5.0), Math.Sqrt(3.0 / 5.0), 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, 5.0 / 9.0 * (155 + Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(Math.Sqrt(3.0 / 5.0), (6 + Math.Sqrt(15)) / 24.0, (6 + Math.Sqrt(15)) / 24.0, 5.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(Math.Sqrt(3.0 / 5.0), 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, (6 + Math.Sqrt(15)) / 24.0, 5.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0),
			new GaussPoint3D(Math.Sqrt(3.0 / 5.0), (6 + Math.Sqrt(15)) / 24.0, 1 - 2 * (6 + Math.Sqrt(15)) / 24.0, 5.0 / 9.0 * (155 - Math.Sqrt(15)) / 2400.0)
		);

		private WedgeQuadrature(params GaussPoint3D[] points)
	    {
		    this.IntegrationPoints = new List<GaussPoint3D>(points);
	    }

	    public IReadOnlyList<GaussPoint3D> IntegrationPoints { get; }
	}
}
