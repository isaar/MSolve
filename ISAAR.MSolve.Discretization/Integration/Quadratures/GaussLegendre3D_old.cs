using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
	/// <summary>
	/// Enum class with the 3D Gauss-Legendre integration rules of varying orders. Tensor product of 
	/// <see cref="GaussLegendre1D_old"/> rules along each axis.
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public sealed class GaussLegendre3D_old : IQuadrature3D
	{
		public static readonly GaussLegendre3D_old Order1x1x1 = new GaussLegendre3D_old(GaussLegendre1D_old.Order1, GaussLegendre1D_old.Order1, GaussLegendre1D_old.Order1);
		public static readonly GaussLegendre3D_old Order2x2x2 = new GaussLegendre3D_old(GaussLegendre1D_old.Order2, GaussLegendre1D_old.Order2, GaussLegendre1D_old.Order2);
		public static readonly GaussLegendre3D_old Order3x3x3 = new GaussLegendre3D_old(GaussLegendre1D_old.Order3, GaussLegendre1D_old.Order3, GaussLegendre1D_old.Order3);
		public static readonly GaussLegendre3D_old Order4x4x4 = new GaussLegendre3D_old(GaussLegendre1D_old.Order4, GaussLegendre1D_old.Order4, GaussLegendre1D_old.Order4);
		public static readonly GaussLegendre3D_old Order5x5x5 = new GaussLegendre3D_old(GaussLegendre1D_old.Order5, GaussLegendre1D_old.Order5, GaussLegendre1D_old.Order5);
		public static readonly GaussLegendre3D_old Order6x6x6 = new GaussLegendre3D_old(GaussLegendre1D_old.Order6, GaussLegendre1D_old.Order6, GaussLegendre1D_old.Order6);
		public static readonly GaussLegendre3D_old Order7x7x7 = new GaussLegendre3D_old(GaussLegendre1D_old.Order7, GaussLegendre1D_old.Order7, GaussLegendre1D_old.Order7);
		public static readonly GaussLegendre3D_old Order8x8x8 = new GaussLegendre3D_old(GaussLegendre1D_old.Order8, GaussLegendre1D_old.Order8, GaussLegendre1D_old.Order8);
		public static readonly GaussLegendre3D_old Order9x9x9 = new GaussLegendre3D_old(GaussLegendre1D_old.Order9, GaussLegendre1D_old.Order9, GaussLegendre1D_old.Order9);
		public static readonly GaussLegendre3D_old Order10x10x10 = new GaussLegendre3D_old(GaussLegendre1D_old.Order10, GaussLegendre1D_old.Order10, GaussLegendre1D_old.Order10);


		private GaussLegendre3D_old(GaussLegendre1D_old ruleXi, GaussLegendre1D_old ruleEta, GaussLegendre1D_old ruleZeta)
		{
			// Combine the integration rules of each axis: 
			// WARNING: Do not change their order (Xi major, Eta minor). Other classes, such as ExtrapolationGaussLegendre2x2 
			//          depend on it.
			var points3D = new List<GaussPoint3D>();
			foreach (var pointZeta in ruleZeta.IntegrationPoints)
			{
				foreach (var pointEta in ruleEta.IntegrationPoints)
				{
					foreach (var pointXi in ruleXi.IntegrationPoints)
					{
						points3D.Add(new GaussPoint3D(pointXi.Xi, pointEta.Xi, pointZeta.Xi,
							pointXi.Weight * pointEta.Weight * pointZeta.Weight));
					}
				}
			}
			
			this.IntegrationPoints = points3D;
		}

		/// <summary>
		/// The integration points are sorted based on an order strictly defined for each quadrature.
		/// </summary>
		public IReadOnlyList<GaussPoint3D> IntegrationPoints { get; }
	}
}
