using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.FEM.Elements
{
	/// <summary>
	/// Creates isoparametric continuum 3D elements. Abstracts the interpolations, integrations,
	/// extrapolations and any other strategies that differentiate the elements(e.g. Hexa8, Hexa20)
	/// It is also very convenient when the material properties are the same throughout the whole domain or a region.
	/// Authors: Dimitris Tsapetis
	/// </summary>
    public class ContinuumElement3DFactory
	{
		private static readonly IReadOnlyDictionary<CellType3D, IGaussPointExtrapolation3D> extrapolations;
		private static readonly IReadOnlyDictionary<CellType3D, IQuadrature3D> integrationsForStiffness;
		private static readonly IReadOnlyDictionary<CellType3D, IQuadrature3D> integrationsForMass;
		private static readonly IReadOnlyDictionary<CellType3D, IIsoparametricInterpolation3D> interpolations;

		private ElasticMaterial commonMaterial;
		private DynamicMaterial commonDynamicProperties;

		static ContinuumElement3DFactory()
		{
			var interpolations=new Dictionary<CellType3D,IIsoparametricInterpolation3D>();
			var integrationsForStiffness = new Dictionary<CellType3D, IQuadrature3D>();
			var integrationsForMass= new Dictionary<CellType3D,IQuadrature3D>();
			var extrapolations = new Dictionary<CellType3D, IGaussPointExtrapolation3D>();

			// Tet4
			// TODO: implementations for Tet4
			interpolations.Add(CellType3D.Tet4, InterpolationTet4.UniqueInstance);

			// Tet10
			// TODO: implementations for Tet10
			interpolations.Add(CellType3D.Tet10, InterpolationTet10.UniqueInstance);

			// Hexa8
			interpolations.Add(CellType3D.Hexa8, InterpolationHexa8.UniqueInstance);
			integrationsForStiffness.Add(CellType3D.Hexa8, GaussLegendre3D.Order2x2x2);
			integrationsForMass.Add(CellType3D.Hexa8, GaussLegendre3D.Order2x2x2);
			extrapolations.Add(CellType3D.Hexa8, ExtrapolationGaussLegendre2x2x2.UniqueInstance);

			// Hexa20
			// TODO: implementations for Hexa20
			interpolations.Add(CellType3D.Hexa20, InterpolationHexa20.UniqueInstance);

			// Hexa27
			// TODO: implementations for Hexa27
			interpolations.Add(CellType3D.Hexa27, InterpolationHexa27.UniqueInstance);

			// Wedge6
			// TODO: implementations for Wedge6
			interpolations.Add(CellType3D.Wedge6, InterpolationWedge6.UniqueInstance);

			// Wedge15
			// TODO: implementations for Wedge15
			interpolations.Add(CellType3D.Wedge15, InterpolationWedge15.UniqueInstance);

			// Wedge18
			// TODO: implementations for Wedge18
			interpolations.Add(CellType3D.Wedge18, InterpolationWedge18.UniqueInstance);

			// Pyra5
			// TODO: implementations for Pyra5
			interpolations.Add(CellType3D.Pyra5, InterpolationPyra5.UniqueInstance);

			// Pyra13
			// TODO: implementations for Pyra13
			interpolations.Add(CellType3D.Pyra13, InterpolationPyra13.UniqueInstance);

			// Pyra14
			// TODO: implementations for Pyra14
			interpolations.Add(CellType3D.Pyra14, InterpolationPyra14.UniqueInstance);
		}
	}
}
