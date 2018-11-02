using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.FEM.Elements
{
    public class ThermalElement3DFactory
    {
        private static readonly IReadOnlyDictionary<CellType3D, IGaussPointExtrapolation3D> extrapolations;
        private static readonly IReadOnlyDictionary<CellType3D, IQuadrature3D> integrationsForStiffness;
        private static readonly IReadOnlyDictionary<CellType3D, IQuadrature3D> integrationsForMass;
        private static readonly IReadOnlyDictionary<CellType3D, IIsoparametricInterpolation3D> interpolations;

        private ThermalMaterial commonMaterial;

        static ThermalElement3DFactory()
        {
            var interpolations = new Dictionary<CellType3D, IIsoparametricInterpolation3D>();
            var integrationsForStiffness = new Dictionary<CellType3D, IQuadrature3D>();
            var integrationsForMass = new Dictionary<CellType3D, IQuadrature3D>();
            var extrapolations = new Dictionary<CellType3D, IGaussPointExtrapolation3D>();

            // Tet4
            // TODO: implementations for Tet4
            interpolations.Add(CellType3D.Tet4, InterpolationTet4.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Tet4, TetrahedronQuadrature.Order1Point1);
            integrationsForMass.Add(CellType3D.Tet4, TetrahedronQuadrature.Order2Points4);
            extrapolations.Add(CellType3D.Tet4, null);

            // Tet10
            // TODO: implementations for Tet10
            interpolations.Add(CellType3D.Tet10, InterpolationTet10.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Tet10, TetrahedronQuadrature.Order2Points4);
            integrationsForMass.Add(CellType3D.Tet10, TetrahedronQuadrature.Order5Points15);
            extrapolations.Add(CellType3D.Tet10, null);

            // Hexa8
            interpolations.Add(CellType3D.Hexa8, InterpolationHexa8.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
            integrationsForMass.Add(CellType3D.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));
            extrapolations.Add(CellType3D.Hexa8, ExtrapolationGaussLegendre2x2x2.UniqueInstance);

            // Hexa20
            // TODO: extrapolations for Hexa20
            interpolations.Add(CellType3D.Hexa20, InterpolationHexa20.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Hexa20, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));
            integrationsForMass.Add(CellType3D.Hexa20, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));
            extrapolations.Add(CellType3D.Hexa20, null);

            // Hexa27
            // TODO: extrapolations for Hexa27
            interpolations.Add(CellType3D.Hexa27, InterpolationHexa27.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Hexa27, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));
            integrationsForMass.Add(CellType3D.Hexa27, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));
            extrapolations.Add(CellType3D.Hexa27, null);

            // Wedge6
            // TODO: implementations for Wedge6
            interpolations.Add(CellType3D.Wedge6, InterpolationWedge6.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Wedge6, WedgeQuadrature.Points6);
            integrationsForMass.Add(CellType3D.Wedge6, WedgeQuadrature.Points8);
            extrapolations.Add(CellType3D.Wedge6, null);

            // Wedge15
            // TODO: implementations for Wedge15
            interpolations.Add(CellType3D.Wedge15, InterpolationWedge15.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Wedge15, WedgeQuadrature.Points8);
            integrationsForMass.Add(CellType3D.Wedge15, WedgeQuadrature.Points21);
            extrapolations.Add(CellType3D.Wedge15, null);

            // Wedge18
            // TODO: implementations for Wedge18
            interpolations.Add(CellType3D.Wedge18, InterpolationWedge18.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Wedge18, WedgeQuadrature.Points8);
            integrationsForMass.Add(CellType3D.Wedge18, WedgeQuadrature.Points21);
            extrapolations.Add(CellType3D.Wedge18, null);

            // Pyra5
            // TODO: implementations for Pyra5
            interpolations.Add(CellType3D.Pyra5, InterpolationPyra5.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Pyra5, PyramidQuadrature.Points5);
            integrationsForMass.Add(CellType3D.Pyra5, PyramidQuadrature.Points5);
            extrapolations.Add(CellType3D.Pyra5, null);

            // Pyra13
            // TODO: implementations for Pyra13
            interpolations.Add(CellType3D.Pyra13, InterpolationPyra13.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Pyra13, PyramidQuadrature.Points6);
            integrationsForMass.Add(CellType3D.Pyra13, PyramidQuadrature.Points6);
            extrapolations.Add(CellType3D.Pyra13, null);

            // Pyra14
            // TODO: implementations for Pyra14
            interpolations.Add(CellType3D.Pyra14, InterpolationPyra14.UniqueInstance);
            integrationsForStiffness.Add(CellType3D.Pyra14, PyramidQuadrature.Points6);
            integrationsForMass.Add(CellType3D.Pyra14, PyramidQuadrature.Points6);
            extrapolations.Add(CellType3D.Pyra14, null);

            ThermalElement3DFactory.interpolations = interpolations;
            ThermalElement3DFactory.integrationsForStiffness = integrationsForStiffness;
            ThermalElement3DFactory.integrationsForMass = integrationsForMass;
            ThermalElement3DFactory.extrapolations = extrapolations;
        }

        public ThermalElement3DFactory(ThermalMaterial commonMaterial)
        {
            this.commonMaterial = commonMaterial;
        }

        public ThermalElement3D CreateElement(CellType3D cellType, IReadOnlyList<Node3D> nodes)
        {
            return new ThermalElement3D(nodes, interpolations[cellType],
                integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType], commonMaterial);
        }
    }
}
