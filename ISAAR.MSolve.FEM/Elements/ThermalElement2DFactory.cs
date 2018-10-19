using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.FEM.Elements
{
    public class ThermalElement2DFactory
    {
        private static readonly IReadOnlyDictionary<CellType2D, IGaussPointExtrapolation2D> extrapolations;
        private static readonly IReadOnlyDictionary<CellType2D, IQuadrature2D> integrationsForStiffness;
        private static readonly IReadOnlyDictionary<CellType2D, IQuadrature2D> integrationsForMass;
        private static readonly IReadOnlyDictionary<CellType2D, IIsoparametricInterpolation2D> interpolations;

        private ThermalMaterial commonMaterial;
        private double commonThickness;

        static ThermalElement2DFactory()
        {
            // Mass integrations require as many Gauss points as there are nodes, in order for the consistent mass matrix to be
            // of full rank (and symmetric positive definite)

            // Collections' declarations
            var interpolations = new Dictionary<CellType2D, IIsoparametricInterpolation2D>();
            var integrationsForStiffness = new Dictionary<CellType2D, IQuadrature2D>();
            var integrationsForMass = new Dictionary<CellType2D, IQuadrature2D>();
            var extrapolations = new Dictionary<CellType2D, IGaussPointExtrapolation2D>();

            // Quad4
            interpolations.Add(CellType2D.Quad4, InterpolationQuad4.UniqueInstance);
            integrationsForStiffness.Add(CellType2D.Quad4, GaussLegendre2D.GetQuadratureWithOrder(2, 2));
            integrationsForMass.Add(CellType2D.Quad4, GaussLegendre2D.GetQuadratureWithOrder(2, 2));
            extrapolations.Add(CellType2D.Quad4, ExtrapolationGaussLegendre2x2.UniqueInstance);

            // Quad8
            interpolations.Add(CellType2D.Quad8, InterpolationQuad8.UniqueInstance);
            integrationsForStiffness.Add(CellType2D.Quad8, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
            integrationsForMass.Add(CellType2D.Quad8, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
            extrapolations.Add(CellType2D.Quad8, ExtrapolationGaussLegendre3x3.UniqueInstance);

            // Quad9
            interpolations.Add(CellType2D.Quad9, InterpolationQuad9.UniqueInstance);
            integrationsForStiffness.Add(CellType2D.Quad9, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
            integrationsForMass.Add(CellType2D.Quad9, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
            extrapolations.Add(CellType2D.Quad9, ExtrapolationGaussLegendre3x3.UniqueInstance);

            // Tri3
            interpolations.Add(CellType2D.Tri3, InterpolationTri3.UniqueInstance);
            integrationsForStiffness.Add(CellType2D.Tri3, TriangleQuadratureSymmetricGaussian.Order1Point1);
            integrationsForMass.Add(CellType2D.Tri3, TriangleQuadratureSymmetricGaussian.Order2Points3);
            extrapolations.Add(CellType2D.Tri3, ExtrapolationGaussTriangular1Point.UniqueInstance);

            // Tri 6
            interpolations.Add(CellType2D.Tri6, InterpolationTri6.UniqueInstance);
            // see https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf, p. 24-13, paragraph "options"
            integrationsForStiffness.Add(CellType2D.Tri6, TriangleQuadratureSymmetricGaussian.Order2Points3);
            integrationsForMass.Add(CellType2D.Tri6, TriangleQuadratureSymmetricGaussian.Order4Points6);
            extrapolations.Add(CellType2D.Tri6, ExtrapolationGaussTriangular3Points.UniqueInstance);

            // Static field assignments
            ThermalElement2DFactory.interpolations = interpolations;
            ThermalElement2DFactory.integrationsForStiffness = integrationsForStiffness;
            ThermalElement2DFactory.integrationsForMass = integrationsForMass;
            ThermalElement2DFactory.extrapolations = extrapolations;
        }

        public ThermalElement2DFactory(double commonThickness, ThermalMaterial commonMaterial)
        {
            this.commonThickness = commonThickness;
            this.commonMaterial = commonMaterial;
        }

        public ThermalElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes)
        {
            return new ThermalElement2D(commonThickness, nodes, interpolations[cellType],
               integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType],
               commonMaterial);
        }

        //public ThermalElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes)
        //{
        //    return CreateElement(cellType, nodes, commonThickness, commonMaterial);
        //}

        //public ThermalElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes, double thickness, ThermalMaterial material)
        //{
        //    var materialsAtGaussPoints = new Dictionary<GaussPoint2D, ThermalMaterial>();
        //    foreach (GaussPoint2D gaussPoint in integrationsForStiffness[cellType].IntegrationPoints)
        //    {
        //        materialsAtGaussPoints[gaussPoint] = material.Clone();
        //    }
        //    return CreateElement(cellType, nodes, thickness, materialsAtGaussPoints);
        //}

        //public ThermalElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes, double thickness,
        //    Dictionary<GaussPoint2D, ThermalMaterial> materialsAtGaussPoints)
        //{
        //    //TODO: check if nodes - interpolation and Gauss points - materials match
        //    return new ThermalElement2D(thickness, nodes, interpolations[cellType],
        //        integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType],
        //        materialsAtGaussPoints);
        //}
    }
}
