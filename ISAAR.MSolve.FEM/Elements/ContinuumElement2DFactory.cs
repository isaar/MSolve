using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Meshes;
using ISAAR.MSolve.Materials;

//TODO: Materials should be passed in the constructor, not each method
//TODO: not sure about the quadratures for mass. At least as many Gauss points as there are nodes are necessary 
//      (for positive definite mass matrices), but perhaps I need more for exact integration of the polynomials
namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Creates isoparametric continuum elements with uniform thickness. Abstracts the interpolations, integrations,
    /// extrapolations and any other strategies that differentiate the elements (e.g. Quad4 from Tri6). It is also very 
    /// convenient when the thickness and material properties are the same throughout the whole domain or a region.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ContinuumElement2DFactory
    {
        private static readonly IReadOnlyDictionary<CellType2D, IGaussPointExtrapolation2D> extrapolations;
        private static readonly IReadOnlyDictionary<CellType2D, IQuadrature2D> integrationsForStiffness;
        private static readonly IReadOnlyDictionary<CellType2D, IQuadrature2D> integrationsForMass;
        private static readonly IReadOnlyDictionary<CellType2D, IIsoparametricInterpolation2D> interpolations;

        private ElasticMaterial2D commonMaterial;
        private DynamicMaterial commonDynamicProperties;
        private double commonThickness;

        static ContinuumElement2DFactory()
        {
            // Interpolations
            var interpolations = new Dictionary<CellType2D, IIsoparametricInterpolation2D>();
            interpolations.Add(CellType2D.Quad4, InterpolationQuad4.UniqueInstance);
            interpolations.Add(CellType2D.Tri3, InterpolationTri3.UniqueInstance);
            ContinuumElement2DFactory.interpolations = interpolations;

            // Integration rules for stiffness
            var integrationsForStiffness = new Dictionary<CellType2D, IQuadrature2D>();
            integrationsForStiffness.Add(CellType2D.Quad4, GaussLegendre2D.Order2x2);
            integrationsForStiffness.Add(CellType2D.Tri3, GaussQuadratureForTriangles.Order1Point1);
            ContinuumElement2DFactory.integrationsForStiffness = integrationsForStiffness;

            // Integration rules for mass
            var integrationsForMass = new Dictionary<CellType2D, IQuadrature2D>();
            integrationsForMass.Add(CellType2D.Quad4, GaussLegendre2D.Order2x2);
            integrationsForMass.Add(CellType2D.Tri3, GaussQuadratureForTriangles.Order2Points3);
            ContinuumElement2DFactory.integrationsForMass = integrationsForMass;

            // Extrapolations from integration points
            var extrapolations = new Dictionary<CellType2D, IGaussPointExtrapolation2D>();
            extrapolations.Add(CellType2D.Quad4, ExtrapolationGaussLegendre2x2.UniqueInstance);
            extrapolations.Add(CellType2D.Tri3, ExtrapolationGaussTriangular1Point.UniqueInstance);
            ContinuumElement2DFactory.extrapolations = extrapolations;
        }

        public ContinuumElement2DFactory(double commonThickness, ElasticMaterial2D commonMaterial, 
            DynamicMaterial commonDynamicProperties)
        {
            this.commonThickness = commonThickness;
            this.commonMaterial = commonMaterial;
            this.commonDynamicProperties = commonDynamicProperties;
        }

        public ContinuumElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes)
        {
            return CreateElement(cellType, nodes, commonThickness, commonMaterial, commonDynamicProperties);
        }

        public ContinuumElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes, double thickness,
            ElasticMaterial2D material, DynamicMaterial dynamicProperties)
        {
            var materialsAtGaussPoints = new Dictionary<GaussPoint2D, ElasticMaterial2D>();
            foreach (GaussPoint2D gaussPoint in integrationsForStiffness[cellType].IntegrationPoints)
            {
                materialsAtGaussPoints[gaussPoint] = material.Clone();
            }
            return CreateElement(cellType, nodes, thickness, materialsAtGaussPoints, dynamicProperties);
        }

        public ContinuumElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes, double thickness,
            Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
        {
            //TODO: check if nodes - interpolation and Gauss points - materials match
            return new ContinuumElement2D(thickness, nodes, interpolations[cellType],
                integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType],
                materialsAtGaussPoints, dynamicProperties);
        }

        #region obsolete
        public ContinuumElement2D CreateQuad4(IReadOnlyList<Node2D> nodes)
        {
            return CreateQuad4(nodes, commonThickness, commonMaterial, commonDynamicProperties);
        }

        public ContinuumElement2D CreateQuad4(IReadOnlyList<Node2D> nodes, double thickness, ElasticMaterial2D material,
            DynamicMaterial dynamicProperties)
        {
            var materialsAtGaussPoints = new Dictionary<GaussPoint2D, ElasticMaterial2D>();
            foreach (GaussPoint2D gaussPoint in GaussLegendre2D.Order2x2.IntegrationPoints)
            {
                materialsAtGaussPoints[gaussPoint] = material.Clone();
            }
            return CreateQuad4(nodes, thickness, materialsAtGaussPoints, dynamicProperties);
        }

        public ContinuumElement2D CreateQuad4(IReadOnlyList<Node2D> nodes, double thickness,
            Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
        {
            //TODO: check if nodes - interpolation and Gauss points - materials match
            return new ContinuumElement2D(thickness, nodes, InterpolationQuad4.UniqueInstance,
                GaussLegendre2D.Order2x2, GaussLegendre2D.Order2x2,
                ExtrapolationGaussLegendre2x2.UniqueInstance, materialsAtGaussPoints, dynamicProperties);
        }

        public ContinuumElement2D CreateTri3(IReadOnlyList<Node2D> nodes)
        {
            return CreateTri3(nodes, commonThickness, commonMaterial, commonDynamicProperties);
        }

        public ContinuumElement2D CreateTri3(IReadOnlyList<Node2D> nodes, double thickness, ElasticMaterial2D material,
            DynamicMaterial dynamicProperties)
        {
            var materialsAtGaussPoints = new Dictionary<GaussPoint2D, ElasticMaterial2D>();
            foreach (GaussPoint2D gaussPoint in GaussQuadratureForTriangles.Order1Point1.IntegrationPoints)
            {
                materialsAtGaussPoints[gaussPoint] = material.Clone();
            }
            return CreateTri3(nodes, thickness, materialsAtGaussPoints, dynamicProperties);
        }

        public ContinuumElement2D CreateTri3(IReadOnlyList<Node2D> nodes, double thickness,
            Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
        {
            //TODO: check if nodes - interpolation and Gauss points - materials match
            return new ContinuumElement2D(thickness, nodes, InterpolationTri3.UniqueInstance, 
                GaussQuadratureForTriangles.Order1Point1, GaussQuadratureForTriangles.Order2Points3,
                ExtrapolationGaussTriangular1Point.UniqueInstance, materialsAtGaussPoints, dynamicProperties);
        }
        #endregion
    }
}
