using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Integration.Points;
using ISAAR.MSolve.FEM.Integration.Quadratures;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Materials;

//TODO: Materials should be passed in the constructor, not each method
//TODO: not sure about the quadratures for mass. At least as many Gauss points as there are nodes are necessary 
//      (for positive definite mass matrices), but perhaps I need more for exact integration of the polynomials
namespace ISAAR.MSolve.FEM.Elements
{
    public class ContinuumElement2DFactory
    {
        private double commonThickness;
        private ElasticMaterial2D commonMaterial;
        private DynamicMaterial commonDynamicProperties;

        public ContinuumElement2DFactory(double commonThickness, ElasticMaterial2D commonMaterial, 
            DynamicMaterial commonDynamicProperties)
        {
            this.commonThickness = commonThickness;
            this.commonMaterial = commonMaterial;
            this.commonDynamicProperties = commonDynamicProperties;
        }

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
    }
}
