using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Integration.Points;
using ISAAR.MSolve.FEM.Integration.Quadratures;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.FEM.Elements
{
    public class ContinuumElement2DFactory
    {
        public ContinuumElement2D CreateQuad4(IReadOnlyList<Node2D> nodes, 
            Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints)
        {
            //TODO: check if nodes - interpolation and Gauss points - materials match
            return new ContinuumElement2D(nodes, InterpolationQuad4.UniqueInstance, GaussLegendre2D.Order2x2,
                ExtrapolationGaussLegendre2x2.UniqueInstance, materialsAtGaussPoints);
        }

        public ContinuumElement2D CreateQuad4(IReadOnlyList<Node2D> nodes, ElasticMaterial2D material)
        {
            var materialsAtGaussPoints = new Dictionary<GaussPoint2D, ElasticMaterial2D>();
            foreach (GaussPoint2D gaussPoint in GaussLegendre2D.Order2x2.IntegrationPoints)
            {
                materialsAtGaussPoints[gaussPoint] = material.Clone();
            }
            return CreateQuad4(nodes, materialsAtGaussPoints);
        }

        public ContinuumElement2D CreateTri3(IReadOnlyList<Node2D> nodes,
            Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints)
        {
            //TODO: check if nodes - interpolation and Gauss points - materials match
            return new ContinuumElement2D(nodes, InterpolationTri3.UniqueInstance, GaussQuadratureForTriangles.Order1Point1,
                ExtrapolationGaussTriangular1Point.UniqueInstance, materialsAtGaussPoints);
        }

        public ContinuumElement2D CreateTri3(IReadOnlyList<Node2D> nodes, ElasticMaterial2D material)
        {
            var materialsAtGaussPoints = new Dictionary<GaussPoint2D, ElasticMaterial2D>();
            foreach (GaussPoint2D gaussPoint in GaussQuadratureForTriangles.Order1Point1.IntegrationPoints)
            {
                materialsAtGaussPoints[gaussPoint] = material.Clone();
            }
            return CreateTri3(nodes, materialsAtGaussPoints);
        }
    }
}
