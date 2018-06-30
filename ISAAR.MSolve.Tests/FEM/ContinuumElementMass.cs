using System;
using System.Collections.Generic;
using System.Text;
using Xunit;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Integration.Points;
using ISAAR.MSolve.FEM.Integration.Quadratures;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ContinuumElementMass
    {
        private static double thickness = 1.0;

        private static readonly ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
        {
            YoungModulus = 2e6,
            PoissonRatio = 0.3
        };

        private static readonly DynamicMaterial dynamicMaterial = new DynamicMaterial(78.5, 0, 0);

        private static readonly ContinuumElement2DFactory factory = 
            new ContinuumElement2DFactory(thickness, material, dynamicMaterial);


        public static readonly IReadOnlyList<Node2D> tri3NodeSet1 = new Node2D[]
        {
            new Node2D(0, 0.0,  0.0),
            new Node2D(1, 1.0,  0.0),
            new Node2D(2, 0.0, 1.0)
        };

        private static readonly IReadOnlyList<Node2D> tri3NodeSet2 = new Node2D[]
        {
            new Node2D(0, 2.0,  1.0),
            new Node2D(1, 6.0,  3.0),
            new Node2D(2, 1.8, 12.0)
        };

        [Fact]
        private static void TestTri3ConsistentMass()
        {
            TestTri3ConsistentMassParametric(tri3NodeSet1, GaussQuadratureForTriangles.Order1Point1, true);
            TestTri3ConsistentMassParametric(tri3NodeSet2, GaussQuadratureForTriangles.Order2Points3, false);
        }

        private static void TestTri3ConsistentMassParametric(IReadOnlyList<Node2D> nodeSet, IQuadrature2D quadratureForMass, 
            bool reducedQuadrature)
        {
            var materialsAtGaussPoints = new Dictionary<GaussPoint2D, ElasticMaterial2D>();
            foreach (GaussPoint2D gaussPoint in quadratureForMass.IntegrationPoints)
            {
                materialsAtGaussPoints[gaussPoint] = material.Clone();
            }
            var tri3 = new ContinuumElement2D(thickness, nodeSet, InterpolationTri3.UniqueInstance, 
                GaussQuadratureForTriangles.Order1Point1, quadratureForMass, 
                ExtrapolationGaussTriangular1Point.UniqueInstance,
                materialsAtGaussPoints, dynamicMaterial);
            IMatrix2D M = tri3.BuildConsistentMassMatrix();

            // Reference: http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf, (eq 31.27) 
            Matrix2D expectedM = new Matrix2D(new double[,]
            {
                { 2, 0, 1, 0, 1, 0 },
                { 0, 2, 0, 1, 0, 1 },
                { 1, 0, 2, 0, 1, 0 },
                { 0, 1, 0, 2, 0, 1 },
                { 1, 0, 1, 0, 2, 0 },
                { 0, 1, 0, 1, 0, 2 }
            });
            double area = CalcTriangleArea(nodeSet);
            double scalar = dynamicMaterial.Density * thickness * area / 12.0;
            expectedM.Scale(scalar);

            if (reducedQuadrature) Assert.False(Utilities.AreMatricesEqual(M, expectedM, 1e-10));
            else Assert.True(Utilities.AreMatricesEqual(M, expectedM, 1e-10));

        }

        private static double CalcTriangleArea(IReadOnlyList<Node2D> nodes)
        {
            return (nodes[0].X * (nodes[1].Y - nodes[2].Y)
                + nodes[1].X * (nodes[2].Y - nodes[0].Y)
                + nodes[2].X * (nodes[0].Y - nodes[1].Y)) / 2.0;
        }
    }
}
