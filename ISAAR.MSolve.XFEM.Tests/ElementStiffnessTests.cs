using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Materials;
using Xunit;

namespace ISAAR.MSolve.XFEM.Tests
{
    public static class ElementStiffnessTests
    {
        private const double E = 2e6;
        private const double v = 0.3;
        private static XNode[] Nodes
        {
            get => new XNode[]
            {
                new XNode(0, 20.0, 0.0),
                new XNode(1, 40.0, 0.0),
                new XNode(2, 40.0, 20.0),
                new XNode(3, 20.0, 20.0),

                new XNode(4, 20.0, -40.0),
                new XNode(5, 40.0, -40.0),
                new XNode(6, 40.0, -20.0),
                new XNode(7, 20.0, -20.0)
            };
        }

        [Fact]
        public static void TestElementWithStrongDiscontinuity()
        {
            double equalityTolerance = 1E-13;
            Matrix expectedKss = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                { +1.154, +0.481, -0.769, +0.096, -0.577, -0.481, +0.192, -0.096 },
                { +0.481, +1.154, -0.096, +0.192, -0.481, -0.577, +0.096, -0.769 },
                { -0.769, -0.096, +1.154, -0.481, +0.192, +0.096, -0.577, +0.481 },
                { +0.096, +0.192, -0.481, +1.154, -0.096, -0.769, +0.481, -0.577 },
                { -0.577, -0.481, +0.192, -0.096, +1.154, +0.481, -0.769, +0.096 },
                { -0.481, -0.577, +0.096, -0.769, +0.481, +1.154, -0.096, +0.192 },
                { +0.192, +0.096, -0.577, +0.481, -0.769, -0.096, +1.154, -0.481 },
                { -0.096, -0.769, +0.481, -0.577, +0.096, +0.192, -0.481, +1.154 }
            });

            Matrix expectedKse = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                { +0.962, +0.240, +0.769, +0.144, +0.577, +0.433, +0.385, -0.048 },
                { +0.240, +0.481, +0.337, -0.192, +0.529, +0.577, +0.048, -0.096 },
                { -0.769, +0.144, -0.962, +0.240, -0.385, -0.048, -0.577, +0.433 },
                { +0.337, +0.192, +0.240, -0.481, +0.048, +0.096, +0.529, -0.577 },
                { -0.577, -0.433, -0.385, +0.048, -0.962, -0.240, -0.769, -0.144 },
                { -0.529, -0.577, -0.048, +0.096, -0.240, -0.481, -0.337, +0.192 },
                { +0.385, +0.048, +0.577, -0.433, +0.769, -0.144, +0.962, -0.240 },
                { -0.048, -0.096, -0.529, +0.577, -0.337, -0.192, -0.240, +0.481 }
            });

            Matrix expectedKes = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                { +0.962, +0.240, -0.769, +0.337, -0.577, -0.529, +0.385, -0.048 },
                { +0.240, +0.481, +0.144, +0.192, -0.433, -0.577, +0.048, -0.096 },
                { +0.769, +0.337, -0.962, +0.240, -0.385, -0.048, +0.577, -0.529 },
                { +0.144, -0.192, +0.240, -0.481, +0.048, +0.096, -0.433, +0.577 },
                { +0.577, +0.529, -0.385, +0.048, -0.962, -0.240, +0.769, -0.337 },
                { +0.433, +0.577, -0.048, +0.096, -0.240, -0.481, -0.144, -0.192 },
                { +0.385, +0.048, -0.577, +0.529, -0.769, -0.337, +0.962, -0.240 },
                { -0.048, -0.096, +0.433, -0.577, -0.144, +0.192, -0.240, +0.481 }
            });

            Matrix expectedKee = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                { +1.923, +0.481, +0.000, +0.000, +0.000, +0.000, +0.769, -0.096 },
                { +0.481, +0.962, +0.000, +0.000, +0.000, +0.000, +0.096, -0.192 },
                { +0.000, +0.000, +1.923, -0.481, +0.769, +0.096, +0.000, +0.000 },
                { +0.000, +0.000, -0.481, +0.962, -0.096, -0.192, +0.000, +0.000 },
                { +0.000, +0.000, +0.769, -0.096, +1.923, +0.481, +0.000, +0.000 },
                { +0.000, +0.000, +0.096, -0.192, +0.481, +0.962, +0.000, +0.000 },
                { +0.769, +0.096, +0.000, +0.000, +0.000, +0.000, +1.923, -0.481 },
                { -0.096, -0.192, +0.000, +0.000, +0.000, +0.000, -0.481, +0.962 }
            });

            var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);

            var integrationStrategy = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.GetQuadratureWithOrder(2, 2)),
                new SimpleIntegration2D());
            //var integrationStrategy = new IntegrationForCrackPropagation2D(GaussLegendre2D.GetQuadratureWithOrder(2, 2),
            //  new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.GetQuadratureWithOrder(2, 2)));
            var factory = new XContinuumElement2DFactory(integrationStrategy, integrationStrategy, material);
            var bodyElement = factory.CreateElement(CellType.Quad4, new XNode[] { Nodes[0], Nodes[1], Nodes[2], Nodes[3] });
            var blendingElement = factory.CreateElement(CellType.Quad4, new XNode[] { Nodes[7], Nodes[6], Nodes[1], Nodes[0] });
            var tipElement = factory.CreateElement(CellType.Quad4, new XNode[] { Nodes[4], Nodes[5], Nodes[6], Nodes[7] });
            var boundary = new Rectangular2DBoundary(20.0, 40.0, -40.0, 20.0);
            var mesh = new SimpleMesh2D<XNode, XContinuumElement2D>(Nodes,
                new XContinuumElement2D[] { bodyElement, blendingElement, tipElement }, boundary);

            var crack = new BasicExplicitCrack2D();
            crack.Mesh = mesh;
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);

            var point1 = new CartesianPoint(30.0, 20.0);
            var point2 = new CartesianPoint(30.0, -30.0);
            crack.InitializeGeometry(point1, point2);
            crack.UpdateEnrichments();

            Matrix Kss = bodyElement.BuildStandardStiffnessMatrix();
            (Matrix Kee1, Matrix Kes) = bodyElement.BuildEnrichedStiffnessMatricesLower();
            (Matrix Kee2, Matrix Kse) = bodyElement.BuildEnrichedStiffnessMatricesUpper();

            Assert.True(expectedKss.Equals(Kss.DoToAllEntries(x => 1E6 * Math.Round(x * 1E-6, 3)), equalityTolerance));
            Assert.True(expectedKee.Equals(Kee1.DoToAllEntries(x => 1E6 * Math.Round(x * 1E-6, 3)), equalityTolerance));
            Assert.True(expectedKee.Equals(Kee2.DoToAllEntries(x => 1E6 * Math.Round(x * 1E-6, 3)), equalityTolerance));
            Assert.True(expectedKes.Equals(Kes.DoToAllEntries(x => 1E6 * Math.Round(x * 1E-6, 3)), equalityTolerance));
            Assert.True(expectedKse.Equals(Kse.DoToAllEntries(x => 1E6 * Math.Round(x * 1E-6, 3)), equalityTolerance));
        }
    }
}
