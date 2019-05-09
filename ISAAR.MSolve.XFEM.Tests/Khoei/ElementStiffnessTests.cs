using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration;
using ISAAR.MSolve.XFEM.Materials;
using Xunit;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    /// <summary>
    /// Tests taken from "Extended Finite Element Method: Theory and Applications, Amir R. Khoei, 2015", sections 2.8.3 & 3.8.4.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class ElementStiffnessTests
    {
        private static readonly Func<double, double> round = x => 1E6 * Math.Round(x * 1E-6, 3);

        [Fact]
        public static void TestJointStiffnessMatrixOfCrackBodyElement()
        {
            double equalityTolerance = 1E-13;
            Matrix expectedStiffnessNodeMajor = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                {  1.154,  0.481,  0.962,  0.240, -0.769,  0.096,  0.769,  0.144, -0.577, -0.481,  0.577,  0.433,  0.192, -0.096,  0.385, -0.048 },
                {  0.481,  1.154,  0.240,  0.481, -0.096,  0.192,  0.337, -0.192, -0.481, -0.577,  0.529,  0.577,  0.096, -0.769,  0.048, -0.096 },
                {  0.962,  0.240,  1.923,  0.481, -0.769,  0.337,  0.000,  0.000, -0.577, -0.529,  0.000,  0.000,  0.385, -0.048,  0.769, -0.096 },
                {  0.240,  0.481,  0.481,  0.962,  0.144,  0.192,  0.000,  0.000, -0.433, -0.577,  0.000,  0.000,  0.048, -0.096,  0.096, -0.192 },
                { -0.769, -0.096, -0.769,  0.144,  1.154, -0.481, -0.962,  0.240,  0.192,  0.096, -0.385, -0.048, -0.577,  0.481, -0.577,  0.433 },
                {  0.096,  0.192,  0.337,  0.192, -0.481,  1.154,  0.240, -0.481, -0.096, -0.769,  0.048,  0.096,  0.481, -0.577,  0.529, -0.577 },
                {  0.769,  0.337,  0.000,  0.000, -0.962,  0.240,  1.923, -0.481, -0.385, -0.048,  0.769,  0.096,  0.577, -0.529,  0.000,  0.000 },
                {  0.144, -0.192,  0.000,  0.000,  0.240, -0.481, -0.481,  0.962,  0.048,  0.096, -0.096, -0.192, -0.433,  0.577,  0.000,  0.000 },
                { -0.577, -0.481, -0.577, -0.433,  0.192, -0.096, -0.385,  0.048,  1.154,  0.481, -0.962, -0.240, -0.769,  0.096, -0.769, -0.144 },
                { -0.481, -0.577, -0.529, -0.577,  0.096, -0.769, -0.048,  0.096,  0.481,  1.154, -0.240, -0.481, -0.096,  0.192, -0.337,  0.192 },
                {  0.577,  0.529,  0.000,  0.000, -0.385,  0.048,  0.769, -0.096, -0.962, -0.240,  1.923,  0.481,  0.769, -0.337,  0.000,  0.000 },
                {  0.433,  0.577,  0.000,  0.000, -0.048,  0.096,  0.096, -0.192, -0.240, -0.481,  0.481,  0.962, -0.144, -0.192,  0.000,  0.000 },
                {  0.192,  0.096,  0.385,  0.048, -0.577,  0.481,  0.577, -0.433, -0.769, -0.096,  0.769, -0.144,  1.154, -0.481,  0.962, -0.240 },
                { -0.096, -0.769, -0.048, -0.096,  0.481, -0.577, -0.529,  0.577,  0.096,  0.192, -0.337, -0.192, -0.481,  1.154, -0.240,  0.481 },
                {  0.385,  0.048,  0.769,  0.096, -0.577,  0.529,  0.000,  0.000, -0.769, -0.337,  0.000,  0.000,  0.962, -0.240,  1.923, -0.481 },
                { -0.048, -0.096, -0.096, -0.192,  0.433, -0.577,  0.000,  0.000, -0.144,  0.192,  0.000,  0.000, -0.240,  0.481, -0.481,  0.962 }
            });

            Matrix expectedStiffnessStdFirst = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                {  1.154,  0.481, -0.769,  0.096, -0.577, -0.481,  0.192, -0.096,  0.962,  0.240,  0.769,  0.144,  0.577,  0.433,  0.385, -0.048 },
                {  0.481,  1.154, -0.096,  0.192, -0.481, -0.577,  0.096, -0.769,  0.240,  0.481,  0.337, -0.192,  0.529,  0.577,  0.048, -0.096 },
                { -0.769, -0.096,  1.154, -0.481,  0.192,  0.096, -0.577,  0.481, -0.769,  0.144, -0.962,  0.240, -0.385, -0.048, -0.577,  0.433 },
                {  0.096,  0.192, -0.481,  1.154, -0.096, -0.769,  0.481, -0.577,  0.337,  0.192,  0.240, -0.481,  0.048,  0.096,  0.529, -0.577 },
                { -0.577, -0.481,  0.192, -0.096,  1.154,  0.481, -0.769,  0.096, -0.577, -0.433, -0.385,  0.048, -0.962, -0.240, -0.769, -0.144 },
                { -0.481, -0.577,  0.096, -0.769,  0.481,  1.154, -0.096,  0.192, -0.529, -0.577, -0.048,  0.096, -0.240, -0.481, -0.337,  0.192 },
                {  0.192,  0.096, -0.577,  0.481, -0.769, -0.096,  1.154, -0.481,  0.385,  0.048,  0.577, -0.433,  0.769, -0.144,  0.962, -0.240 },
                { -0.096, -0.769,  0.481, -0.577,  0.096,  0.192, -0.481,  1.154, -0.048, -0.096, -0.529,  0.577, -0.337, -0.192, -0.240,  0.481 },
                {  0.962,  0.240, -0.769,  0.337, -0.577, -0.529,  0.385, -0.048,  1.923,  0.481,  0.000,  0.000,  0.000,  0.000,  0.769, -0.096 },
                {  0.240,  0.481,  0.144,  0.192, -0.433, -0.577,  0.048, -0.096,  0.481,  0.962,  0.000,  0.000,  0.000,  0.000,  0.096, -0.192 },
                {  0.769,  0.337, -0.962,  0.240, -0.385, -0.048,  0.577, -0.529,  0.000,  0.000,  1.923, -0.481,  0.769,  0.096,  0.000,  0.000 },
                {  0.144, -0.192,  0.240, -0.481,  0.048,  0.096, -0.433,  0.577,  0.000,  0.000, -0.481,  0.962, -0.096, -0.192,  0.000,  0.000 },
                {  0.577,  0.529, -0.385,  0.048, -0.962, -0.240,  0.769, -0.337,  0.000,  0.000,  0.769, -0.096,  1.923,  0.481,  0.000,  0.000 },
                {  0.433,  0.577, -0.048,  0.096, -0.240, -0.481, -0.144, -0.192,  0.000,  0.000,  0.096, -0.192,  0.481,  0.962,  0.000,  0.000 },
                {  0.385,  0.048, -0.577,  0.529, -0.769, -0.337,  0.962, -0.240,  0.769,  0.096,  0.000,  0.000,  0.000,  0.000,  1.923, -0.481 },
                { -0.048, -0.096,  0.433, -0.577, -0.144,  0.192, -0.240,  0.481, -0.096, -0.192,  0.000,  0.000,  0.000,  0.000, -0.481,  0.962 }
            });

            XContinuumElement2D element = CreateCrackBodyElement();
            IMatrix stiffnessNodeMajor = element.JoinStiffnessesNodeMajor();
            IMatrix stiffnessStdFirst = element.JoinStiffnessesStandardFirst();

            Assert.True(expectedStiffnessNodeMajor.Equals(
                stiffnessNodeMajor.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedStiffnessStdFirst.Equals(
                stiffnessStdFirst.DoToAllEntries(round), equalityTolerance));
        }

        [Fact]
        public static void TestStiffnessSubmatricesOfCrackBodyElement()
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

            XContinuumElement2D element = CreateCrackBodyElement();
            Matrix Kss = element.BuildStandardStiffnessMatrix();
            (Matrix Kee1, Matrix Kes) = element.BuildEnrichedStiffnessMatricesLower();
            (Matrix Kee2, Matrix Kse) = element.BuildEnrichedStiffnessMatricesUpper();

            Assert.True(expectedKss.Equals(Kss.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedKee.Equals(Kee1.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedKee.Equals(Kee2.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedKes.Equals(Kes.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedKse.Equals(Kse.DoToAllEntries(round), equalityTolerance));
        }

        [Fact]
        public static void TestStiffnessSubmatricesOfMaterialInterfaceElement()
        {
            double equalityTolerance = 1E-13;
            Matrix expectedKss = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                { +0.913, +0.421, -0.577, +0.012, -0.433, -0.349, +0.096, -0.084 },
                { +0.421, +1.034, -0.132, +0.144, -0.373, -0.433, +0.084, -0.745 },
                { -0.577, -0.132, +0.817, -0.300, +0.192, +0.060, -0.433, +0.373 },
                { +0.012, +0.144, -0.300, +0.697, -0.060, -0.409, +0.349, -0.433 },
                { -0.433, -0.373, +0.192, -0.060, +0.817, +0.300, -0.577, +0.132 },
                { -0.349, -0.433, +0.060, -0.409, +0.300, +0.697, -0.012, +0.144 },
                { +0.096, +0.084, -0.433, +0.349, -0.577, -0.012, +0.913, -0.421 },
                { -0.084, -0.745, +0.373, -0.433, +0.132, +0.144, -0.421, +1.034 }
            });

            Matrix expectedKse = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                { +1.242, +0.080, +1.643, +0.160, +1.723, +2.083, +2.123, +2.484 },
                { +1.122, -2.865, +1.042, -1.462, +2.324, +2.424, +2.724, +3.826 },
                { -2.845, +0.881, -2.925, +0.801, -0.441, -1.122, -0.521, -1.522 },
                { +0.321, -2.744, +0.401, -3.025, -0.881, +2.063, -1.282, +1.783 },
                { -0.521, +1.522, -0.441, +1.122, -2.925, -0.801, -2.845, -0.881 },
                { +1.282, +1.783, +0.881, +2.063, -0.401, -3.025, -0.321, -2.744 },
                { +2.123, -2.484, +1.723, -2.083, +1.643, -0.160, +1.242, -0.080 },
                { -2.724, +3.826, -2.324, +2.424, -1.042, -1.462, -1.122, -2.865 }
            });

            Matrix expectedKes = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                { +1.242, +1.122, -2.845, +0.321, -0.521, +1.282, +2.123, -2.724 },
                { +0.080, -2.865, +0.881, -2.744, +1.522, +1.783, -2.484, +3.826 },
                { +1.643, +1.042, -2.925, +0.401, -0.441, +0.881, +1.723, -2.324 },
                { +0.160, -1.462, +0.801, -3.025, +1.122, +2.063, -2.083, +2.424 },
                { +1.723, +2.324, -0.441, -0.881, -2.925, -0.401, +1.643, -1.042 },
                { +2.083, +2.424, -1.122, +2.063, -0.801, -3.025, -0.160, -1.462 },
                { +2.123, +2.724, -0.521, -1.282, -2.845, -0.321, +1.242, -1.122 },
                { +2.484, +3.826, -1.522, +1.783, -0.881, -2.744, -0.080, -2.865 }
            });

            Matrix expectedKee = 1E6 * Matrix.CreateFromArray(new double[,]
            {
                {  95.753,  -6.010,  49.279, -3.606,  18.029, -10.817, 38.862,  1.202 },
                {  -6.010,  46.675,  -8.413, 28.245, -13.221,  -9.014, -1.202, -8.213 },
                {  49.279,  -8.413,  94.151, -6.010,  40.465,   1.202, 18.029, 13.221 },
                {  -3.606,  28.245,  -6.010, 41.066,  -1.202,  -2.604, 10.817, -9.014 },
                {  18.029, -13.221,  40.465, -1.202,  94.151,   6.010, 49.279,  8.413 },
                { -10.817,  -9.014,  1.2020, -2.604,   6.010,  41.066,  3.606, 28.245 },
                {  38.862,  -1.202,  18.029, 10.817,  49.279,   3.606, 95.753,  6.010 },
                {  1.2020,  -8.213,  13.221, -9.014,   8.413,  28.245,  6.010, 46.675 }
            });

            XContinuumElement2D element = CreateMaterialInterfaceElement();
            Matrix Kss = element.BuildStandardStiffnessMatrix();
            (Matrix Kee1, Matrix Kes) = element.BuildEnrichedStiffnessMatricesLower();
            (Matrix Kee2, Matrix Kse) = element.BuildEnrichedStiffnessMatricesUpper();

            
            Assert.True(expectedKss.Equals(Kss.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedKee.Equals(Kee1.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedKee.Equals(Kee2.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedKes.Equals(Kes.DoToAllEntries(round), equalityTolerance));
            Assert.True(expectedKse.Equals(Kse.DoToAllEntries(round), equalityTolerance));
        }

        private static XContinuumElement2D CreateCrackBodyElement()
        {
            XNode[] nodes = new XNode[]
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

            double E = 2e6;
            double v = 0.3;
            var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);

            var integrationStrategy = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.GetQuadratureWithOrder(2, 2)),
                new SimpleIntegration2D());
            //var integrationStrategy = new IntegrationForCrackPropagation2D(GaussLegendre2D.GetQuadratureWithOrder(2, 2),
            //  new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.GetQuadratureWithOrder(2, 2)));
            var factory = new XContinuumElement2DFactory(integrationStrategy, integrationStrategy, material);
            var bodyElement = factory.CreateElement(0, CellType.Quad4, new XNode[] { nodes[0], nodes[1], nodes[2], nodes[3] });
            var blendingElement = factory.CreateElement(1, CellType.Quad4, new XNode[] { nodes[7], nodes[6], nodes[1], nodes[0] });
            var tipElement = factory.CreateElement(2, CellType.Quad4, new XNode[] { nodes[4], nodes[5], nodes[6], nodes[7] });
            var boundary = new Rectangular2DBoundary(20.0, 40.0, -40.0, 20.0);
            var mesh = new SimpleMesh2D<XNode, XContinuumElement2D>(nodes,
                new XContinuumElement2D[] { bodyElement, blendingElement, tipElement }, boundary);

            var crack = new BasicExplicitCrack2D();
            crack.Mesh = mesh;
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);

            var point1 = new CartesianPoint(30.0, 20.0);
            var point2 = new CartesianPoint(30.0, -30.0);
            crack.InitializeGeometry(point1, point2);
            crack.UpdateEnrichments();

            return bodyElement;
        }

        private static XContinuumElement2D CreateMaterialInterfaceElement()
        {
            XNode[] nodes = new XNode[]
            {
                new XNode(0, 20.0, 0.0),
                new XNode(1, 40.0, 0.0),
                new XNode(2, 40.0, 20.0),
                new XNode(3, 20.0, 20.0),
            };

            var point1 = new CartesianPoint(30.0, 0.0);
            var point2 = new CartesianPoint(30.0, 20.0);
            PolyLine2D discontinuity = new PolyLine2D(point1, point2);
            MaterialInterface2D enrichmentItem = new MaterialInterface2D(discontinuity);

            double E1 = 2E6;
            double E2 = 0.5 * E1;
            double v = 0.3;
            var material = BiElasticMaterial2D.CreateMaterialForPlainStrain(E1, v, E2, v, enrichmentItem);

            var integrationStrategy = new RectangularSubgridIntegration2D<XContinuumElement2D>(
                2, GaussLegendre2D.GetQuadratureWithOrder(2, 2));

            //var integrationStrategy = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
            //    new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.Order2x2));
            var factory = new XContinuumElement2DFactory(integrationStrategy, integrationStrategy, material);
            var element = factory.CreateElement(0, CellType.Quad4, nodes);

            //discontinuity.ElementIntersections.Add(element, new CartesianPoint2D[] { point1, point2 });

            //OBSOLETE: Elements access their enrichments from nodes now.
            //enrichmentItem.EnrichElement(element);
            enrichmentItem.EnrichNode(nodes[0]);
            enrichmentItem.EnrichNode(nodes[1]);
            enrichmentItem.EnrichNode(nodes[2]);
            enrichmentItem.EnrichNode(nodes[3]);

            return element;
        }
    }
}
