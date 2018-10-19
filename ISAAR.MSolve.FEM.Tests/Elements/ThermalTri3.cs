using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.FEM.Tests.Elements
{
    public static class ThermalTri3
    {
        private static double thickness = 1.0;
        private static double thermalConductivity = 1.0;
        private static double density = 12.0;
        private static double specialHeatCoeff = 2.0;

        /// <summary>
        /// Random shape, not too distorted.
        /// </summary>
        private static readonly IReadOnlyList<Node2D> nodeSet0 = new Node2D[]
        {
            new Node2D(0, 0.0,  0.0),
            new Node2D(1, 1.0,  0.0),
            new Node2D(2, 0.0,  1.0)
        };

        [Fact]
        private static void TestCapacity()
        {
            var factory = new ThermalElement2DFactory(thickness, new ThermalMaterial(density, specialHeatCoeff, thermalConductivity));
            ThermalElement2D element = factory.CreateElement(CellType2D.Tri3, nodeSet0);
            IMatrix2D M = element.BuildCapacityMatrix();

            var expectedM = new Matrix2D(new double[,]
            {
                {2, 1, 1 },
                {1, 2, 1 },
                {1, 1, 2 }

            });

            Assert.True(Utilities.AreMatricesEqual(M, expectedM, 1e-10));
        }

        [Fact]
        private static void TestConductivity()
        {
            var factory = new ThermalElement2DFactory(thickness, new ThermalMaterial(density, specialHeatCoeff, thermalConductivity));
            ThermalElement2D element = factory.CreateElement(CellType2D.Tri3, nodeSet0);
            IMatrix2D K = element.BuildConductivityMatrix();

            double[,] expectedK = new double[,]
            {
                { 1.0, -0.5, -0.5 },
                { -0.5, 0.5,  0.0 },
                { -0.5, 0.0,  0.5 }
            };
            Assert.True(Utilities.AreMatricesEqual(K, new Matrix2D(expectedK), 1e-10));
        }
    }
}
