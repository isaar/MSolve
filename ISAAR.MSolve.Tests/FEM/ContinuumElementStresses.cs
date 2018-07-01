using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Materials;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ContinuumElementStresses
    {
        private static double thickness = 1.0;

        private static readonly IReadOnlyList<Node2D> quad4NodeSet1 = new Node2D[]
        {
            new Node2D(0, -1.0, -1.0),
            new Node2D(1, +1.0, -1.0),
            new Node2D(2, +1.0, +1.0),
            new Node2D(3, -1.0, +1.0)
        };

        private static readonly IReadOnlyList<Node2D> quad4NodeSet2 = new Node2D[]
        {
            new Node2D(0, 0.2, 0.3),
            new Node2D(1, 2.2, 1.5),
            new Node2D(2, 3.0, 2.7),
            new Node2D(3, 0.7, 2.0)
        };

        private static readonly ElasticMaterial2D material = new ElasticMaterial2D(StressState2D.PlaneStress)
        {
            YoungModulus = 210000,
            PoissonRatio = 0.3
        };

        [Fact]
        public static void TestQuad4StrainsStresses()
        {
            var factory = new ContinuumElement2DFactory(thickness, material, null);
            ContinuumElement2D quad4 = factory.CreateQuad4(quad4NodeSet1);

            // Abaqus results
            double[] displacements =
            {
                  0.0,          0.0,          // Node 1
                  0.0,          0.0,          // Node 2
                  1.0,       -499.614E-03,    // Node 4  
                986.100E-03,    1.0           // Node 3
            };

            // The Gauss points in Abaqus have the same order as in GaussLegendre2D.Order2x2: Xi major, Eta minor
            double[][] expectedStrainsAtGPs =
            {
                new double[] { 1.46867E-03,  341.547E-03,  336.066E-03 },  // Gauss point 1
                new double[] { 1.46867E-03,  -91.3541E-03, 340.078E-03 },  // Gauss point 2
                new double[] { 5.48114E-03,  341.547E-03,  -96.8352E-03 }, // Gauss point 3
                new double[] { 5.48114E-03,  -91.3541E-03, -92.8228E-03 }  // Gauss point 4
            };
            double[][] expectedStressesAtGPs =
            {
                new double[] { 23.9845E+03,   78.9202E+03,  27.1438E+03 },  // Gauss point 1
                new double[] { -5.98559E+03, -20.9800E+03,  27.4679E+03 },  // Gauss point 2
                new double[] { 24.9104E+03,   79.1980E+03,  -7.82131E+03 }, // Gauss point 3
                new double[] { -5.05964E+03, -20.7023E+03,  -7.49722E+03 }  // Gauss point 4
            };

            // However the order of nodes is different, since in Abaqus they are numbered by the preprocessor, while Quad4  
            // nodes here have a counter-clockwise order.
            double[][] expectedStrainsAtNodes =
            {
                new double[] { 58.2077E-12,  500.000E-03,  493.050E-03 },  // Node 1
                new double[] {  0.0,        -249.807E-03,  500.0E-03 },    // Node 2
                new double[] { 6.94981E-03, -249.807E-03, -249.807E-03 },  // Node 4
                new double[] { 6.94981E-03,  500.000E-03, -256.757E-03 }   // Node 3
            };
            double[][] expectedStressesAtNodes =
            {
                new double[] {  34.6154E+03, 115.385E+03,   39.8233E+03 },  // Node 1
                new double[] { -17.2943E+03, -57.6478E+03,  40.3846E+03 },  // Node 2
                new double[] { -15.6905E+03, -57.1666E+03, -20.1767E+03 },  // Node 4
                new double[] {  36.2192E+03, 115.866E+03,  -20.7380E+03 }   // Node 3
            };

            (IReadOnlyList<double[]> strainsAtGPs, IReadOnlyList<double[]> stressesAtGPs) = 
                quad4.UpdateStrainsStressesAtGaussPoints(displacements);
            IReadOnlyList<double[]> strainsAtNodes =
                quad4.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, quad4.Interpolation);
            IReadOnlyList<double[]> stressesAtNodes =
                quad4.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, quad4.Interpolation);

            Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtGPs, strainsAtGPs, 1e-4));
            Assert.True(Utilities.AreTensorsEqual(expectedStressesAtGPs, stressesAtGPs, 1e-4));
            Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtNodes, strainsAtNodes, 1e-4));
            Assert.True(Utilities.AreTensorsEqual(expectedStressesAtNodes, stressesAtNodes, 1e-4));
        }
    }
}
