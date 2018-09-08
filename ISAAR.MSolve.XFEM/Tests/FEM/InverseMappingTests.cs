using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Interpolation.InverseMappings;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests.FEM
{
    static class InverseMappingTests
    {
        /*
         * a1 = x0- x1 +x2 -x3
         * a2 = y0 -y1 +y2 -y3
         * b1 = -x0 +x1 +x2 -x3
         * b2 = -y0 +y1 +y2 -y3
         * c1 = -x0 -x1 +x2 +x3
         * c2 = -y0 -y1 +y2 +y3
         * ab = a1*b2 - a2*b1
         * ac = a1*c2 - a2*c1
         */

        /// <summary>
        /// a1=0, a2=0. Case 6.
        /// </summary>
        public static Node2D[] nodes_1_1 = {
            new Node2D(0, -2.0, -3.0),
            new Node2D(1, 2.0, -3.0),
            new Node2D(2, 2.0, 3.0),
            new Node2D(3, -2.0, 3.0)
        };

        /// <summary>
        /// a1=0, a2!=0, c1=0. Case 6
        /// </summary>
        public static Node2D[] nodes_1_2_1 = {
            new Node2D(0, 1.0, -3.0),
            new Node2D(1, 2.0, -3.0),
            new Node2D(2, 2.0, 5.0),
            new Node2D(3, 1.0, 3.0)
        };

        /// <summary>
        /// a1=0, a2!=0, c1!=0. Case 2
        /// </summary>
        public static Node2D[] nodes_1_2_2 = {
            new Node2D(0, -4.0, -3.0),
            new Node2D(1, -1.0, -3.0),
            new Node2D(2, 4.0, 3.0),
            new Node2D(3, 1.0, 5.0)
        };

        /// <summary>
        /// a1!=0, a2!=0, ab!=0, ac!=0. Case 1
        /// </summary>
        public static Node2D[] nodes_2_1_1_1 = {
            new Node2D(0, -4.0, -3.0),
            new Node2D(1, 0.0, -3.5),
            new Node2D(2, 4.0, 3.0),
            new Node2D(3, 1.0, 5.0)
        };

        /// <summary>
        /// a1!=0, a2!=0, ab!=0, ac=0. Case 5
        /// </summary>
        public static Node2D[] nodes_2_1_1_2 = {
            new Node2D(0, -4.0, -3.0),
            new Node2D(1, 0.0, -3.4),
            new Node2D(2, 4.0, 3.0),
            new Node2D(3, 1.0, 5.0)
        };

        /// <summary>
        /// a1!=0, a2!=0, ab=0, (if ab=0 then ac!=0). Case 4
        /// </summary>
        public static Node2D[] nodes_2_1_2 = {
            new Node2D(0, -4.0, -2.0),
            new Node2D(1, 0.0, -4.66666666666667),
            new Node2D(2, 4.0, 3.0),
            new Node2D(3, 1.0, 5.0)
        };

        /// <summary>
        /// a1!=0, a2=0, b2=0. Case 6
        /// </summary>
        public static Node2D[] nodes_2_2_1 = {
            new Node2D(0, -4.0, -3.0),
            new Node2D(1, 0.0, -3),
            new Node2D(2, 4.0, 3.0),
            new Node2D(3, 1.0, 3.0)
        };

        /// <summary>
        /// a1!=0, a2=0, b2!=0. Case 3
        /// </summary>
        public static Node2D[] nodes_2_2_2 = {
            new Node2D(0, -4.0, -3.0),
            new Node2D(1, 0.0, -4.0),
            new Node2D(2, 4.0, 3.0),
            new Node2D(3, 1.0, 4.0)
        };

        private static void PrintCoefficients(Node2D[] nodes)
        {
            double x1 = nodes[0].X, x2 = nodes[1].X, x3 = nodes[2].X, x4 = nodes[3].X;
            double y1 = nodes[0].Y, y2 = nodes[1].Y, y3 = nodes[2].Y, y4 = nodes[3].Y;

            // Calculate coefficients. TODO: cover all node ordering cases
            double sum1 = x1 + x2 + x3 + x4;
            double a1 = x1 - x2 + x3 - x4;
            double b1 = -x1 + x2 + x3 - x4;
            double c1 = -x1 - x2 + x3 + x4;

            double sum2 = y1 + y2 + y3 + y4;
            double a2 = y1 - y2 + y3 - y4;
            double b2 = -y1 + y2 + y3 - y4;
            double c2 = -y1 - y2 + y3 + y4;

            // Calculate determinants
            double ab = a1 * b2 - a2 * b1;
            double ac = a1 * c2 - a2 * c1;
            double bc = b1 * c2 - b2 * c1;

            Console.WriteLine("sum1 = " + sum1);
            Console.WriteLine("a1 = " + a1);
            Console.WriteLine("b1 = " + b1);
            Console.WriteLine("c1 = " + c1);
            Console.WriteLine("sum2 = " + sum2);
            Console.WriteLine("a2 = " + a2);
            Console.WriteLine("b2 = " + b2);
            Console.WriteLine("c2 = " + c2);
            Console.WriteLine("ab = " + ab);
            Console.WriteLine("ac = " + ac);
            Console.WriteLine("bc = " + bc);
            Console.WriteLine();
        }

        private static void Run(Node2D[] nodes)
        {
            PrintCoefficients(nodes);
            bool allCorrect = true;

            double E = 1;
            double v = 0.25;
            double t = 1.0;
            var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, t);
            var element = new ContinuumElement2D(IsoparametricElementType2D.Quad4, nodes, 
                new SimpleIntegration2D(), material);

            IsoparametricInterpolation2D interpolation = element.Interpolation;
            IInverseMapping2D inverseMapping = interpolation.CreateInverseMappingFor(nodes);
            

            const int pointsCount = 20;
            double naturalStep = 2.0 / pointsCount;

            for (int i = 0; i <= pointsCount; ++i)
            {
                double xi = -1.0 + i * naturalStep;
                if ((xi < -1.0) || (xi > 1.0)) continue;
                for (int j = 0; j <= pointsCount; ++j)
                {
                    //if (!((i == 0) && (j == 6))) continue; // For debugging purposes only
                    double eta = -1.0 + j * naturalStep;
                    if ((eta < -1.0) || (eta > 1.0)) continue;

                    NaturalPoint2D naturalOriginal = new NaturalPoint2D(xi, eta);
                    ICartesianPoint2D cartesian =
                        interpolation.EvaluateAt(nodes, naturalOriginal).TransformPointNaturalToGlobalCartesian(naturalOriginal);
                    INaturalPoint2D naturalReconstructed = inverseMapping.TransformCartesianToNatural(cartesian);

                    Console.WriteLine("Point (" + i + "," + j + "):");
                    Console.WriteLine("Original natural point: " + naturalOriginal);
                    Console.WriteLine("Cartesian point: " + cartesian);
                    Console.WriteLine("Reconstructed natural point: " + naturalReconstructed);


                    // There should be a method in the Point class for this.
                    double tolerance = 1e-8;
                    double dXi = Math.Abs(naturalOriginal.Xi - naturalReconstructed.Xi);
                    double dEta = Math.Abs(naturalOriginal.Eta - naturalReconstructed.Eta);
                    if (dXi <= tolerance && dEta <= tolerance)
                    {
                        Console.WriteLine("Correct mapping!\n");
                    }
                    else
                    {
                        allCorrect = false;
                        Console.WriteLine("Wrong mapping!\n");
                    }
                }
            }

            if (allCorrect) Console.WriteLine("----------------- All correct! ------------------\n");
            else Console.WriteLine("----------------- There is an error! ------------------\n");
        
           
        }

        public static void Main()
        {
            //Run(nodes_1_1);
            //Run(nodes_1_2_1);
            //Run(nodes_1_2_2);
            //Run(nodes_2_1_1_1);
            //Run(nodes_2_1_1_2);
            //Run(nodes_2_2_1);
            Run(nodes_2_2_2);
        }
    }
}
