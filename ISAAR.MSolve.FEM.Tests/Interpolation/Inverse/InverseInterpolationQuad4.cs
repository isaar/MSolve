using System;
using System.Collections.Generic;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Numerical.Commons;
using Xunit;

namespace ISAAR.MSolve.FEM.Tests.Interpolation.Inverse
{
    public static class InverseInterpolationQuad4
    {
        private static readonly ValueComparer comparer = new ValueComparer(1E-10);

        /// <summary>
        /// Random shape, not too distorted.
        /// </summary>
        private static readonly IReadOnlyList<Node_v2> nodeSet = new Node_v2[]
        {
            new Node_v2 { ID = 0, X = 0.7, Y = 2.0 },
            new Node_v2 { ID = 1, X = 0.2, Y = 0.3 },
            new Node_v2 { ID = 2, X = 2.0, Y = 0.9 },
            new Node_v2 { ID = 3, X = 3.0, Y = 2.7 }
        };
        
        private static bool Coincide(NaturalPoint2D point1, NaturalPoint2D point2)
            => comparer.AreEqual(point1.Xi, point2.Xi) && comparer.AreEqual(point1.Eta, point2.Eta);

        /// <summary>
        /// Reorders the nodes such that the 1st one becomes the 2nd, the 2nd one becomes the 3rd, etc.
        /// </summary>
        /// <param name="originalOrder"></param>
        private static Node_v2[] CycleNodes(IReadOnlyList<Node_v2> originalOrder)
        {
            var cycled = new Node_v2[originalOrder.Count];
            cycled[0] = originalOrder[originalOrder.Count - 1];
            for (int i = 0; i < originalOrder.Count - 1; ++i) cycled[i + 1] = originalOrder[i];
            return cycled;
        }

        /// <summary>
        /// Generates random points in the square: -1 &lt;= xi &lt; 1 , -1 &lt;= eta &lt; 1
        /// </summary>
        /// <returns></returns>
        private static NaturalPoint2D[] GenerateRandomPointsInSquare(int numRandomPoints)
        {
            var rand = new Random();
            var randomPoints = new NaturalPoint2D[numRandomPoints];
            for (int i = 0; i < numRandomPoints; ++i)
            {
                double xi = -1.0 + rand.NextDouble() * 2.0;
                double eta = -1.0 + rand.NextDouble() * 2.0;
                randomPoints[i] = new NaturalPoint2D(xi, eta);
            }
            return randomPoints;
        }

        [Fact]
        private static void TestInverseMapping()
        {
            var directMapping = InterpolationQuad4.UniqueInstance;
            int numRandomPoints = 10;
            NaturalPoint2D[] naturalPoints = GenerateRandomPointsInSquare(numRandomPoints);
            IReadOnlyList<Node_v2> elementNodes = nodeSet;

            for (int i = 0; i < 4; ++i)
            {
                IInverseInterpolation2D inverseMapping = directMapping.CreateInverseMappingFor(elementNodes);
                foreach (NaturalPoint2D originalPoint in naturalPoints)
                {
                    CartesianPoint2D cartesianPoint = directMapping.TransformNaturalToCartesian(elementNodes, originalPoint);
                    NaturalPoint2D remappedPoint = inverseMapping.TransformPointCartesianToNatural(cartesianPoint);
                    Assert.True(Coincide(originalPoint, remappedPoint));
                }

                elementNodes = CycleNodes(elementNodes); // The next iteration will use a different node order
            }
        }
    }
}
