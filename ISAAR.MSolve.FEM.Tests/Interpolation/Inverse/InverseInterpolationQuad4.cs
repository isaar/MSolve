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
        private static readonly IReadOnlyList<Node2D> nodeSet = new Node2D[]
        {
            new Node2D(0, 0.7, 2.0),
            new Node2D(1, 0.2, 0.3),
            new Node2D(2, 2.0, 0.9),
            new Node2D(3, 3.0, 2.7)
        };
        
        private static bool Coincide(NaturalPoint2D point1, NaturalPoint2D point2)
            => comparer.AreEqual(point1.Xi, point2.Xi) && comparer.AreEqual(point1.Eta, point2.Eta);

        /// <summary>
        /// Reorders the nodes such that the 1st one becomes the 2nd, the 2nd one becomes the 3rd, etc.
        /// </summary>
        /// <param name="originalOrder"></param>
        private static Node2D[] CycleNodes(IReadOnlyList<Node2D> originalOrder)
        {
            var cycled = new Node2D[originalOrder.Count];
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
            IReadOnlyList<Node2D> elementNodes = nodeSet;

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
