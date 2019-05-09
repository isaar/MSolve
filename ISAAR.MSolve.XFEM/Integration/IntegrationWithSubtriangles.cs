using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Integration
{
    // TODO: Needs option to refine the mesh for the J-integral and tip blending elements.
    public class IntegrationWithSubtriangles: IIntegrationStrategy2D<XContinuumElement2D>
    {
        private readonly TriangleQuadratureSymmetricGaussian triangleIntegrationRule;
        private readonly ITriangulator2D<NaturalPoint> triangulator;
        private readonly ISingleCrack crack;
        private readonly double triangleOverElementArea;

        public IntegrationWithSubtriangles(TriangleQuadratureSymmetricGaussian triangleIntegrationRule, ISingleCrack crack, 
            ITriangulator2D<NaturalPoint> triangulator, double triangleOverElementArea = double.PositiveInfinity)
        {
            this.triangleIntegrationRule = triangleIntegrationRule;
            this.crack = crack;
            this.triangulator = triangulator;
            this.triangleOverElementArea = triangleOverElementArea;
        }

        public IReadOnlyList<GaussPoint> GenerateIntegrationPoints(XContinuumElement2D element)
        {
            SortedSet<CartesianPoint> cartesianDelaunyPoints = crack.FindTriangleVertices(element);
            IReadOnlyList<NaturalPoint> naturalDelaunyPoints = 
                FindNaturalPointsForTriangulation(element, cartesianDelaunyPoints);
            IReadOnlyList<Triangle2D<NaturalPoint>> subtriangles;
            if (double.IsPositiveInfinity(triangleOverElementArea))
            {
                subtriangles = triangulator.CreateMesh(naturalDelaunyPoints);
            }
            else
            {
                double elementArea = (ConvexPolygon2D.CreateUnsafe(element.Nodes)).ComputeArea();
                subtriangles = triangulator.CreateMesh(naturalDelaunyPoints, triangleOverElementArea * elementArea);
            }

            var integrationPoints = new List<GaussPoint>();
            foreach (Triangle2D<NaturalPoint> triangle in subtriangles)
            {
                integrationPoints.AddRange(GenerateIntegrationPointsOfTriangle(triangle));
            }
            return integrationPoints;
        }

        // These triangles are output by the delauny triangulation and the order of their nodes might be 
        // counter-clockwise or clockwise. In the second case the jacobian will be negative, 
        // but it doesn't matter otherwise. 
        private IReadOnlyList<GaussPoint> GenerateIntegrationPointsOfTriangle(Triangle2D<NaturalPoint> triangle)
        {
            // Coordinates of the triangle's nodes in the natural system of the element
            double xi1 = triangle.Vertices[0].Xi;
            double eta1 = triangle.Vertices[0].Eta;
            double xi2 = triangle.Vertices[1].Xi;
            double eta2 = triangle.Vertices[1].Eta;
            double xi3 = triangle.Vertices[2].Xi;
            double eta3 = triangle.Vertices[2].Eta;

            // Determinant of the Jacobian of the linear mapping from the natural system of the triangle  to the 
            // natural system of the element. If the triangle's nodes are in clockwise order, the determinant will be 
            // negative. It doesn't matter since its absolute value is used for integration with change of variables.
            double jacobian = Math.Abs(xi1 * (eta2 - eta3) + xi2 * (eta3 - eta1) + xi3 * (eta1 - eta2));

            var triangleGaussPoints = triangleIntegrationRule.IntegrationPoints;
            var elementGaussPoints = new GaussPoint[triangleGaussPoints.Count];
            for (int i = 0; i < triangleGaussPoints.Count; ++i)
            {
                GaussPoint triangleGP = triangleGaussPoints[i];

                // Linear shape functions evaluated at the Gauss point's coordinates in the triangle's natural system.
                double N1 = 1.0 - triangleGP.Xi - triangleGP.Eta;
                double N2 = triangleGP.Xi;
                double N3 = triangleGP.Eta;

                // Coordinates of the same gauss point in the element's natural system
                double elementXi = N1 * xi1 + N2 * xi2 + N3 * xi3;
                double elementEta = N1 * eta1 + N2 * eta2 + N3 * eta3;

                // The integral would need to be multiplicated with |detJ|. 
                // It is simpler for the caller to have it already included in the weight.
                double elementWeight = triangleGP.Weight * jacobian;

                elementGaussPoints[i] = new GaussPoint(elementXi, elementEta, elementWeight);
            }

            return elementGaussPoints;
        }

        private static IReadOnlyList<NaturalPoint> FindNaturalPointsForTriangulation(XContinuumElement2D element,
            IEnumerable<CartesianPoint> cartesianDelaunyPoints)
        {
            IInverseInterpolation2D inverseMapping = element.Interpolation.CreateInverseMappingFor(element.Nodes);
            var naturalDelaunyPoints = new List<NaturalPoint>();
            foreach (var cartesianPoint in cartesianDelaunyPoints)
            {
                naturalDelaunyPoints.Add(inverseMapping.TransformPointCartesianToNatural(cartesianPoint));
            }
            return naturalDelaunyPoints;
        }
    }
}
