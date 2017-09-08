using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Geometry.Triangulation;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Interpolation.InverseMappings;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    class IntegrationWithSubtriangles: IIntegrationStrategy2D<XContinuumElement2D>
    {
        private readonly GaussQuadratureForTriangle triangleIntegrationRule;
        private readonly ITriangulator2D mesher;
        private readonly IExteriorCrack crack;

        public IntegrationWithSubtriangles(GaussQuadratureForTriangle triangleIntegrationRule, ITriangulator2D mesher, 
            IExteriorCrack crack)
        {
            this.triangleIntegrationRule = triangleIntegrationRule;
            this.mesher = mesher;
        }
        
        public IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints(XContinuumElement2D element)
        {
            IEnumerable<ICartesianPoint2D> cartesianDelaunyPoints = FindCartesianPointsForTriangulation(element);
            IReadOnlyList<INaturalPoint2D> naturalDelaunyPoints = 
                FindNaturalPointsForTriangulation(element, cartesianDelaunyPoints);
            IReadOnlyList<Triangle2D> subtriangles = mesher.CreateMesh(naturalDelaunyPoints);

            var integrationPoints = new List<GaussPoint2D>();
            foreach (Triangle2D triangle in subtriangles)
            {
                integrationPoints.AddRange(GenerateIntegrationPointsOfTriangle(triangle));
            }
            return integrationPoints;
        }

        // These triangles are output by the delauny triangulation and the order of their nodes might be 
        // counter-clockwise or clockwise. In the second case the jacobian will be negative, 
        // but it doesn't matter otherwise. 
        private IReadOnlyList<GaussPoint2D> GenerateIntegrationPointsOfTriangle(Triangle2D triangle)
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
            var elementGaussPoints = new GaussPoint2D[triangleGaussPoints.Count];
            for (int i = 0; i < triangleGaussPoints.Count; ++i)
            {
                GaussPoint2D triangleGP = triangleGaussPoints[i];

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

                elementGaussPoints[i] = new GaussPoint2D(elementXi, elementEta, elementWeight);
            }

            return elementGaussPoints;
        }

        private static IReadOnlyList<INaturalPoint2D> FindNaturalPointsForTriangulation(XContinuumElement2D element,
            IEnumerable<ICartesianPoint2D> cartesianDelaunyPoints)
        {
            IInverseMapping2D inverseMapping = element.Interpolation.CreateInverseMappingFor(element.Nodes);
            var naturalDelaunyPoints = new List<INaturalPoint2D>();
            foreach (var cartesianPoint in cartesianDelaunyPoints)
            {
                naturalDelaunyPoints.Add(inverseMapping.TransformCartesianToNatural(cartesianPoint));
            }
            return naturalDelaunyPoints;
        }

        // TODO: perhaps the triangulation should be done at the global system
        private static IEnumerable<ICartesianPoint2D> FindCartesianPointsForTriangulation(XContinuumElement2D element)
        {
            if (element.EnrichmentItems.Count == 0)
            {
                throw new Exception("Should not be used for standard elements");
            }
            else if (element.EnrichmentItems.Count > 1)
            {
                throw new NotImplementedException("I must also find the intersection points of the 2 enrichment items' curves");
            }
            
            var pointsOfInterest = new List<ICartesianPoint2D>(element.Nodes); // TODO: better to use a set
            foreach (IEnrichmentItem2D enrichment in element.EnrichmentItems)
            {
                pointsOfInterest.AddRange(enrichment.IntersectionPointsForIntegration(element));
            }
            return pointsOfInterest.Distinct(); // TODO: Provide a comparer that checks (x,y) proximity
        }
    }
}
