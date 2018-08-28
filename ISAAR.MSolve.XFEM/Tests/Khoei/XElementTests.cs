using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Boundaries;
using ISAAR.MSolve.XFEM.Geometry.Mesh;
using ISAAR.MSolve.XFEM.Geometry.Descriptions;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.CrackGeometry.Explicit;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class XElementTests
    {
        private static XNode2D[] nodes = {
            new XNode2D(0, 20.0, 0.0),
            new XNode2D(1, 40.0, 0.0),
            new XNode2D(2, 40.0, 20.0),
            new XNode2D(3, 20.0, 20.0),

            new XNode2D(4, 20.0, -40.0),
            new XNode2D(5, 40.0, -40.0),
            new XNode2D(6, 40.0, -20.0),
            new XNode2D(7, 20.0, -20.0)
        };

        private static void IsoparametricQuad4WithCrackTest()
        {
            double E = 2e6;
            double v = 0.3;
            var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStrain(E, v);

            var integrationStrategy = new IntegrationForCrackPropagation2D(
                new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2),
                new XSimpleIntegration2D());
            //var integrationStrategy = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
            //new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.Order2x2));
            var bodyElement = new XContinuumElement2D(IsoparametricElementType2D.Quad4, 
                new XNode2D[] { nodes[0], nodes[1], nodes[2], nodes[3]}, material, integrationStrategy);
            var blendingElement = new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                new XNode2D[] { nodes[7], nodes[6], nodes[1], nodes[0] }, material, integrationStrategy);
            var tipElement = new XContinuumElement2D(IsoparametricElementType2D.Quad4,
                new XNode2D[] { nodes[4], nodes[5], nodes[6], nodes[7] }, material, integrationStrategy);
            var boundary = new RectangularBoundary(20.0, 40.0, -40.0, 20.0);
            var mesh = new SimpleMesh2D<XNode2D, XContinuumElement2D>(nodes, 
                new XContinuumElement2D[] { bodyElement, blendingElement, tipElement }, boundary);

            var crack = new BasicExplicitCrack2D();
            crack.Mesh = mesh;
            crack.CrackBodyEnrichment = new CrackBodyEnrichment2D(crack);
            crack.CrackTipEnrichments = new CrackTipEnrichments2D(crack, CrackTipPosition.Single);

            var point1 = new CartesianPoint2D(30.0, 20.0);
            var point2 = new CartesianPoint2D(30.0, -30.0);
            crack.InitializeGeometry(point1, point2);
            crack.UpdateEnrichments();

            Matrix Kss = bodyElement.BuildStandardStiffnessMatrix();
            (Matrix Kee, Matrix Kes) = bodyElement.BuildEnrichedStiffnessMatricesLower();

            var writer = new FullMatrixWriter();
            Console.WriteLine("Quad4 standard-standard stiffness matrix = ");
            writer.WriteToConsole(Kss.Scale(1e-6));
            Console.WriteLine("Quad4 enriched-standard stiffness matrix = ");
            writer.WriteToConsole(Kes.Scale(1e-6));
            Console.WriteLine("Quad4 enriched-enriched stiffness matrix = ");
            writer.WriteToConsole(Kee.Scale(1e-6));
        }

        private static void IsoparametricQuad4BimaterialTest(XNode2D[] nodes)
        {
            double E = 2e6;
            double v = 0.3;

            var point1 = new CartesianPoint2D(30.0, 20.0);
            var point2 = new CartesianPoint2D(30.0, 0.0);
            Polyline2DOLD discontinuity = new Polyline2DOLD(point1, point2);
            MaterialInterface2D enrichmentItem = new MaterialInterface2D(discontinuity);
            var material = BiElasticMaterial2D.CreateMaterialForPlainStrain(0.5 * E, v, E, v, enrichmentItem);

            var integrationStrategy = new IntegrationForCrackPropagation2D(new XSimpleIntegration2D(),
                    new RectangularSubgridIntegration2D<XContinuumElement2D>(8, GaussLegendre2D.Order2x2));
            //var integrationStrategy = new IntegrationForCrackPropagation2D(GaussLegendre2D.Order2x2,
            //    new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.Order2x2));
            var element = new XContinuumElement2D(IsoparametricElementType2D.Quad4, nodes, material, 
                integrationStrategy);

            discontinuity.ElementIntersections.Add(element, new CartesianPoint2D[] { point1, point2 });
            enrichmentItem.EnrichElement(element);
            enrichmentItem.EnrichNode(nodes[0]);
            enrichmentItem.EnrichNode(nodes[1]);
            enrichmentItem.EnrichNode(nodes[2]);
            enrichmentItem.EnrichNode(nodes[3]);


            Matrix Kss = element.BuildStandardStiffnessMatrix();
            (Matrix Kee, Matrix Kes) = element.BuildEnrichedStiffnessMatricesLower();

            var writer = new FullMatrixWriter();
            Console.WriteLine("Quad4 standard-standard stiffness matrix = ");
            writer.WriteToConsole(Kss.Scale(1e-6));
            Console.WriteLine("Quad4 enriched-standard stiffness matrix = ");
            writer.WriteToConsole(Kes.Scale(1e-6));
            Console.WriteLine("Quad4 enriched-enriched stiffness matrix = ");
            writer.WriteToConsole(Kee.Scale(1e-6));
        }

        static void Main(string[] args)
        {
            IsoparametricQuad4WithCrackTest();
            //IsoparametricQuad4BimaterialTest(NodeSets.nodeSet7);
        }
    }
}
