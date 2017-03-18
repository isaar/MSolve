using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Geometry.Descriptions;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests
{
    class XElementTests
    {
        private static void IsoparametricQuad4WithCrackTest(XNode2D[] nodes)
        {
            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStrain.Create(E, v, t);

            var point1 = new CartesianPoint2D(30.0, 20.0);
            var point2 = new CartesianPoint2D(30.0, 0.0);
            Polyline2D discontinuity = new Polyline2D(point1, point2);
            IEnrichmentItem2D enrichmentItem = new CrackBody2D(discontinuity);

            var integrationRule = new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.Order2x2);
            var element = new XContinuumElement2D(IsoparametricElementType2D.Quad4, nodes,
                new HomogeneousIntegration2D(integrationRule, material));

            discontinuity.ElementIntersections.Add(element, new CartesianPoint2D[] { point1, point2 });
            enrichmentItem.AffectElement(element);
            enrichmentItem.EnrichNodes();

            SymmetricMatrix2D<double> stiffnessStd = element.BuildStandardStiffnessMatrix();
            Matrix2D<double> stiffnessStdEnr;
            SymmetricMatrix2D<double> stiffnessEnr;
            element.BuildEnrichedStiffnessMatrices(out stiffnessStdEnr, out stiffnessEnr);

            stiffnessStd.Scale(1e-6);
            stiffnessStdEnr.Scale(1e-6);
            stiffnessEnr.Scale(1e-6);

            Console.WriteLine("Quad4 standard-standard stiffness matrix = ");
            Console.WriteLine(stiffnessStd);
            Console.WriteLine("Quad4 standard-enriched stiffness matrix = ");
            Console.WriteLine(stiffnessStdEnr);
            Console.WriteLine("Quad4 enriched-enriched stiffness matrix = ");
            Console.WriteLine(stiffnessEnr);
        }

        //private static void IsoparametricQuad4WithTipTest(XNode2D[] nodes)
        //{
        //    double E = 2e6;
        //    double v = 0.3;
        //    double t = 1.0;
        //    var material = ElasticMaterial2DPlainStrain.Create(E, v, t);

        //    ICartesianPoint2D tip = new CartesianPoint2D(30.0, 50.0);
        //    double angle = Math.PI / 2.0;
        //    CrackTip2D enrichmentItem = new CrackTip2D(tip, angle);

        //    var integrationRule = new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.Order2x2);
        //    var element = new XContinuumElement2D(IsoparametricElementType2D.Quad4, nodes,
        //        new HomogeneousIntegration2D(integrationRule, material));

        //    enrichmentItem.AffectElement(element);
        //    enrichmentItem.EnrichNodes();

        //    SymmetricMatrix2D<double> stiffnessStd = element.BuildStandardStiffnessMatrix();
        //    Matrix2D<double> stiffnessStdEnr;
        //    SymmetricMatrix2D<double> stiffnessEnr;
        //    element.BuildEnrichedStiffnessMatrices(out stiffnessStdEnr, out stiffnessEnr);

        //    stiffnessStd.Scale(1e-6);
        //    stiffnessStdEnr.Scale(1e-6);
        //    stiffnessEnr.Scale(1e-6);

        //    Console.WriteLine("Quad4 standard-standard stiffness matrix = ");
        //    Console.WriteLine(stiffnessStd);
        //    Console.WriteLine("Quad4 standard-enriched stiffness matrix = ");
        //    Console.WriteLine(stiffnessStdEnr);
        //    Console.WriteLine("Quad4 enriched-enriched stiffness matrix = ");
        //    Console.WriteLine(stiffnessEnr);
        //    Console.WriteLine();
        //}

        private static void IsoparametricQuad4BimaterialTest(XNode2D[] nodes)
        {
            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var materialLeft = ElasticMaterial2DPlainStrain.Create(E, v, t);
            var materialRight = ElasticMaterial2DPlainStrain.Create(E / 2.0, v, t);

            var point1 = new CartesianPoint2D(30.0, 20.0);
            var point2 = new CartesianPoint2D(30.0, 0.0);
            Polyline2D discontinuity = new Polyline2D(point1, point2);
            MaterialInterface2D enrichmentItem = new MaterialInterface2D(discontinuity, materialLeft, materialRight);

            var integrationRule = new RectangularSubgridIntegration2D<XContinuumElement2D>(2, GaussLegendre2D.Order2x2);
            var element = new XContinuumElement2D(
               IsoparametricElementType2D.Quad4, nodes, new BimaterialIntegration2D(integrationRule, enrichmentItem));

            discontinuity.ElementIntersections.Add(element, new CartesianPoint2D[] { point1, point2 });
            enrichmentItem.AffectElement(element);
            enrichmentItem.EnrichNodes();

            SymmetricMatrix2D<double> stiffnessStd = element.BuildStandardStiffnessMatrix();
            Matrix2D<double> stiffnessStdEnr;
            SymmetricMatrix2D<double> stiffnessEnr;
            element.BuildEnrichedStiffnessMatrices(out stiffnessStdEnr, out stiffnessEnr);

            stiffnessStd.Scale(1e-6);
            stiffnessStdEnr.Scale(1e-6);
            stiffnessEnr.Scale(1e-6);

            Console.WriteLine("Quad4 standard-standard stiffness matrix = ");
            Console.WriteLine(stiffnessStd);
            Console.WriteLine("Quad4 standard-enriched stiffness matrix = ");
            Console.WriteLine(stiffnessStdEnr);
            Console.WriteLine("Quad4 enriched-enriched stiffness matrix = ");
            Console.WriteLine(stiffnessEnr);
        }

        static void Main(string[] args)
        {
            IsoparametricQuad4WithCrackTest(NodeSets.nodeSet7);
            //IsoparametricQuad4WithTipTest(NodeSets.nodeSet8);
            //IsoparametricQuad4BimaterialTest(NodeSets.nodeSet7);
        }
    }
}
