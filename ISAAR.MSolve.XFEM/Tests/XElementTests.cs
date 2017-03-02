using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Utilities;

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

            ICurve2D discontinuity = new Line2D(new CartesianPoint2D(30.0, 0.0), new CartesianPoint2D(30.0, 20.0));
            IEnrichmentItem2D enrichmentItem = new CrackBody2D(discontinuity);

            var element = XElement2D.CreateHomogeneous(new IsoparametricQuad4(nodes), material);

            enrichmentItem.AffectElement(element);
            enrichmentItem.EnrichNodes();

            SymmetricMatrix2D<double> stiffnessStd = element.BuildStdStiffnessMatrix();
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

        private static void IsoparametricQuad4WithTipTest(XNode2D[] nodes)
        {
            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStrain.Create(E, v, t);

            ICartesianPoint2D tip = new CartesianPoint2D(30.0, 50.0);
            double angle = Math.PI / 2.0;
            CrackTip2D enrichmentItem = new CrackTip2D(tip, angle);

            var element = XElement2D.CreateHomogeneous(new IsoparametricQuad4(nodes), material);

            enrichmentItem.AffectElement(element);
            enrichmentItem.EnrichNodes();

            SymmetricMatrix2D<double> stiffnessStd = element.BuildStdStiffnessMatrix();
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
            Console.WriteLine();

            //// Check enrichments at Nodes
            //for (int n = 0; n < 4; ++n)
            //{
            //    PolarPoint2D polar = enrichmentItem.TipSystem.TransformCartesianGlobalToPolarLocal(nodes[n]);
            //    Console.Write("Node " + n + ": ");
            //    Console.Write("(x,y) = (" + nodes[n].X + " , " + nodes[n].Y + ") - ");
            //    Console.WriteLine("(r, theta) = (" + polar.R + " , " + polar.Theta + ")");
            //    for (int enr = 0; enr < 4; ++enr)
            //    {
            //        IEnrichmentFunction2D func = enrichmentItem.EnrichmentFunctions[enr];
            //        EvaluatedFunction2D H = func.EvaluateAllAt(nodes[n]);
            //        Console.Write("Enrichment H" + enr + ": ");
            //        Console.WriteLine("H = " + H.Value + " , dH/dx = " + H.CartesianDerivatives.Item1 
            //            + " , dH/dy = " + H.CartesianDerivatives.Item2);
            //    }
            //    Console.WriteLine();
            //}
            //Console.WriteLine();

            // Check enrichments at GP
            //for (int i = 0; i < GaussPointsForTipTest.gaussPoints.Count; ++i)
            //{
            //    GaussPoint2D gp = GaussPointsForTipTest.gaussPoints[i];
            //    Console.Write("GP " + i + ": \n");
            //    var interpolation = element.stdFiniteElement.Interpolation.EvaluateAt(element.stdFiniteElement.Nodes, gp);
            //    ICartesianPoint2D cartesianGP = interpolation.TransformNaturalToCartesian(gp);
            //    for (int enr = 0; enr < 4; ++enr)
            //    {
            //        IEnrichmentFunction2D func = enrichmentItem.EnrichmentFunctions[enr];
            //        EvaluatedFunction2D H = func.EvaluateAllAt(cartesianGP);
            //        Console.Write("Enrichment H" + enr + ": ");
            //        Console.WriteLine("H = " + H.Value + " , dH/dx = " + H.CartesianDerivatives.Item1
            //            + " , dH/dy = " + H.CartesianDerivatives.Item2);
            //    }
            //    Console.WriteLine();
            //}
            //Console.WriteLine();

            //// Check B at GP
            //int[] gpIndices = { 2, 27, 35, 49, 53, 61 };
            //for (int i = 0; i < gpIndices.Length; ++i)
            //{
            //    GaussPoint2D gp = GaussPointsForTipTest.gaussPoints[gpIndices[i]];
            //    Console.WriteLine("GP " + i + ": (x,y,w) = (" + gp.Xi + " , " + gp.Eta + " , " + gp.Weight + ")");
            //    int artDofsCount = element.CountArtificialDofs();
            //    EvaluatedInterpolation2D evaluatedInterpolation =
            //        element.stdFiniteElement.Interpolation.EvaluateAt(element.Nodes, gp);
            //    Matrix2D<double> B = element.CalculateEnrichedDeformationMatrix(artDofsCount, gp, evaluatedInterpolation);
            //    Console.WriteLine(B);
            //    Console.WriteLine();
            //}

            //// Check E, dV at GP
            //for (int i = 0; i < gpIndices.Length; ++i)
            //{
            //    GaussPoint2D gp = GaussPointsForTipTest.gaussPoints[gpIndices[i]];
            //    Console.WriteLine("GP: (x,y,w) = (" + gp.Xi + " , " + gp.Eta + " , " + gp.Weight + ")");
            //    int artDofsCount = element.CountArtificialDofs();
            //    IFiniteElementMaterial2D mat = element.stdFiniteElement.MaterialsOfGaussPoints[gp];
            //    Console.WriteLine(" E = " + mat.CalculateConstitutiveMatrix());
            //    Console.WriteLine();
            //    EvaluatedInterpolation2D evaluatedInterpolation =
            //        element.stdFiniteElement.Interpolation.EvaluateAt(element.Nodes, gp);
            //    double dV = mat.Thickness * evaluatedInterpolation.Jacobian.Determinant * gp.Weight;
            //    Console.WriteLine("dV = " + dV);
            //    Console.WriteLine();
            //}
        }

        private static void IsoparametricQuad4BimaterialTest(XNode2D[] nodes)
        {
            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var materialLeft = ElasticMaterial2DPlainStrain.Create(E, v, t);
            var materialRight = ElasticMaterial2DPlainStrain.Create(E / 2.0, v, t);

            ICurve2D discontinuity = new Line2D(new CartesianPoint2D(30.0, 0.0), new CartesianPoint2D(30.0, 20.0));
            IEnrichmentItem2D enrichmentItem = new MaterialInterface2D(discontinuity);

            var element = XElement2D.CreateBimaterial(new IsoparametricQuad4(nodes), materialLeft, materialRight);

            enrichmentItem.AffectElement(element);
            enrichmentItem.EnrichNodes();

            SymmetricMatrix2D<double> stiffnessStd = element.BuildStdStiffnessMatrix();
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
            //IsoparametricQuad4WithCrackTest(NodeSets.nodeSet7);
            IsoparametricQuad4WithTipTest(NodeSets.nodeSet8);
            //IsoparametricQuad4BimaterialTest(NodeSets.nodeSet7);
        }
    }
}
