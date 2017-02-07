using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments;
using ISAAR.MSolve.XFEM.Enrichments.Jump;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM
{
    class Program
    {
        private static void Quad4Test()
        {
            Node2D[] nodes = new Node2D[4];

            //nodes[0] = new Node2D(0, -1.0, -1.0);
            //nodes[1] = new Node2D(1, 1.0, -1.0);
            //nodes[2] = new Node2D(2, 1.0, 1.0);
            //nodes[3] = new Node2D(3, -1.0, 1.0);

            //nodes[0] = new Node2D(0, -2.0, -2.0);
            //nodes[1] = new Node2D(1, 2.0, -2.0);
            //nodes[2] = new Node2D(2, 2.0, 2.0);
            //nodes[3] = new Node2D(3, -2.0, 2.0);

            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 2.0, 0.0);
            //nodes[2] = new Node2D(2, 2.0, 2.0);
            //nodes[3] = new Node2D(3, 0.0, 2.0);

            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 4.0, 0.0);
            //nodes[2] = new Node2D(2, 4.0, 4.0);
            //nodes[3] = new Node2D(3, 0.0, 4.0);

            nodes[0] = new Node2D(0, 0.0, 0.0);
            nodes[1] = new Node2D(1, 4.0, 0.0);
            nodes[2] = new Node2D(2, 4.0, 3.0);
            nodes[3] = new Node2D(3, 0.0, 3.0);

            double E = 1.0;
            double v = 0.25;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStress.Create(E, v, t);

            var element = new Quad4(nodes, material);
            SymmetricMatrix2D<double> k = element.BuildStiffnessMatrix();
            Console.WriteLine("Quad4 stiffness matrix = ");
            Console.WriteLine(k);
        }

        private static void IsoparametricQuad4Test()
        {
            Node2D[] nodes = new Node2D[4];

            //nodes[0] = new Node2D(0, -1.0, -1.0);
            //nodes[1] = new Node2D(1, 1.0, -1.0);
            //nodes[2] = new Node2D(2, 1.0, 1.0);
            //nodes[3] = new Node2D(3, -1.0, 1.0);

            //nodes[0] = new Node2D(0, -2.0, -2.0);
            //nodes[1] = new Node2D(1, 2.0, -2.0);
            //nodes[2] = new Node2D(2, 2.0, 2.0);
            //nodes[3] = new Node2D(3, -2.0, 2.0);

            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 2.0, 0.0);
            //nodes[2] = new Node2D(2, 2.0, 2.0);
            //nodes[3] = new Node2D(3, 0.0, 2.0);

            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 4.0, 0.0);
            //nodes[2] = new Node2D(2, 4.0, 4.0);
            //nodes[3] = new Node2D(3, 0.0, 4.0);

            nodes[0] = new Node2D(0, 0.0, 0.0);
            nodes[1] = new Node2D(1, 4.0, 0.0);
            nodes[2] = new Node2D(2, 4.0, 3.0);
            nodes[3] = new Node2D(3, 0.0, 3.0);

            //nodes[0] = new Node2D(0, 0.2, 0.3);
            //nodes[1] = new Node2D(1, 2.2, 1.5);
            //nodes[2] = new Node2D(2, 3.0, 2.7);
            //nodes[3] = new Node2D(3, 0.7, 2.0);

            double E = 1.0;
            double v = 0.25;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStress.Create(E, v, t);

            var element = new IsoparametricQuad4(nodes, material);
            SymmetricMatrix2D<double> k = element.BuildStiffnessMatrix();
            Console.WriteLine("Isoparametric Quad4 stiffness matrix = ");
            Console.WriteLine(k);
        }

        private static void Quad4WithCrackTest()
        {
            Node2D[] nodes = new Node2D[4];

            nodes[0] = new Node2D(0, 0.0, 0.0);
            nodes[1] = new Node2D(1, 20.0, 0.0);
            nodes[2] = new Node2D(2, 20.0, 20.0);
            nodes[3] = new Node2D(3, 0.0, 20.0);

            //nodes[0] = new Node2D(0, 20.0, 0.0);
            //nodes[1] = new Node2D(1, 40.0, 0.0);
            //nodes[2] = new Node2D(2, 40.0, 20.0);
            //nodes[3] = new Node2D(3, 20.0, 20.0);

            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStrain.Create(E, v, t);

            var element = new Quad4(nodes, material);
            SymmetricMatrix2D<double> k = element.BuildStiffnessMatrix();
            k.Scale(1e-6);
            Console.WriteLine("Quad4 stiffness matrix = ");
            Console.WriteLine(k);
        }

        private static void IsoparametricQuad4WithCrackTest()
        {
            Node2D[] nodes = new Node2D[4];

            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 20.0, 0.0);
            //nodes[2] = new Node2D(2, 20.0, 20.0);
            //nodes[3] = new Node2D(3, 0.0, 20.0);

            nodes[0] = new Node2D(0, 20.0, 0.0);
            nodes[1] = new Node2D(1, 40.0, 0.0);
            nodes[2] = new Node2D(2, 40.0, 20.0);
            nodes[3] = new Node2D(3, 20.0, 20.0);

            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStrain.Create(E, v, t);

            ICurve2D discontinuity = new Line2D(new Point2D(30.0, 0.0), new Point2D(30.0, 20.0));
            IInterfaceEnrichment enrichmentFunction = new Sign();
            //IInterfaceEnrichment enrichmentFunction = new PolynomialSmoothedHeaviside(0.1);

            var element = new IsoparametricQuad4WithDiscontinuity(nodes, material, discontinuity, enrichmentFunction);

            SymmetricMatrix2D<double> stiffnessStd = element.BuildStdStiffnessMatrix();
            Matrix2D<double> stiffnessStdEnr;
            SymmetricMatrix2D<double> stiffnessEnr;
            element.BuildEnrichedStiffnessMatrices(out stiffnessStdEnr, out stiffnessEnr);

            stiffnessStd.Scale(1e-6);
            stiffnessStdEnr.Scale(1e-6);
            stiffnessEnr.Scale(1e-6);

            Console.WriteLine("Quad4 standard-standard stiffness matrix = ");
            Console.WriteLine(stiffnessStd);
            Console.WriteLine("Quad4 stdandard-enriched stiffness matrix = ");
            Console.WriteLine(stiffnessStdEnr);
            Console.WriteLine("Quad4 enriched-enriched stiffness matrix = ");
            Console.WriteLine(stiffnessEnr);
        }

        private static void IsoparametricQuad4BimaterialTest()
        {
            Node2D[] nodes = new Node2D[4];

            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 20.0, 0.0);
            //nodes[2] = new Node2D(2, 20.0, 20.0);
            //nodes[3] = new Node2D(3, 0.0, 20.0);

            nodes[0] = new Node2D(0, 20.0, 0.0);
            nodes[1] = new Node2D(1, 40.0, 0.0);
            nodes[2] = new Node2D(2, 40.0, 20.0);
            nodes[3] = new Node2D(3, 20.0, 20.0);

            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var materialLeft = ElasticMaterial2DPlainStrain.Create(E, v, t); //TODO: need a way to supply materials with predicates concerning (x, y) or (ksi, eta)
            var materialRight = ElasticMaterial2DPlainStrain.Create(E / 2.0, v, t);

            ICurve2D discontinuity = new Line2D(new Point2D(30.0, 0.0), new Point2D(30.0, 20.0));
            IInterfaceEnrichment enrichmentFunction = new Ramp();


            var element = new IsoparametricQuad4WithDiscontinuity(nodes, materialLeft, materialRight, discontinuity, enrichmentFunction);

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

        private static void EnrichedIsoparametricQuad4WithCrackTest()
        {
            Node2D[] nodes = new Node2D[4];

            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 20.0, 0.0);
            //nodes[2] = new Node2D(2, 20.0, 20.0);
            //nodes[3] = new Node2D(3, 0.0, 20.0);

            nodes[0] = new Node2D(0, 20.0, 0.0);
            nodes[1] = new Node2D(1, 40.0, 0.0);
            nodes[2] = new Node2D(2, 40.0, 20.0);
            nodes[3] = new Node2D(3, 20.0, 20.0);

            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStrain.Create(E, v, t);

            ICurve2D discontinuity = new Line2D(new Point2D(30.0, 0.0), new Point2D(30.0, 20.0));
            IEnrichmentFunction2D enrichmentFunction = new HeavisideEnrichment2D(discontinuity);

            var element = new IsoparametricQuad4Enriched(nodes, material, enrichmentFunction);

            SymmetricMatrix2D<double> stiffnessStd = element.BuildStdStiffnessMatrix();
            Matrix2D<double> stiffnessStdEnr;
            SymmetricMatrix2D<double> stiffnessEnr;
            element.BuildEnrichedStiffnessMatrices(out stiffnessStdEnr, out stiffnessEnr);

            stiffnessStd.Scale(1e-6);
            stiffnessStdEnr.Scale(1e-6);
            stiffnessEnr.Scale(1e-6);

            Console.WriteLine("Quad4 standard-standard stiffness matrix = ");
            Console.WriteLine(stiffnessStd);
            Console.WriteLine("Quad4 stdandard-enriched stiffness matrix = ");
            Console.WriteLine(stiffnessStdEnr);
            Console.WriteLine("Quad4 enriched-enriched stiffness matrix = ");
            Console.WriteLine(stiffnessEnr);
        }

        private static void EnrichedIsoparametricQuad4BimaterialTest()
        {
            Node2D[] nodes = new Node2D[4];

            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 20.0, 0.0);
            //nodes[2] = new Node2D(2, 20.0, 20.0);
            //nodes[3] = new Node2D(3, 0.0, 20.0);

            nodes[0] = new Node2D(0, 20.0, 0.0);
            nodes[1] = new Node2D(1, 40.0, 0.0);
            nodes[2] = new Node2D(2, 40.0, 20.0);
            nodes[3] = new Node2D(3, 20.0, 20.0);

            double E = 2e6;
            double v = 0.3;
            double t = 1.0;
            var materialLeft = ElasticMaterial2DPlainStrain.Create(E, v, t); //TODO: need a way to supply materials with predicates concerning (x, y) or (ksi, eta)
            var materialRight = ElasticMaterial2DPlainStrain.Create(E / 2.0, v, t);

            ICurve2D discontinuity = new Line2D(new Point2D(30.0, 0.0), new Point2D(30.0, 20.0));
            IEnrichmentFunction2D enrichmentFunction = new RampEnrichment2D(discontinuity);


            var element = new IsoparametricQuad4Enriched(nodes, materialLeft, materialRight, enrichmentFunction);

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
            //Quad4Test();
            //IsoparametricQuad4Test();

            //Quad4WithCrackTest();
            //IsoparametricQuad4WithCrackTest();
            //EnrichedIsoparametricQuad4WithCrackTest();
            IsoparametricQuad4BimaterialTest();
            EnrichedIsoparametricQuad4BimaterialTest();
        }
    }
}
