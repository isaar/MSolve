using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Integration.Strategies;

namespace ISAAR.MSolve.XFEM.Tests
{
    static class FEMTests
    {
        private static void Quad4Test(XNode2D[] nodes)
        {
            double E = 1.0;
            double v = 0.25;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStress.Create(E, v, t);

            var element = new Quad4(nodes, material);
            SymmetricMatrix2D<double> k = element.BuildStiffnessMatrix();
            Console.WriteLine("Quad4 stiffness matrix = ");
            Console.WriteLine(k);
        }

        private static void IsoparametricQuad4Test(XNode2D[] nodes)
        {
            double E = 1;
            double v = 0.25;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStress.Create(E, v, t);
            IIntegrationStrategyFactory2D integrationFactory = new SimpleIntegration2D.Factory(material);

            var element = new ContinuumElement2D(IsoparametricElementType2D.Quad4, nodes, integrationFactory);
            SymmetricMatrix2D<double> k = element.BuildStiffnessMatrix();
            Console.WriteLine("Isoparametric Quad4 stiffness matrix = ");
            Console.WriteLine(k);
        }

        public static void Main()
        {
            //Quad4Test(NodeSets.nodeSet1);
            IsoparametricQuad4Test(NodeSets.nodeSet1);
        }
    }
}
