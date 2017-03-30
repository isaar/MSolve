using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.LinearAlgebra;

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
            SymmetricMatrix2D k = element.BuildStiffnessMatrix();
            Console.WriteLine("Quad4 stiffness matrix = ");
            MatrixUtilities.PrintDense(k);
        }

        private static void IsoparametricQuad4Test(XNode2D[] nodes)
        {
            double E = 1;
            double v = 0.25;
            double t = 1.0;
            var material = ElasticMaterial2DPlainStress.Create(E, v, t);
            var element = new ContinuumElement2D(IsoparametricElementType2D.Quad4, nodes, new SimpleIntegration2D(material));
            SymmetricMatrix2D k = element.BuildStiffnessMatrix();
            Console.WriteLine("Isoparametric Quad4 stiffness matrix = ");
            MatrixUtilities.PrintDense(k);
        }

        public static void Main()
        {
            //Quad4Test(NodeSets.nodeSet1);
            IsoparametricQuad4Test(NodeSets.nodeSet1);
        }
    }
}
