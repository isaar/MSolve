using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Materials;
using ISAAR.MSolve.XFEM.Integration.Strategies;

namespace ISAAR.MSolve.XFEM.Tests.FEM
{
    static class FEMTests
    {
        public static XNode2D[] nodeSet1 = {
            new XNode2D(0, -1.0, -1.0),
            new XNode2D(1, 1.0, -1.0),
            new XNode2D(2, 1.0, 1.0),
            new XNode2D(3, -1.0, 1.0)
        };

        public static XNode2D[] nodeSet2 = {
            new XNode2D(0, -2.0, -2.0),
            new XNode2D(1, 2.0, -2.0),
            new XNode2D(2, 2.0, 2.0),
            new XNode2D(3, -2.0, 2.0)
        };

        public static XNode2D[] nodeSet3 = {
            new XNode2D(0, 0.0, 0.0),
            new XNode2D(1, 2.0, 0.0),
            new XNode2D(2, 2.0, 2.0),
            new XNode2D(3, 0.0, 2.0)
        };

        public static XNode2D[] nodeSet4 = {
            new XNode2D(0, 0.0, 0.0),
            new XNode2D(1, 4.0, 0.0),
            new XNode2D(2, 4.0, 4.0),
            new XNode2D(3, 0.0, 4.0)
        };

        public static XNode2D[] nodeSet5 = {
            new XNode2D(0, 0.0, 0.0),
            new XNode2D(1, 4.0, 0.0),
            new XNode2D(2, 4.0, 3.0),
            new XNode2D(3, 0.0, 3.0)
        };

        public static XNode2D[] nodeSet6A = {
            new XNode2D(0, 0.2, 0.3),
            new XNode2D(1, 2.2, 1.5),
            new XNode2D(2, 3.0, 2.7),
            new XNode2D(3, 0.7, 2.0)
        };

        public static XNode2D[] nodeSet6B = {
            new XNode2D(0, 0.2, 0.3),
            new XNode2D(1, 2.2, 0.9),
            new XNode2D(2, 3.0, 2.7),
            new XNode2D(3, 0.7, 2.0)
        };

        public static XNode2D[] nodeSet7 = {
            new XNode2D(0, 20.0, 0.0),
            new XNode2D(1, 40.0, 0.0),
            new XNode2D(2, 40.0, 20.0),
            new XNode2D(3, 20.0, 20.0)
        };

        public static XNode2D[] nodeSet8 = {
            new XNode2D(0, 20.0, 40.0),
            new XNode2D(1, 40.0, 40.0),
            new XNode2D(2, 40.0, 60.0),
            new XNode2D(3, 20.0, 60.0)
        };

        private static void IsoparametricQuad4Test(XNode2D[] nodes)
        {
            double E = 1;
            double v = 0.25;
            double t = 1.0;
            var material = HomogeneousElasticMaterial2D.CreateMaterialForPlaneStress(E, v, t);
            var element = new ContinuumElement2D(IsoparametricElementType2D.Quad4, nodes, 
                new SimpleIntegration2D(), material);
            Matrix k = element.BuildStiffnessMatrix();
            Console.WriteLine("Isoparametric Quad4 stiffness matrix = ");
            (new FullMatrixWriter()).WriteToConsole(k);
        }

        public static void Main()
        {
            IsoparametricQuad4Test(nodeSet6A);
        }
    }
}
