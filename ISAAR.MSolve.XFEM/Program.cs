using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM
{
    class Program
    {
        private static void Quad4Test()
        {
            Node2D[] nodes = new Node2D[4];
            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 4.0, 0.0);
            //nodes[2] = new Node2D(2, 4.0, 3.0);
            //nodes[3] = new Node2D(3, 0.0, 3.0);

            nodes[0] = new Node2D(0, -1.0, -1.0);
            nodes[1] = new Node2D(1, 1.0, -1.0);
            nodes[2] = new Node2D(2, 1.0, 1.0);
            nodes[3] = new Node2D(3, -1.0, 1.0);

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
            //nodes[0] = new Node2D(0, 0.0, 0.0);
            //nodes[1] = new Node2D(1, 4.0, 0.0);
            //nodes[2] = new Node2D(2, 4.0, 3.0);
            //nodes[3] = new Node2D(3, 0.0, 3.0);

            nodes[0] = new Node2D(0, -1.0, -1.0);
            nodes[1] = new Node2D(1, 1.0, -1.0);
            nodes[2] = new Node2D(2, 1.0, 1.0);
            nodes[3] = new Node2D(3, -1.0, 1.0);

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

        static void Main(string[] args)
        {
            Quad4Test();
            IsoparametricQuad4Test();
        }
    }
}
