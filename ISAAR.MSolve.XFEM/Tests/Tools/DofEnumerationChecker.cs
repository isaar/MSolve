using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class DofEnumerationChecker
    {
        public static void PrintEnumeration(Model2D model)
        {
            DOFEnumerator enumerator = model.DofEnumerator;

            Console.WriteLine("Standards dofs:");
            foreach (XNode2D node in model.Nodes)
            {
                Console.Write("Node " + node.ID + ": ");
                foreach (int dof in enumerator.GetStandardDofsOf(node))
                {
                    Console.Write(dof + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            Console.WriteLine("Enriched dofs:");
            foreach (XNode2D node in model.Nodes)
            {
                Console.Write("Node " + node.ID + ": ");
                foreach (int dof in enumerator.GetArtificialDofsOf(node))
                {
                    Console.Write(dof + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }
    }
}
