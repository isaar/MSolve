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
        public static void PrintEnumeration(IReadOnlyList<XNode2D> nodes, IDOFEnumerator dofEnumerator)
        {

            Console.WriteLine("Standards dofs:");
            foreach (XNode2D node in nodes)
            {
                Console.Write("Node " + node.ID + ": ");
                foreach (int dof in dofEnumerator.GetFreeDofsOf(node))
                {
                    Console.Write(dof + " ");
                }
                foreach (int dof in dofEnumerator.GetConstrainedDofsOf(node))
                {
                    Console.Write("c" + dof + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            Console.WriteLine("Enriched dofs:");
            foreach (XNode2D node in nodes)
            {
                Console.Write("Node " + node.ID + ": ");
                foreach (int dof in dofEnumerator.GetArtificialDofsOf(node))
                {
                    Console.Write(dof + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }
    }
}
