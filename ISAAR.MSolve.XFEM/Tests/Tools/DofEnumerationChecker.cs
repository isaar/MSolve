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
        public static void PrintEnumeration(IReadOnlyList<XNode2D> nodes, IDOFEnumerator dofEnumerator, 
            IEnumerable<DisplacementDOF> standardDofs, IEnumerable<EnrichedDOF> enrichedDofs)
        {
            Console.WriteLine("Standards dofs:");
            foreach (XNode2D node in nodes)
            {
                Console.Write("Node " + node.ID + ": ");
                foreach (DisplacementDOF dof in standardDofs)
                {
                    Console.Write(dofEnumerator.GetFreeDofOf(node, dof) + " ");
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
                foreach (EnrichedDOF dof in enrichedDofs)
                {
                    Console.Write(dofEnumerator.GetEnrichedDofOf(node, dof) + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }
    }
}
