using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class DofOrderingChecker
    {
        public static void PrintEnumeration(IReadOnlyList<XNode2D> nodes, IDofOrderer dofOrderer, 
            IEnumerable<DisplacementDof> standardDofs, IEnumerable<EnrichedDof> enrichedDofs)
        {
            Console.WriteLine("Standards dofs:");
            foreach (XNode2D node in nodes)
            {
                Console.Write("Node " + node.ID + ": ");
                foreach (DisplacementDof dof in standardDofs)
                {
                    Console.Write(dofOrderer.GetStandardDofOf(node, dof) + " ");
                }
                foreach (int dof in dofOrderer.GetConstrainedDofsOf(node))
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
                foreach (EnrichedDof dof in enrichedDofs)
                {
                    Console.Write(dofOrderer.GetEnrichedDofOf(node, dof) + " ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }
    }
}
