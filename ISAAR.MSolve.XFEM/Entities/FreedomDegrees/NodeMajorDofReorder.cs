using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Entities.FreedomDegrees
{
    static class NodeMajorDofReorder
    {
        /// <summary>
        /// Node major order. Constraint dofs remain labeled as -1.
        /// </summary>
        /// <param name=""></param>
        /// <returns></returns>
        public static int[] OldToNewDofs(Node2D[] nodes, IDOFEnumerator enumerator)
        {
            int totalDofsCount = enumerator.FreeDofsCount + enumerator.ArtificialDofsCount;

            // In case some dof is mistakenly skipped, initialize the array to an illegal value, such as -2 (-1 is taken)
            List<int> permutation = Enumerable.Repeat(-2, totalDofsCount).ToList();

            Array.Sort(nodes); // Nodes must be sorted to produce a node major numbering

            int newIndex = 0;
            foreach (XNode2D node in nodes)
            {
                var stdDofs = enumerator.GetFreeDofsOf(node);
                var enrDofs = enumerator.GetArtificialDofsOf(node);
                foreach (int dof in stdDofs)
                {
                    if (dof !=- 1) permutation[dof] = newIndex++;
                }
                foreach (int dof in enrDofs) permutation[dof] = newIndex++;
            }
            return permutation.GetRange(0, newIndex).ToArray();
        }
    }
}
