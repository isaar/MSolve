using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class DofReorder
    {
        /// <summary>
        /// Node major order. Constraint dofs remain labeled as -1.
        /// </summary>
        /// <param name=""></param>
        /// <returns></returns>
        public static int[] OldToNewDofs(Model2D model, int[][] expectedNodalDofs)
        {
            DOFEnumerator enumerator = model.DofEnumerator;
            
            XNode2D[] nodes = model.Nodes.ToArray();
            Array.Sort(nodes); // Nodes must be sorted to produce a node major numbering
            int[][] currentNodalDofs = new int[nodes.Length][];

            for (int n = 0; n < nodes.Length; ++n)
            {
                var dofsOfNode = new List<int>();
                foreach (int dof in enumerator.GetStandardDofsOf(nodes[n]))
                {
                    if (dof != -1) dofsOfNode.Add(dof);
                    else throw new Exception("Does not work if there are constraints");
                }
                foreach (int dof in enumerator.GetArtificialDofsOf(nodes[n])) dofsOfNode.Add(dof);
                currentNodalDofs[n] = dofsOfNode.ToArray();
            }

            int totalDofsCount = 0;
            for (int i = 0; i < expectedNodalDofs.Length; ++i) totalDofsCount += expectedNodalDofs[i].Length;
            int[] oldToNewIndices = new int[totalDofsCount];

            for (int node = 0; node < nodes.Length; ++node)
            {
                Debug.Assert(currentNodalDofs[node].Length == expectedNodalDofs[node].Length, 
                    "Wrong number of dofs for node " + node);
                for (int dof = 0; dof < currentNodalDofs[node].Length; ++dof)
                {
                    oldToNewIndices[currentNodalDofs[node][dof]] = expectedNodalDofs[node][dof];
                }
            }
            return oldToNewIndices;
        }
    }
}
