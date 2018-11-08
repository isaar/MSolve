using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    class DofReorder
    {
        /// <summary>
        /// Node major order. Constraint dofs remain labeled as -1.
        /// </summary>
        /// <param name=""></param>
        /// <returns></returns>
        public static int[] OldToNewDofs(Model2D model, int[][] expectedNodalDofs, IDofOrderer dofOrderer)
        {
            XNode2D[] nodes = model.Nodes.ToArray();
            Array.Sort(nodes); // Nodes must be sorted to produce a node major numbering
            int[][] currentNodalDofs = new int[nodes.Length][];

            for (int n = 0; n < nodes.Length; ++n)
            {
                var dofsOfNode = new List<int>();
                foreach (int dof in dofOrderer.GetStandardDofsOf(nodes[n]))
                {
                    dofsOfNode.Add(dof);
                    //if (dof != -1) dofsOfNode.Add(dof);
                    //else throw new Exception("Does not work if there are constraints");
                }
                foreach (int dof in dofOrderer.GetEnrichedDofsOf(nodes[n])) dofsOfNode.Add(dof);
                currentNodalDofs[n] = dofsOfNode.ToArray();
            }

            int freeDofsCount = CountStandardDofs(expectedNodalDofs);
            int[] oldToNewIndices = new int[freeDofsCount];

            for (int node = 0; node < nodes.Length; ++node)
            {
                Debug.Assert(currentNodalDofs[node].Length == expectedNodalDofs[node].Length, 
                    "Wrong number of dofs for node " + node);
                for (int d = 0; d < currentNodalDofs[node].Length; ++d)
                {
                    int currentDof = currentNodalDofs[node][d];
                    int expectedDof = expectedNodalDofs[node][d];
                    if ((currentDof < 0) || (expectedDof < 0))
                    {
                        if (currentDof != expectedDof)
                        {
                            throw new Exception("The " + d + "th dof of node " + node +
                                " is constrained only on one of the 2 models that are compared.");
                        }
                    }
                    else oldToNewIndices[currentDof] = expectedDof;
                }
            }
            return oldToNewIndices;
        }

        private static int CountStandardDofs(int[][] expectedNodalDofs)
        {
            int count = 0;
            foreach (int[] nodeDofs in expectedNodalDofs)
            {
                foreach (int dof in nodeDofs)
                {
                    if (dof >= 0) ++count;
                }
            }
            return count;
        }
    }
}
