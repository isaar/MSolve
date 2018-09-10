using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;

//TODO: remove duplicate code between this class and FreeDofOrderer.
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// This class is responsible for ordering the freedom degrees and mapping between element and global instances of the same
    /// freedom degrees. Constrained dofs are explicitly ordered, in order to build free and constrained global vectors/matrices.
    /// Even if there are different dof types (e.g. displacements & temperature), they are all treated the same.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class AllDofOrderer
    {
        private AllDofOrderer(int numFreeDofs, DofTable<IDof> freeDofs, int numConstrainedDofs, DofTable<IDof> constrainedDofs)
        {
            this.NumFreeDofs = numFreeDofs;
            this.FreeDofs = freeDofs;
            this.NumConstrainedDofs = numConstrainedDofs;
            this.ConstrainedDofs = constrainedDofs;
        }

        public DofTable<IDof> ConstrainedDofs { get; }
        public DofTable<IDof> FreeDofs { get; }
        public int NumConstrainedDofs { get; }
        public int NumFreeDofs { get; }

        /// <summary>
        /// Faster than <see cref="CreateWithNodeMajorFreeDofOrder"/>. Use it if the ordering does not matter, e.g. with an
        /// iterative solver (you should probably check that the performance gain of a good ordering is minimal).
        /// </summary>
        public static AllDofOrderer CreateWithElementMajorFreeDofOrder(IEnumerable<ContinuumElement2D> elements,
            ITable<IDiscretePoint, IDof, double> constraints)
        {
            var freeDofs = new DofTable<IDof>();
            int dofCounter = 0;
            foreach (var element in elements)
            {
                DofTable<IDof> elementDofs = element.GetNodalDofs();
                foreach (var node in element.Nodes)
                {
                    foreach (var dof in elementDofs.GetColumnsOfRow(node))
                    {
                        if (!constraints.Contains(node, dof)) freeDofs[node, dof] = dofCounter++;
                    }
                }
            }

            (int numConstrainedDofs, DofTable<IDof> constrainedDofs) = OrderConstrainedDofsNodeMajor(constraints);
            return new AllDofOrderer(dofCounter, freeDofs, numConstrainedDofs, constrainedDofs);
        }

        /// <summary>
        /// The node major ordering is decent for direct solvers, at least as a starting point for reordering algorithms. 
        /// </summary>
        public static AllDofOrderer CreateWithNodeMajorFreeDofOrder(IEnumerable<ContinuumElement2D> elements,
            IReadOnlyList<IDiscretePoint> sortedNodes, ITable<IDiscretePoint, IDof, double> constraints)
        {
            var orderer = CreateWithElementMajorFreeDofOrder(elements, constraints);
            ReorderFreeDofsNodeMajor(orderer.FreeDofs, sortedNodes);
            return orderer;
        }

        /// <summary>
        /// Use this only when all nodes have the exact dofs. Its is faster than the other alternatives.
        /// </summary>
        public static AllDofOrderer CreateWithNodeMajorUniformFreeDofOrder(IReadOnlyList<IDiscretePoint> sortedNodes,
            IReadOnlyList<IDof> dofsPerNode, ITable<IDiscretePoint, IDof, double> constraints)
        {
            var freeDofs = new DofTable<IDof>();
            int dofCounter = 0;
            foreach (var node in sortedNodes)
            {
                foreach (var dof in dofsPerNode)
                {
                    if (!constraints.Contains(node, dof)) freeDofs[node, dof] = dofCounter++;
                }
            }

            (int numConstrainedDofs, DofTable<IDof> constrainedDofs) = OrderConstrainedDofsNodeMajor(constraints);
            return new AllDofOrderer(dofCounter, freeDofs, numConstrainedDofs, constrainedDofs);
        }

        public Vector ExtractVectorElementFromGlobal(ContinuumElement2D element, Vector globalFreeVector,
            Vector globalConstrainedVector)
        {
            DofTable<IDof> elementDofs = element.GetNodalDofs();
            double[] elementVector = new double[elementDofs.EntryCount];
            foreach (Tuple<IDiscretePoint, IDof, int> nodeDofNumber in elementDofs)
            {
                bool isFree = FreeDofs.TryGetValue(nodeDofNumber.Item1, nodeDofNumber.Item2, out int globalFreeDof);
                if (isFree) elementVector[nodeDofNumber.Item3] = globalFreeVector[globalFreeDof];
                else
                {
                    int globalConstrainedDof = ConstrainedDofs[nodeDofNumber.Item1, nodeDofNumber.Item2];
                    elementVector[nodeDofNumber.Item3] = globalConstrainedVector[globalConstrainedDof];
                }
            }
            return Vector.CreateFromArray(elementVector);
        }

        public (IReadOnlyDictionary<int, int> freeDofsMap, IReadOnlyDictionary<int, int> constrainedDofsMap)
            MapDofsElementToGlobal(ContinuumElement2D element)
        {
            DofTable<IDof> elementDofs = element.GetNodalDofs();
            var freeDofMap = new Dictionary<int, int>();
            var constrainedDofMap = new Dictionary<int, int>();
            foreach (Tuple<IDiscretePoint, IDof, int> nodeDofNumber in elementDofs)
            {
                bool isFree = FreeDofs.TryGetValue(nodeDofNumber.Item1, nodeDofNumber.Item2, out int globalDofNumber);
                if (isFree) freeDofMap[nodeDofNumber.Item3] = globalDofNumber;
                else constrainedDofMap[nodeDofNumber.Item3] = globalDofNumber;
            }
            return (freeDofMap, constrainedDofMap);
        }

        private static (int numConstrainedDofs, DofTable<IDof> constrainedDofs) OrderConstrainedDofsNodeMajor(
            ITable<IDiscretePoint, IDof, double> constraints)
        {
            var constrainedDofs = new DofTable<IDof>();
            int dofCounter = 0;
            foreach (Tuple<IDiscretePoint, IDof, double> nodeDofValue in constraints)
            {
                constrainedDofs[nodeDofValue.Item1, nodeDofValue.Item2] = dofCounter++;
            }
            return (dofCounter, constrainedDofs);
        }

        private static void ReorderFreeDofsNodeMajor(DofTable<IDof> freeDofs, IReadOnlyList<IDiscretePoint> sortedNodes)
        {
            int dofCounter = 0;
            foreach (var node in sortedNodes)
            {
                foreach (var dof in freeDofs.GetColumnsOfRow(node))
                {
                    freeDofs[node, dof] = dofCounter++;
                }
            }
        }
    }
}
