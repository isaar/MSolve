using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;

//TODO: The count should also be calculated when ordering them. Perhaps that should be done by the table when an entry is added.
//TODO: Lazy evaluation for constrained dofs or in another class.
//TODO: This assumes dofs are the same. If there are different dof types, more sophisticated Orderers and polymorphism are needed.
//TODO: Use an interface that only defines dofs and interconnection for elements instead of IElementType or ContinuumElement2D
//TODO: While Tables are much clearer than using lists based on the position of nodes or dictionaries of dictionaries, the
//      perfarmance suffers. Benchmark it and find a more efficient way if necessary.
//TODO: Perhaps expose a view of the tables and do not use DofTable.
//TODO: Perhaps the element dof tables or the dof maps should be cached. They are used twice per analysis step.
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// This class is responsible for ordering the freedom degrees and mapping between element and global instances of the same
    /// freedom degrees. Constrained dofs are not explicitly ordered, but assumed to correspond to 0 values for the quantities 
    /// of interest. Even if there are different dof types (e.g. displacements & temperature), they are all treated the same.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class FreeDofOrderer
    {
        private FreeDofOrderer(int numFreeDofs, DofTable<IDof> freeDofs)
        {
            this.NumFreeDofs = numFreeDofs;
            this.FreeDofs = freeDofs;
        }

        public DofTable<IDof> FreeDofs { get; }
        public int NumFreeDofs { get; }

        /// <summary>
        /// Faster than <see cref="CreateWithNodeMajorFreeDofOrder"/>. Use it if the ordering does not matter, e.g. with an
        /// iterative solver (you should probably check that the performance gain of a good ordering is minimal).
        /// </summary>
        public static FreeDofOrderer CreateWithElementMajorFreeDofOrder(IEnumerable<ContinuumElement2D> elements,
            ITable<IDiscretePoint, IDof, double> constraints)
        {
            var freeDofs = new DofTable<IDof>();
            int dofCounter = 0;
            foreach (var element in elements)
            {
                IReadOnlyList<IReadOnlyList<IDof>> elementDofs = element.GetNodalDofs();
                for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
                {
                    IDiscretePoint node = element.Nodes[nodeIdx];
                    for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                    {
                        IDof dof = elementDofs[nodeIdx][dofIdx];
                        if (!constraints.Contains(node, dof)) freeDofs[node, dof] = dofCounter++;
                    }
                }

                // These are for using tables instead of nested lists. They should be purged.
                //DofTable<IDof> elementDofs = element.GetNodalDofsTable();
                //foreach (var node in element.Nodes)
                //{
                //    foreach (var dof in elementDofs.GetColumnsOfRow(node))
                //    {
                //        if (!constraints.Contains(node, dof)) freeDofs[node, dof] = dofCounter++;
                //    }
                //}
            }
            return new FreeDofOrderer(dofCounter, freeDofs);
        }

        /// <summary>
        /// The node major ordering is decent for direct solvers, at least as a starting point for reordering algorithms. 
        /// </summary>
        public static FreeDofOrderer CreateWithNodeMajorFreeDofOrder(IEnumerable<ContinuumElement2D> elements,
            IReadOnlyList<IDiscretePoint> sortedNodes, ITable<IDiscretePoint, IDof, double> constraints)
        {
            var orderer = CreateWithElementMajorFreeDofOrder(elements, constraints);
            ReorderFreeDofsNodeMajor(orderer.FreeDofs, sortedNodes);
            return orderer;
        }

        /// <summary>
        /// Use this only when all nodes have the exact dofs. Its is faster than the other alternatives.
        /// </summary>
        public static FreeDofOrderer CreateWithNodeMajorUniformFreeDofOrder(IReadOnlyList<IDiscretePoint> sortedNodes,
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
            return new FreeDofOrderer(dofCounter, freeDofs);
        }

        public Vector ExtractVectorElementFromGlobal(ContinuumElement2D element, Vector globalVector)
        {
            double[] elementVector = new double[element.GetNodalDofsCount()];
            IReadOnlyList<IReadOnlyList<IDof>> elementDofs = element.GetNodalDofs();
            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
            {
                IDiscretePoint node = element.Nodes[nodeIdx];
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(node, elementDofs[nodeIdx][dofIdx], out int globalDofIdx);
                    if (isFree) elementVector[elementDofIdx] = globalVector[globalDofIdx];
                    // Else, the quantity of interest is 0.0 at all constrained dofs.
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }

            // These are for using tables instead of nested lists. They should be purged.
            //DofTable<IDof> elementDofs = element.GetNodalDofsTable();
            //double[] elementVector = new double[elementDofs.EntryCount];
            //foreach (Tuple<IDiscretePoint, IDof, int> nodeDofNumber in elementDofs)
            //{
            //    bool isFree = FreeDofs.TryGetValue(nodeDofNumber.Item1, nodeDofNumber.Item2, out int globalDof);
            //    if (isFree) elementVector[nodeDofNumber.Item3] = globalVector[globalDof];
            //    // Else, the quantity of interest is 0.0 at all constrained dofs.
            //}

            return Vector.CreateFromArray(elementVector);
        }

        public IReadOnlyDictionary<int, int> MapFreeDofsElementToGlobal(ContinuumElement2D element)
        {
            var dofMap = new Dictionary<int, int>();
            IReadOnlyList<IReadOnlyList<IDof>> elementDofs = element.GetNodalDofs();
            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < element.Nodes.Count; ++nodeIdx)
            {
                IDiscretePoint node = element.Nodes[nodeIdx];
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(node, elementDofs[nodeIdx][dofIdx], out int globalDofIdx);
                    if (isFree) dofMap[elementDofIdx] = globalDofIdx;
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }

            // These are for using tables instead of nested lists. They should be purged.
            //DofTable<IDof> elementDofs = element.GetNodalDofsTable();
            //foreach (Tuple<IDiscretePoint, IDof, int> nodeDofNumber in elementDofs)
            //{
            //    bool isFree = FreeDofs.TryGetValue(nodeDofNumber.Item1, nodeDofNumber.Item2, out int globalDofNumber);
            //    if (isFree) dofMap[nodeDofNumber.Item3] = globalDofNumber;
            //}

            return dofMap;
        }

        /// <summary>
        /// Perhaps this should be called directly, instead of using a dedicated static factory method.
        /// </summary>
        /// <param name="freeDofs"></param>
        /// <param name="sortedNodes"></param>
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
