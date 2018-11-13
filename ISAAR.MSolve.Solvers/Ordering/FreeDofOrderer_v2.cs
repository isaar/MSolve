using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

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
    public class FreeDofOrderer_v2
    {
        private FreeDofOrderer_v2(int numFreeDofs, DofTable_v2 freeDofs)
        {
            this.NumFreeDofs = numFreeDofs;
            this.FreeDofs = freeDofs;
        }

        public DofTable_v2 FreeDofs { get; }
        public int NumFreeDofs { get; }

        /// <summary>
        /// Faster than <see cref="CreateWithNodeMajorFreeDofOrder"/>. Use it if the ordering does not matter, e.g. with an
        /// iterative solver (you should probably check that the performance gain of a good ordering is minimal).
        /// </summary>
        public static FreeDofOrderer_v2 CreateWithElementMajorFreeDofOrder(IEnumerable<IElement> elements, 
            Dictionary<int, Dictionary<DOFType, double>> constraints)
        {
            var freeDofs = new DofTable_v2();
            int dofCounter = 0;
            foreach (IElement element in elements)
            {
                //IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element); //this is wrong
                IList<INode> elementNodes = element.INodes;
                IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element);
                for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
                {
                    bool isNodeConstrained = constraints.TryGetValue(elementNodes[nodeIdx].ID, 
                        out Dictionary<DOFType, double> constraintsOfNode);
                    for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                    {
                        DOFType dofType = elementDofs[nodeIdx][dofIdx];
                        bool isDofConstrained = isNodeConstrained ? constraintsOfNode.ContainsKey(dofType) : false;
                        if (!isDofConstrained)
                        {
                            bool isNewDof = freeDofs.TryAdd(elementNodes[nodeIdx], dofType, dofCounter);
                            if (isNewDof) ++dofCounter;
                        }
                    }
                }
            }
            return new FreeDofOrderer_v2(dofCounter, freeDofs);
        }

        /// <summary>
        /// The node major ordering is decent for direct solvers, at least as a starting point for reordering algorithms. 
        /// </summary>
        public static FreeDofOrderer_v2 CreateWithNodeMajorFreeDofOrder(IEnumerable<IElement> elements,
            IReadOnlyList<INode> sortedNodes, Dictionary<int, Dictionary<DOFType, double>> constraints)
        {
            var orderer = CreateWithElementMajorFreeDofOrder(elements, constraints);
            ReorderFreeDofsNodeMajor(orderer.FreeDofs, sortedNodes);
            return orderer;
        }

        /// <summary>
        /// Use this only when all nodes have the exact dofs. Its is faster than the other alternatives.
        /// </summary>
        public static FreeDofOrderer_v2 CreateWithNodeMajorUniformFreeDofOrder(IReadOnlyList<INode> sortedNodes,
            IReadOnlyList<DOFType> dofsPerNode, Dictionary<int, Dictionary<DOFType, double>> constraints)
        {
            var freeDofs = new DofTable_v2();
            int dofCounter = 0;
            foreach (INode node in sortedNodes)
            {
                bool isNodeConstrained = constraints.TryGetValue(node.ID, out Dictionary<DOFType, double> constraintsOfNode);
                foreach (DOFType dof in dofsPerNode)
                {
                    bool isDofConstrained = isNodeConstrained ? constraintsOfNode.ContainsKey(dof) : false;
                    if (!isDofConstrained) freeDofs[node, dof] = dofCounter++;
                }
            }
            return new FreeDofOrderer_v2(dofCounter, freeDofs);
        }

        public Vector ExtractVectorElementFromGlobal(IElement element, Vector globalVector)
        {
            IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
            IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element);
            int numElementDofs = 0;
            for (int nodeIdx = 0; nodeIdx < elementDofs.Count; ++nodeIdx) numElementDofs += elementDofs[nodeIdx].Count;

            double[] elementVector = new double[numElementDofs];
            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx], 
                        out int globalDofIdx);
                    if (isFree) elementVector[elementDofIdx] = globalVector[globalDofIdx];
                    // Else, the quantity of interest is 0.0 at all constrained dofs.
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
            return Vector.CreateFromArray(elementVector);
        }

        public IReadOnlyDictionary<int, int> MapFreeDofsElementToGlobal(IElement element)
        {
            IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
            //IList<INode> elementNodes = element.INodes;
            //IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element);
            IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypes(element);

            var dofMap = new Dictionary<int, int>();
            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx], 
                        out int globalDofIdx);
                    if (isFree) dofMap[elementDofIdx] = globalDofIdx;
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
            return dofMap;
        }

        /// <summary>
        /// Perhaps this should be called directly, instead of using a dedicated static factory method.
        /// </summary>
        /// <param name="freeDofs"></param>
        /// <param name="sortedNodes"></param>
        private static void ReorderFreeDofsNodeMajor(DofTable_v2 freeDofs, IReadOnlyList<INode> sortedNodes)
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
