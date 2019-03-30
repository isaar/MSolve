using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: This works for other dual solvers, not only FETI-1
//TODO: This should only calculate them. Another object should manage them.
//TODO: The enumation code is quite complex and error prone. It should be simplified and decomposed into smaller methods, as
//      much as possible without sacrificing performance.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    /// <summary>
    /// Calculates the signed boolean matrices of the equations that enforce continuity between the multiple instances of 
    /// boundary dofs.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class LagrangeMultipliersEnumerator
    {
        private readonly ICrosspointStrategy crosspointStrategy;

        internal LagrangeMultipliersEnumerator(ICrosspointStrategy crosspointStrategy)
        {
            this.crosspointStrategy = crosspointStrategy;
        }

        public Dictionary<int, SignedBooleanMatrix> BooleanMatrices { get; private set; }

        //TODO: I am not too thrilled about objects with properties that may or may not be null
        /// <summary>
        /// WARNING: This property will be null in homogeneous problems.
        /// Associates each lagrange multiplier with the instances of the boundary dof, for which continuity is enforced. 
        /// </summary>
        internal LagrangeMultiplier[] LagrangeMultipliers { get; private set; }

        public int NumLagrangeMultipliers { get; private set; }

        /// <summary>
        /// Only <see cref="BooleanMatrices"/> will be explicitly created. <see cref="LagrangeMultipliers"/> will not.
        /// For use in homogeneous problems, where we do not need that much info about lagrange multipliers and boundary dofs.
        /// </summary>
        /// <param name="model"></param>
        public void DefineBooleanMatrices(IStructuralModel_v2 model, IDofSeparator dofSeparator)
        {
            List<(INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)> boundaryNodeData 
                = DefineBoundaryDofConstraints(dofSeparator);

            InitializeBooleanMatrices(model);

            // Fill the boolean matrices: node major, subdomain medium, dof minor. TODO: not sure about this order.
            int lag = 0; // Lagrange multiplier index
            foreach ((INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)
                in boundaryNodeData)
            {
                int numSubdomainCombos = subdomainsPlus.Length;
                for (int c = 0; c < numSubdomainCombos; ++c)
                {
                    //TODO: each subdomain appears in many combinations. It would be faster to cache its indices for the specific dofs.
                    SignedBooleanMatrix booleanPlus = BooleanMatrices[subdomainsPlus[c].ID];
                    SignedBooleanMatrix booleanMinus = BooleanMatrices[subdomainsMinus[c].ID];
                    IReadOnlyDictionary<DOFType, int> dofsPlus = subdomainsPlus[c].FreeDofOrdering.FreeDofs.GetDataOfRow(node);
                    IReadOnlyDictionary<DOFType, int> dofsMinus = subdomainsMinus[c].FreeDofOrdering.FreeDofs.GetDataOfRow(node);

                    foreach (DOFType dof in dofs)
                    {
                        booleanPlus.AddEntry(lag, dofsPlus[dof], true);
                        booleanMinus.AddEntry(lag, dofsMinus[dof], false);
                        ++lag;
                    }
                }
            }
        }

        /// <summary>
        /// Creates both <see cref="BooleanMatrices"/> and <see cref="LagrangeMultipliers"/>. For use in heterogeneous problems.
        /// </summary>
        /// <param name="model"></param>
        public void DefineLagrangesAndBooleanMatrices(IStructuralModel_v2 model, IDofSeparator dofSeparator)
        {
            List<(INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)> boundaryNodeData
                = DefineBoundaryDofConstraints(dofSeparator);

            InitializeBooleanMatrices(model);
            LagrangeMultipliers = new LagrangeMultiplier[NumLagrangeMultipliers];

            // Fill the boolean matrices and lagrange multiplier data: node major, subdomain medium, dof minor. TODO: not sure about this order.
            int lag = 0; // Lagrange multiplier index
            foreach ((INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)
                in boundaryNodeData)
            {
                int numSubdomainCombos = subdomainsPlus.Length;
                for (int c = 0; c < numSubdomainCombos; ++c)
                {
                    //TODO: each subdomain appears in many combinations. It would be faster to cache its indices for the specific dofs.
                    SignedBooleanMatrix booleanPlus = BooleanMatrices[subdomainsPlus[c].ID];
                    SignedBooleanMatrix booleanMinus = BooleanMatrices[subdomainsMinus[c].ID];

                    //TODO: The dof indices have already been accessed. Reuse it if possible.
                    IReadOnlyDictionary<DOFType, int> dofsPlus = subdomainsPlus[c].FreeDofOrdering.FreeDofs.GetDataOfRow(node);
                    IReadOnlyDictionary<DOFType, int> dofsMinus = subdomainsMinus[c].FreeDofOrdering.FreeDofs.GetDataOfRow(node);

                    foreach (DOFType dof in dofs)
                    {
                        booleanPlus.AddEntry(lag, dofsPlus[dof], true);
                        booleanMinus.AddEntry(lag, dofsMinus[dof], false);
                        LagrangeMultipliers[lag] = new LagrangeMultiplier(node, dof, subdomainsPlus[c], subdomainsMinus[c]);
                        ++lag;
                    }
                }
            }
        }
        
        //TODO: Perhaps this should return the number of lagranges, instead of setting it. It is unexpected.
        private List<(INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)> 
            DefineBoundaryDofConstraints(IDofSeparator dofSeparator)
        {
            // Find boundary nodes and dofs
            var boundaryNodeData = new List<(INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, 
                ISubdomain_v2[] subdomainsMinus)>(dofSeparator.GlobalBoundaryDofs.Count);

            // Find continuity equations.
            NumLagrangeMultipliers = 0;
            foreach (var nodeDofsPair in dofSeparator.GlobalBoundaryDofs)
            {
                INode node = nodeDofsPair.Key;
                DOFType[] dofsOfNode = nodeDofsPair.Value;
                IEnumerable <ISubdomain_v2> nodeSubdomains = node.SubdomainsDictionary.Values;
                (ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus) =
                    crosspointStrategy.FindSubdomainCombinations(nodeSubdomains.Count(), nodeSubdomains);

                boundaryNodeData.Add((node, dofsOfNode, subdomainsPlus, subdomainsMinus));
                NumLagrangeMultipliers += dofsOfNode.Length * subdomainsPlus.Length;
            }

            //// Find boundary nodes and continuity equations.
            //var boundaryNodeData =
            //    new List<(INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)>();
            //NumLagrangeMultipliers = 0;
            //foreach (INode node in model.Nodes) //TODO: this probably doesn't work if there are embedded nodes. It is time to isolate the embedded nodes.
            //{
            //    int nodeMultiplicity = node.SubdomainsDictionary.Count;
            //    if (nodeMultiplicity > 1)
            //    {
            //        DOFType[] dofsOfNode = model.GlobalDofOrdering.GlobalFreeDofs.GetColumnsOfRow(node).ToArray(); //TODO: interacting with the table can be optimized.

            //        // If all dofs of this node are constrained, then it is not considered boundary.
            //        if (dofsOfNode.Length == 0) continue;

            //        IEnumerable<ISubdomain_v2> nodeSubdomains = node.SubdomainsDictionary.Values;
            //        (ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus) =
            //            crosspointStrategy.FindSubdomainCombinations(nodeMultiplicity, nodeSubdomains);

            //        boundaryNodeData.Add((node, dofsOfNode, subdomainsPlus, subdomainsMinus));
            //        NumLagrangeMultipliers += dofsOfNode.Length * subdomainsPlus.Length;
            //    }
            //}

            return boundaryNodeData;
        }

        private void InitializeBooleanMatrices(IStructuralModel_v2 model)
        {
            // Create the signed boolean matrices.
            BooleanMatrices = new Dictionary<int, SignedBooleanMatrix>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                BooleanMatrices[subdomain.ID] =
                    new SignedBooleanMatrix(NumLagrangeMultipliers, subdomain.FreeDofOrdering.NumFreeDofs);
            }
        }

        ////TODO: Perhaps optimizations with arrays or a dedicated DTO class are possible
        ///// <summary>
        ///// WARNING: This property will be null in homogeneous problems.
        ///// Associates each lagrange multiplier with the instances of the boundary dof, for which continuity is enforced. These 
        ///// instances are described as a Dictionary with keys the subdomains containing the dof and values the indices of the dof
        ///// in the vectors and matrices of that subdomain. 
        ///// node i
        ///// E.g. to sum the stiffnesses of the dofs associated with 4th lagrange multiplier:
        ///// double totalStiffness = 0.0;
        ///// foreach (var subdomainDofIdxPair in <see cref="LagrangeMultipliersDofRelations"/>[3]) 
        ///// { 
        ///// ISubdomain subdomain = subdomainDofIdxPair.Key;
        ///// int dofIdx = subdomainDofIdxPair.Value;
        ///// totalStiffness += linearSystems[subdomain.ID].Matrix[dofIdx, dofIdx];
        ///// }
        ///// </summary>
        //public Dictionary<ISubdomain_v2, int>[] LagrangeMultipliersDofRelations { get; private set; }

        //private void FillHeterogeneousBooleanMatricesAndAssociateLagrangesWithDofs(
        //    List<(INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)> boundaryNodeData)
        //{
        //    // For each boundary dofs, aggregate the subdomains in which there is an instance and its index in those subdomains' 
        //    // vectors / matrices
        //    var dofAggregation = new Table<INode, DOFType, Dictionary<ISubdomain_v2, int>>();
        //    foreach ((INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)
        //        in boundaryNodeData)
        //    {
        //        IEnumerable<ISubdomain_v2> nodeSubdomains = node.SubdomainsDictionary.Values;
        //        foreach (DOFType dofType in dofs)
        //        {
        //            var dofSubdomainIndices = new Dictionary<ISubdomain_v2, int>();
        //            foreach (ISubdomain_v2 subdomain in nodeSubdomains)
        //            {
        //                dofSubdomainIndices[subdomain] = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
        //            }
        //            dofAggregation[node, dofType] = dofSubdomainIndices;
        //        }
        //    }

        //    // Fill the boolean matrices and lagrange multiplier data: node major, subdomain medium, dof minor. TODO: not sure about this order.
        //    LagrangeMultipliersDofRelations = new Dictionary<ISubdomain_v2, int>[NumLagrangeMultipliers];
        //    int equationCounter = 0;
        //    foreach ((INode node, DOFType[] dofs, ISubdomain_v2[] subdomainsPlus, ISubdomain_v2[] subdomainsMinus)
        //        in boundaryNodeData)
        //    {
        //        int numSubdomainCombos = subdomainsPlus.Length;
        //        for (int c = 0; c < numSubdomainCombos; ++c)
        //        {
        //            //TODO: each subdomain appears in many combinations. It would be faster to cache its indices for the specific dofs.
        //            SignedBooleanMatrix booleanPlus = BooleanMatrices[subdomainsPlus[c].ID];
        //            SignedBooleanMatrix booleanMinus = BooleanMatrices[subdomainsMinus[c].ID];

        //            //TODO: The dof indices have already been accessed. Reuse it if possible.
        //            bool nodeIsOrdered = subdomainsPlus[c].FreeDofOrdering.FreeDofs.TryGetDataOfRow(node,
        //                out IReadOnlyDictionary<DOFType, int> dofsPlus);
        //            Debug.Assert(nodeIsOrdered);
        //            nodeIsOrdered = subdomainsMinus[c].FreeDofOrdering.FreeDofs.TryGetDataOfRow(node,
        //                out IReadOnlyDictionary<DOFType, int> dofsMinus);
        //            Debug.Assert(nodeIsOrdered);

        //            foreach (DOFType dof in dofs)
        //            {
        //                booleanPlus.AddEntry(equationCounter, dofsPlus[dof], true);
        //                booleanMinus.AddEntry(equationCounter, dofsMinus[dof], false);
        //                LagrangeMultipliersDofRelations[equationCounter] = dofAggregation[node, dof];
        //                ++equationCounter;
        //            }
        //        }
        //    }
        //}
    }
}
