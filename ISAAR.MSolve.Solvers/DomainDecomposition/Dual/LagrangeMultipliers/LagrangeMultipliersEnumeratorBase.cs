using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;

//TODO: This should only calculate them. Another object should manage them.
//TODO: The enumeration code is quite complex and error prone. It should be simplified and decomposed into smaller methods, as
//      much as possible without sacrificing performance.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers
{
    /// <summary>
    /// Calculates the signed boolean matrices of the equations that enforce continuity between the multiple instances of 
    /// boundary dofs.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class LagrangeMultipliersEnumeratorBase : ILagrangeMultipliersEnumerator
    {
        private readonly ICrosspointStrategy crosspointStrategy;
        private readonly IDofSeparator dofSeparator;

        protected LagrangeMultipliersEnumeratorBase(ICrosspointStrategy crosspointStrategy, IDofSeparator dofSeparator)
        {
            this.crosspointStrategy = crosspointStrategy;
            this.dofSeparator = dofSeparator;
        }

        public Dictionary<int, SignedBooleanMatrixColMajor> BooleanMatrices { get; private set; }

        //TODO: I am not too thrilled about objects with properties that may or may not be null
        /// <summary>
        /// WARNING: This property will be null in homogeneous problems.
        /// Associates each lagrange multiplier with the instances of the boundary dof, for which continuity is enforced. 
        /// </summary>
        public LagrangeMultiplier[] LagrangeMultipliers { get; private set; }

        public int NumLagrangeMultipliers { get; private set; }

        /// <summary>
        /// Only <see cref="BooleanMatrices"/> will be explicitly created. <see cref="LagrangeMultipliers"/> will not.
        /// For use in homogeneous problems, where we do not need that much info about lagrange multipliers and boundary dofs.
        /// </summary>
        /// <param name="model"></param>
        protected void DefineBooleanMatrices(IStructuralModel model, 
            Dictionary<int, int> numRemainderDofs, Dictionary<int, DofTable> remainderDofOrderings) //TODO: Rename the "remainder"
        {
            List<(INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)> boundaryNodeData 
                = DefineBoundaryDofConstraints(dofSeparator);

            // Create the signed boolean matrices.
            BooleanMatrices = new Dictionary<int, SignedBooleanMatrixColMajor>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                BooleanMatrices[subdomain.ID] =
                    new SignedBooleanMatrixColMajor(NumLagrangeMultipliers, numRemainderDofs[subdomain.ID]);
            }

            // Fill the boolean matrices: node major, subdomain medium, dof minor. TODO: not sure about this order.
            int lag = 0; // Lagrange multiplier index
            foreach ((INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)
                in boundaryNodeData)
            {
                int numSubdomainCombos = subdomainsPlus.Length;
                for (int c = 0; c < numSubdomainCombos; ++c)
                {
                    //TODO: each subdomain appears in many combinations. It would be faster to cache its indices for the specific dofs.
                    SignedBooleanMatrixColMajor booleanPlus = BooleanMatrices[subdomainsPlus[c].ID];
                    SignedBooleanMatrixColMajor booleanMinus = BooleanMatrices[subdomainsMinus[c].ID];
                    IReadOnlyDictionary<IDofType, int> dofsPlus = remainderDofOrderings[subdomainsPlus[c].ID].GetDataOfRow(node);
                    IReadOnlyDictionary<IDofType, int> dofsMinus = remainderDofOrderings[subdomainsMinus[c].ID].GetDataOfRow(node);

                    foreach (IDofType dof in dofs)
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
        protected void DefineLagrangesAndBooleanMatrices(IStructuralModel model,
            Dictionary<int, int> numRemainderDofs, Dictionary<int, DofTable> remainderDofOrderings) //TODO: Rename the "remainder"
        {
            List<(INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)> boundaryNodeData
                = DefineBoundaryDofConstraints(dofSeparator);

            // Create the signed boolean matrices.
            BooleanMatrices = new Dictionary<int, SignedBooleanMatrixColMajor>();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                BooleanMatrices[subdomain.ID] =
                    new SignedBooleanMatrixColMajor(NumLagrangeMultipliers, numRemainderDofs[subdomain.ID]);
            }
            LagrangeMultipliers = new LagrangeMultiplier[NumLagrangeMultipliers];

            // Fill the boolean matrices and lagrange multiplier data: node major, subdomain medium, dof minor. TODO: not sure about this order.
            int lag = 0; // Lagrange multiplier index
            foreach ((INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)
                in boundaryNodeData)
            {
                int numSubdomainCombos = subdomainsPlus.Length;
                for (int c = 0; c < numSubdomainCombos; ++c)
                {
                    //TODO: each subdomain appears in many combinations. It would be faster to cache its indices for the specific dofs.
                    SignedBooleanMatrixColMajor booleanPlus = BooleanMatrices[subdomainsPlus[c].ID];
                    SignedBooleanMatrixColMajor booleanMinus = BooleanMatrices[subdomainsMinus[c].ID];

                    //TODO: The dof indices have already been accessed. Reuse it if possible.
                    IReadOnlyDictionary<IDofType, int> dofsPlus = remainderDofOrderings[subdomainsPlus[c].ID].GetDataOfRow(node);
                    IReadOnlyDictionary<IDofType, int> dofsMinus = remainderDofOrderings[subdomainsMinus[c].ID].GetDataOfRow(node);

                    foreach (IDofType dof in dofs)
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
        private List<(INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, ISubdomain[] subdomainsMinus)> 
            DefineBoundaryDofConstraints(IDofSeparator dofSeparator)
        {
            // Find boundary dual nodes and dofs
            var boundaryNodeData = new List<(INode node, IDofType[] dofs, ISubdomain[] subdomainsPlus, 
                ISubdomain[] subdomainsMinus)>(dofSeparator.GlobalBoundaryDofs.Count);

            // Find continuity equations.
            NumLagrangeMultipliers = 0;
            foreach (var nodeDofsPair in dofSeparator.GlobalBoundaryDofs)
            {
                INode node = nodeDofsPair.Key;
                IDofType[] dofsOfNode = nodeDofsPair.Value;
                ISubdomain[] nodeSubdomains = node.SubdomainsDictionary.Values.ToArray();
                ISubdomain[] subdomainsPlus, subdomainsMinus;
                int multiplicity = nodeSubdomains.Length;
                if (multiplicity == 2)
                {
                    subdomainsPlus = new ISubdomain[] { nodeSubdomains[0] };
                    subdomainsMinus = new ISubdomain[] { nodeSubdomains[1] };
                }
                else (subdomainsPlus, subdomainsMinus) = crosspointStrategy.FindSubdomainCombinations(nodeSubdomains);

                boundaryNodeData.Add((node, dofsOfNode, subdomainsPlus, subdomainsMinus));
                NumLagrangeMultipliers += dofsOfNode.Length * subdomainsPlus.Length;
            }

            return boundaryNodeData;
        }
    }
}
