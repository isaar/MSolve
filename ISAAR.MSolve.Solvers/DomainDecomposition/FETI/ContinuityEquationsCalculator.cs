using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.FETI
{
    /// <summary>
    /// Calculates the signed boolean matrices of the equations that enforce continuity between the multiple instances of 
    /// boundary dofs.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    internal class ContinuityEquationsCalculator
    {
        private readonly ICrosspointStrategy crosspointStrategy;

        internal ContinuityEquationsCalculator(ICrosspointStrategy crosspointStrategy)
        {
            this.crosspointStrategy = crosspointStrategy;
        }

        public Dictionary<int, SignedBooleanMatrix> CreateBooleanMatrices(Model_v2 model)
        {
            // Find boundary nodes and continuity equations
            var boundaryNodes = 
                new List<(Node_v2 node, DOFType[] dofs, Subdomain_v2[] subdomainsPlus, Subdomain_v2[] subdomainsMinus)>();
            int numContinuityEquations = 0;
            foreach (Node_v2 node in model.Nodes) //TODO: this probably doesn't work if there are embedded nodes.
            {
                int nodeMultiplicity = node.SubdomainsDictionary.Count;
                if (nodeMultiplicity > 1)
                {
                    IEnumerable<Subdomain_v2> nodeSubdomains = node.SubdomainsDictionary.Values;
                    (Subdomain_v2[] subdomainsPlus, Subdomain_v2[] subdomainsMinus) = 
                        FindSubdomainCombinations(nodeMultiplicity, nodeSubdomains);
                    DOFType[] dofsOfNode = model.GlobalDofOrdering.GlobalFreeDofs.GetColumnsOfRow(node).ToArray(); //TODO: interacting with the table can be optimized.
                    boundaryNodes.Add((node, dofsOfNode, subdomainsPlus, subdomainsMinus));
                    numContinuityEquations += dofsOfNode.Length * subdomainsPlus.Length;
                }
            }

            // Create the signed boolean matrices.
            var booleanMatrices = new Dictionary<int, SignedBooleanMatrix>();
            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                booleanMatrices.Add(subdomain.ID, 
                    new SignedBooleanMatrix(numContinuityEquations, subdomain.FreeDofOrdering.NumFreeDofs));
            }

            // Fill the boolean matrices: node major, subdomain medium, dof minor. TODO: not sure about this order.
            int equationCounter = 0;
            foreach ((Node_v2 node, DOFType[] dofs, Subdomain_v2[] subdomainsPlus, Subdomain_v2[] subdomainsMinus) 
                in boundaryNodes)
            {
                int numSubdomainCombos = subdomainsPlus.Length;
                for (int c = 0; c < numSubdomainCombos; ++c)
                {
                    //TODO: each subdomain appears in many combinations. It would be faster to cache its indices for the specific dofs.
                    SignedBooleanMatrix booleanPlus = booleanMatrices[subdomainsPlus[c].ID];
                    SignedBooleanMatrix booleanMinus = booleanMatrices[subdomainsMinus[c].ID];
                    bool nodeIsOrdered = subdomainsPlus[c].FreeDofOrdering.FreeDofs.TryGetDataOfRow(node,
                        out IReadOnlyDictionary<DOFType, int> dofsPlus);
                    Debug.Assert(nodeIsOrdered);
                    nodeIsOrdered = subdomainsMinus[c].FreeDofOrdering.FreeDofs.TryGetDataOfRow(node,
                        out IReadOnlyDictionary<DOFType, int> dofsMinus);
                    Debug.Assert(nodeIsOrdered);

                    foreach (DOFType dof in dofs)
                    {
                        booleanPlus.AddEntry(equationCounter, dofsPlus[dof], true);
                        booleanMinus.AddEntry(equationCounter, dofsMinus[dof], false);
                        ++equationCounter;
                    }
                }
            }

            return booleanMatrices;
        }

        //TODO: Perhaps a list of pairs is better than a pair of lists.
        private (Subdomain_v2[] subdomainsPlus, Subdomain_v2[] subdomainsMinus) FindSubdomainCombinations(int nodeMultiplicity,
            IEnumerable<Subdomain_v2> nodeSubdomains)
        {
            Debug.Assert(nodeMultiplicity > 1);
            int numNodeCombos = (nodeMultiplicity * (nodeMultiplicity - 1)) / 2; //TODO: not sure about this
            var subdomainsPlus = new Subdomain_v2[numNodeCombos];
            var subdomainsMinus = new Subdomain_v2[numNodeCombos];

            var processedSubdomains = new HashSet<Subdomain_v2>(nodeSubdomains);
            int comboCounter = 0;
            foreach (Subdomain_v2 subdomain1 in nodeSubdomains)
            {
                processedSubdomains.Remove(subdomain1);
                foreach (Subdomain_v2 subdomain2 in processedSubdomains)
                {
                    subdomainsPlus[comboCounter] = subdomain1;
                    subdomainsMinus[comboCounter] = subdomain2;
                    ++comboCounter;
                }
            }
            return (subdomainsPlus, subdomainsMinus);
        }

    }
}
