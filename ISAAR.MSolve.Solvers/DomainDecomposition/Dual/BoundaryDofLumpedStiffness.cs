using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.Commons;

//TODO: there are algebraic expressions for these. E.g. inv(Lb^T * Db * Lb) for the inverse of total stiffness. Should I use 
//      those instead?
//TODO: Could this be used outside FETI?
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual
{
    public class BoundaryDofLumpedStiffness
    {
        internal BoundaryDofLumpedStiffness(Dictionary<ISubdomain_v2, double> subdomainStiffnesses, double totalStiffness)
        {
            this.SubdomainStiffnesses = subdomainStiffnesses;
            this.TotalStiffness = totalStiffness;
        }

        internal Dictionary<ISubdomain_v2, double> SubdomainStiffnesses { get; }
        internal double TotalStiffness { get; }

        //TODO: Is it more efficient to use (INode node, DOFType[] dofTypes)[]? It would reduce the cost of accessing node data?
        public static Table<INode, DOFType, BoundaryDofLumpedStiffness> ExtractBoundaryDofLumpedStiffnesses(
            IDofSeparator dofSeparator, Dictionary<int, IMatrixView> stiffnesses)
        {
            var result = new Table<INode, DOFType, BoundaryDofLumpedStiffness>();
            foreach (var nodeDofsPair in dofSeparator.GlobalBoundaryDofs)
            {
                INode node = nodeDofsPair.Key;
                foreach (DOFType dofType in nodeDofsPair.Value)
                {
                    var subdomainStiffnesses = new Dictionary<ISubdomain_v2, double>();
                    double totalStiffness = 0.0;
                    foreach (ISubdomain_v2 subdomain in node.SubdomainsDictionary.Values)
                    {
                        int dofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                        double stiffness = stiffnesses[subdomain.ID][dofIdx, dofIdx]; //TODO: optimized GetDiagonal(i) method for matrices.
                        subdomainStiffnesses[subdomain] = stiffness;
                        totalStiffness += stiffness;
                    }
                    result[node, dofType] = new BoundaryDofLumpedStiffness(subdomainStiffnesses, totalStiffness);
                }
            }
            return result;
        }
    }
}
