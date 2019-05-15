using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: Also implement the Superlumped smoothening from Rixen, Farhat (1999). It should be equivalent to Fragakis' PhD approach.
//TODO: Use this to map displacements, forces, etc between subdomain - global level.
//TODO: perhaps I should store the stiffness of each boundary dof per subdomain, instead of storing the same data
//      Table<INode, DOFType, BoundaryDofLumpedStiffness>.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public class HeterogeneousStiffnessDistribution : IStiffnessDistribution
    {
        private readonly IDofSeparator dofSeparator;
        private readonly IStructuralModel model;
        private readonly Dictionary<int, IMatrixView> stiffnessMatrices;

        public HeterogeneousStiffnessDistribution(IStructuralModel model, IDofSeparator dofSeparator, 
            Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.stiffnessMatrices = stiffnessMatrices;

            this.BoundaryDofStiffnesses = BoundaryDofLumpedStiffness.ExtractBoundaryDofLumpedStiffnesses(
                dofSeparator.GlobalBoundaryDofs, stiffnessMatrices);
        }

        //TODO: If a solver operation needs this, it is probably better to delegate that operation to this class.
        public Table<INode, IDofType, BoundaryDofLumpedStiffness> BoundaryDofStiffnesses { get; }

        public double[] CalcBoundaryDofCoefficients(ISubdomain subdomain)
        {
            //TODO: Should this be cached? It stores the same info as HeterogeneousStiffnessDistribution.BoundaryDofStiffnesses.
            //      This format is more compact and has more efficient indexing when dealing with only 1 subdomain at a time, 
            //      but is difficult to index when accessing global boundary dofs. Is it possible to only use one of the two 
            //      formats? 
            //TODO: Could this be handled when extracting the lumped boundary stiffnesses? That way we can avoid searching
            //      the subdomain in BoundaryDofLumpedStiffness.SubdomainStiffnesses for each dof.

            int[] boundaryDofIndices = dofSeparator.BoundaryDofIndices[subdomain.ID];
            (INode, IDofType)[] boundaryDofs = dofSeparator.BoundaryDofs[subdomain.ID];
            int numBoundaryDofs = boundaryDofIndices.Length;
            var relativeStiffnesses = new double[numBoundaryDofs];
            for (int i = 0; i < boundaryDofIndices.Length; ++i)
            {
                (INode node, IDofType dofType) = boundaryDofs[i];
                BoundaryDofLumpedStiffness dofStiffness = BoundaryDofStiffnesses[node, dofType];
                double relativeStiffness = dofStiffness.SubdomainStiffnesses[subdomain] / dofStiffness.TotalStiffness;
                relativeStiffnesses[i] = relativeStiffness;
            }
            return relativeStiffnesses;
        }

        public Dictionary<int, double> CalcBoundaryDofCoefficients(INode node, IDofType dofType)
        {
            var coeffs = new Dictionary<int, double>();
            BoundaryDofLumpedStiffness dofStiffness = BoundaryDofStiffnesses[node, dofType];
            foreach (var idSubdomainPair in node.SubdomainsDictionary)
            {
                int id = idSubdomainPair.Key;
                ISubdomain subdomain = idSubdomainPair.Value;
                coeffs[id] = dofStiffness.SubdomainStiffnesses[subdomain] / dofStiffness.TotalStiffness;
            }
            return coeffs;
        }

        public Dictionary<int, Matrix> CalcBoundaryPreconditioningSignedBooleanMatrices(
            ILagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, Matrix> boundarySignedBooleanMatrices)
        {
            // According to Fragakis PhD (e.q. 3.28): 
            // Bpb = Dλ * Bb * inv(Db(s)), Dλ[λ,λ] = K(i)[b,b] * K(j)[b,b] / Sum(K(1)[b,b] + K(2)[b,b] + ...)
            // where K(s)[b,b] is the diagonal entry of (s) subdomain's stiffess matrix corresponding to the boundary dof b 
            // and (i, j) are the subdomains connected via the Lagrange multiplier λ. 

            Matrix Dlambda = BuildDlambda(lagrangeEnumerator); // Common for all subdomains
            var matricesBpb = new Dictionary<int, Matrix>();
            foreach (int subdomainId in boundarySignedBooleanMatrices.Keys)
            {
                Matrix invDb = InvertBoundaryDofStiffnesses(stiffnessMatrices[subdomainId],
                    dofSeparator.BoundaryDofIndices[subdomainId]);
                Matrix Bb = boundarySignedBooleanMatrices[subdomainId];
                matricesBpb[subdomainId] = Bb.MultiplyRight(invDb).MultiplyLeft(Dlambda);
            }
            return matricesBpb;
        }

        private Matrix BuildDlambda(ILagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numLagranges = lagrangeEnumerator.NumLagrangeMultipliers;
            var Dlambda = Matrix.CreateZero(numLagranges, numLagranges);
            for (int i = 0; i < numLagranges; ++i)
            {
                LagrangeMultiplier lagrange = lagrangeEnumerator.LagrangeMultipliers[i];
                BoundaryDofLumpedStiffness boundaryDofStiffness = BoundaryDofStiffnesses[lagrange.Node, lagrange.DofType];
                Dictionary<ISubdomain, double> stiffnessPerSubdomain = boundaryDofStiffness.SubdomainStiffnesses;
                double totalStiffness = boundaryDofStiffness.TotalStiffness;
                Dlambda[i, i] = stiffnessPerSubdomain[lagrange.SubdomainPlus] * stiffnessPerSubdomain[lagrange.SubdomainMinus] 
                    / totalStiffness;
            }
            return Dlambda;
        }

        //TODO: this is also done when distributing the nodal loads. Do it here only and use the inv(Db) matrix there.
        //      Even better that code should be incorporated here, and inv(Db) should be created once and stored.
        //TODO: Kbb is also calculated for most preconditioners. Just take its diagonal and invert.
        //TODO: Can this be done using the class BoundaryDofLumpedStiffness? Should that class manage the Db(s) matrices?
        private Matrix InvertBoundaryDofStiffnesses(IMatrixView stiffness, int[] boundaryDofs)
        {
            var inverse = Matrix.CreateZero(boundaryDofs.Length, boundaryDofs.Length);
            for (int i = 0; i < boundaryDofs.Length; ++i) inverse[i, i] = 1.0 / stiffness[boundaryDofs[i], boundaryDofs[i]];
            return inverse;
        }
    }
}
