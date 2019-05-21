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
//TODO: In FETI-DP, this class should operate using Krr, while it now uses Kff. The problem is that Krr is created after the 
//      global loads are distributed, which is when the IStiffnessDistribution is first created.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution
{
    public class HeterogeneousStiffnessDistribution : IStiffnessDistribution
    {
        //TODO: perhaps it would be faster to have a field Dictionary<int, double[]> boundaryDofStiffnesses, instead of the next
        private readonly Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses;

        private readonly IDofSeparator dofSeparator;
        private readonly IStructuralModel model;

        public HeterogeneousStiffnessDistribution(IStructuralModel model, IDofSeparator dofSeparator,
            Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.boundaryDofStiffnesses = boundaryDofStiffnesses;
        }

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
                BoundaryDofLumpedStiffness dofStiffness = boundaryDofStiffnesses[node, dofType];
                double relativeStiffness = dofStiffness.SubdomainStiffnesses[subdomain] / dofStiffness.TotalStiffness;
                relativeStiffnesses[i] = relativeStiffness;
            }
            return relativeStiffnesses;
        }

        public Dictionary<int, double> CalcBoundaryDofCoefficients(INode node, IDofType dofType)
        {
            var coeffs = new Dictionary<int, double>();
            BoundaryDofLumpedStiffness dofStiffness = boundaryDofStiffnesses[node, dofType];
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
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                Matrix invDb = InvertBoundaryDofStiffnesses(subdomain);
                Matrix Bb = boundarySignedBooleanMatrices[subdomain.ID];
                matricesBpb[subdomain.ID] = Bb.MultiplyRight(invDb).MultiplyLeft(Dlambda);
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
                BoundaryDofLumpedStiffness boundaryDofStiffness = boundaryDofStiffnesses[lagrange.Node, lagrange.DofType];
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
        private Matrix InvertBoundaryDofStiffnesses(ISubdomain subdomain)
        {
            (INode node, IDofType dofType)[] boundaryDofs = dofSeparator.BoundaryDofs[subdomain.ID];
            Matrix Db = Matrix.CreateZero(boundaryDofs.Length, boundaryDofs.Length);
            for (int i = 0; i < boundaryDofs.Length; ++i)
            {
                (INode node, IDofType dofType) = boundaryDofs[i];
                double subdomainStiffness = boundaryDofStiffnesses[node, dofType].SubdomainStiffnesses[subdomain];
                Db[i, i] = 1.0 / subdomainStiffness;
            }
            return Db;
        }

    }
}
