using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.Commons;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;

//TODO: Also implement the Superlumped smoothening from Rixen, Farhat (1999). It should be equivalent to Fragakis' PhD approach.
//TODO: Use this to map displacements, forces, etc between subdomain - global level.
//TODO: perhaps I should store the stiffness of each boundary dof per subdomain, instead of storing the same data
//      Table<INode, DOFType, BoundaryDofLumpedStiffness>.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution
{
    public class HeterogeneousStiffnessDistribution : IFeti1StiffnessDistribution
    {
        //TODO: If a solver operation needs this, it is probably better to delegate that operation to this class.
        private readonly Table<INode, DOFType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses;
        private readonly Feti1DofSeparator dofSeparator;
        private readonly IStructuralModel_v2 model;
        private readonly Dictionary<int, IMatrixView> stiffnessMatrices;

        public HeterogeneousStiffnessDistribution(IStructuralModel_v2 model, Feti1DofSeparator dofSeparator, 
            Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            this.model = model;
            this.dofSeparator = dofSeparator;
            this.stiffnessMatrices = stiffnessMatrices;
            this.boundaryDofStiffnesses = BoundaryDofLumpedStiffness.ExtractBoundaryDofLumpedStiffnesses(
                dofSeparator, stiffnessMatrices);
            this.SubdomainGlobalConversion = new HeterogeneousSubdomainGlobalConversion(model, dofSeparator, 
                FindRelativeBoundaryStiffnesses(model, dofSeparator, boundaryDofStiffnesses));
        }

        public ISubdomainGlobalConversion SubdomainGlobalConversion { get; }

        public Dictionary<int, Matrix> CalcBoundaryPreconditioningSignedBooleanMatrices(
            LagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, Matrix> boundarySignedBooleanMatrices)
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

        private static Dictionary<int, double[]> FindRelativeBoundaryStiffnesses(IStructuralModel_v2 model,
            Feti1DofSeparator dofSeparator, Table<INode, DOFType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses)
        {
            //TODO: this should probably be handled when extracting the lumped boundary stiffnesses, in order to avoid searching
            //      the subdomain in BoundaryDofLumpedStiffness.SubdomainStiffnesses for each dof.
            var relativeBoundaryStiffnesses = new Dictionary<int, double[]>();
            foreach (ISubdomain_v2 subdomain in model.Subdomains)
            {
                int[] boundaryDofIndices = dofSeparator.BoundaryDofIndices[subdomain.ID];
                (INode, DOFType)[] boundaryDofConnectivities = dofSeparator.BoundaryDofConnectivities[subdomain.ID];
                int numBoundaryDofs = boundaryDofIndices.Length;
                var subdomainStiffnesses = new double[numBoundaryDofs];
                for (int i = 0; i < boundaryDofIndices.Length; ++i)
                {
                    (INode node, DOFType dofType) = boundaryDofConnectivities[i];
                    BoundaryDofLumpedStiffness dofStiffness = boundaryDofStiffnesses[node, dofType];
                    double relativeStiffness = dofStiffness.SubdomainStiffnesses[subdomain] / dofStiffness.TotalStiffness;
                    subdomainStiffnesses[i] = relativeStiffness;
                }
                relativeBoundaryStiffnesses.Add(subdomain.ID, subdomainStiffnesses);
            }
            return relativeBoundaryStiffnesses;
        }

        private Matrix BuildDlambda(LagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numLagranges = lagrangeEnumerator.NumLagrangeMultipliers;
            var Dlambda = Matrix.CreateZero(numLagranges, numLagranges);
            for (int i = 0; i < numLagranges; ++i)
            {
                LagrangeMultiplier lagrange = lagrangeEnumerator.LagrangeMultipliers[i];
                BoundaryDofLumpedStiffness boundaryDofStiffness = boundaryDofStiffnesses[lagrange.Node, lagrange.DofType];
                Dictionary<ISubdomain_v2, double> stiffnessPerSubdomain = boundaryDofStiffness.SubdomainStiffnesses;
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
