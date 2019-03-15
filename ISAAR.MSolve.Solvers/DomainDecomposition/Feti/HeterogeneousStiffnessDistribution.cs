using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.Commons;

//TODO: Also implement the Superlumped smoothening from Rixen, Farhat (1999). It should be equivalent to Fragakis' PhD approach.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    internal class HeterogeneousStiffnessDistribution : IStiffnessDistribution
    {
        //TODO: Use these to map displacements, forces, etc between subdomain - global level.
        private Table<INode, DOFType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses;
        private Dictionary<int, IMatrixView> stiffnessMatrices;

        public INodalLoadDistributor NodalLoadDistributor { get; } = new HeterogeneousNodalLoadDistributor();

        public Dictionary<int, Matrix> CalcBoundaryPreconditioningSignedBooleanMatrices(DofSeparator dofSeparator,
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
                    dofSeparator.BoundaryDofs[subdomainId]);
                Matrix Bb = boundarySignedBooleanMatrices[subdomainId];
                matricesBpb[subdomainId] = Bb.MultiplyRight(invDb).MultiplyLeft(Dlambda);
            }
            return matricesBpb;
        }

        //TODO: These should be done in the constructor and the fields should be readonly. However the analyzer-solver
        //      interaction needs modification, because distributing nodal loads is needed before this class can be instantiated.
        //      Since I want to use these stiffnesses for distributing nodal loads, it makes sense to make sure that this class
        //      is instantiated beforehand.
        public void StoreStiffnesses(Dictionary<int, IMatrixView> stiffnessMatrices,
            Table<INode, DOFType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses)
        {
            this.stiffnessMatrices = stiffnessMatrices;
            this.boundaryDofStiffnesses = boundaryDofStiffnesses;
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
            for (int i = 0; i < boundaryDofs.Length; ++i) inverse[i, i] = 1.0 / stiffness[i, i];
            return inverse;
        }
    }
}
