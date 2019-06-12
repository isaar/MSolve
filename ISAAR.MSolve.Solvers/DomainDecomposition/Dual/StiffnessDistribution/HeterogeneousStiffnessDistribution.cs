using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.LinearSystems;

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
        private readonly IDofSeparator dofSeparator;
        private readonly IReadOnlyList<ISubdomain> subdomains;
        private readonly Dictionary<int, DiagonalMatrix> inverseDbMatrices;

        //TODO: Is it more efficient to use (INode node, DOFType[] dofTypes)[]? It would reduce the cost of accessing node data?
        //TODO: perhaps it would be faster to have a field Dictionary<int, double[]> boundaryDofStiffnesses, instead of the next
        private Table<INode, IDofType, BoundaryDofLumpedStiffness> boundaryDofStiffnesses;

        public HeterogeneousStiffnessDistribution(IStructuralModel model, IDofSeparator dofSeparator)
        {
            this.subdomains = model.Subdomains;
            this.dofSeparator = dofSeparator;
            inverseDbMatrices = new Dictionary<int, DiagonalMatrix>();
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

        public Dictionary<int, IMappingMatrix> CalcBoundaryPreconditioningSignedBooleanMatrices(
            ILagrangeMultipliersEnumerator lagrangeEnumerator, 
            Dictionary<int, SignedBooleanMatrixColMajor> boundarySignedBooleanMatrices)
        {
            return ScalingBooleanMatrixImplicit.CreateBpbOfSubdomains(this, lagrangeEnumerator, boundarySignedBooleanMatrices);
        }

        public void Update(Dictionary<int, IMatrixView> stiffnessesFreeFree)
        {
            Table<INode, IDofType, BoundaryDofLumpedStiffness> previousBoundaryDofStiffnesses = this.boundaryDofStiffnesses;
            this.boundaryDofStiffnesses = new Table<INode, IDofType, BoundaryDofLumpedStiffness>();
            foreach (var nodeDofsPair in dofSeparator.GlobalBoundaryDofs)
            {
                INode node = nodeDofsPair.Key;
                foreach (IDofType dofType in nodeDofsPair.Value)
                {
                    var subdomainStiffnesses = new Dictionary<ISubdomain, double>();
                    double totalStiffness = 0.0;
                    foreach (ISubdomain subdomain in node.SubdomainsDictionary.Values)
                    {
                        double stiffness;
                        if (subdomain.StiffnessModified)
                        {
                            int dofIdx = subdomain.FreeDofOrdering.FreeDofs[node, dofType];
                            IMatrixView Kff = stiffnessesFreeFree[subdomain.ID];
                            stiffness = Kff[dofIdx, dofIdx]; //TODO: optimized GetDiagonal(i) method for matrices.
                        }
                        else
                        {
                            // If at the current and previous analysis steps Kff is the same for this subdomain, then Kff might 
                            // be overwritten by its factorization and we cannot access Kff(i, i)
                            stiffness = previousBoundaryDofStiffnesses[node, dofType].SubdomainStiffnesses[subdomain];
                        }
                        subdomainStiffnesses[subdomain] = stiffness;
                        totalStiffness += stiffness;
                    }
                    boundaryDofStiffnesses[node, dofType] = new BoundaryDofLumpedStiffness(subdomainStiffnesses, totalStiffness);
                }
            }
        }

        private DiagonalMatrix BuildDlambda(ILagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numLagranges = lagrangeEnumerator.NumLagrangeMultipliers;
            var Dlambda = new double[numLagranges];
            for (int i = 0; i < numLagranges; ++i)
            {
                LagrangeMultiplier lagrange = lagrangeEnumerator.LagrangeMultipliers[i];
                BoundaryDofLumpedStiffness boundaryDofStiffness = boundaryDofStiffnesses[lagrange.Node, lagrange.DofType];
                Dictionary<ISubdomain, double> stiffnessPerSubdomain = boundaryDofStiffness.SubdomainStiffnesses;
                double totalStiffness = boundaryDofStiffness.TotalStiffness;
                Dlambda[i] = stiffnessPerSubdomain[lagrange.SubdomainPlus] * stiffnessPerSubdomain[lagrange.SubdomainMinus] 
                    / totalStiffness;
            }
            return DiagonalMatrix.CreateFromArray(Dlambda, false);
        }

        //TODO: this is also done when distributing the nodal loads. Do it here only and use the inv(Db) matrix there.
        //      Even better that code should be incorporated here, and inv(Db) should be created once and stored.
        //TODO: Kbb is also calculated for most preconditioners. Just take its diagonal and invert.
        private DiagonalMatrix InvertBoundaryDofStiffnesses(ISubdomain subdomain)
        {
            if (subdomain.StiffnessModified)
            {
                Debug.WriteLine($"{this.GetType().Name}: Calculating inv(diag(Kbb)) of subomain {subdomain.ID}");
                (INode node, IDofType dofType)[] boundaryDofs = dofSeparator.BoundaryDofs[subdomain.ID];
                var invDb = new double[boundaryDofs.Length];
                for (int i = 0; i < boundaryDofs.Length; ++i)
                {
                    (INode node, IDofType dofType) = boundaryDofs[i];
                    double subdomainStiffness = boundaryDofStiffnesses[node, dofType].SubdomainStiffnesses[subdomain];
                    invDb[i] = 1.0 / subdomainStiffness;
                }
                var result = DiagonalMatrix.CreateFromArray(invDb, false);
                inverseDbMatrices[subdomain.ID] = result;
                return result;
            }
            else return inverseDbMatrices[subdomain.ID];
        }

        //TODO: This should be modified to CSR or CSC format and then benchmarked against the implicit alternative.
        /// <summary>
        /// Calculates the product Bpb = Dλ * Bb * inv(Db) explicitly, stores it and uses it for multiplications.
        /// </summary>
        private class ScalingBooleanMatrixExplicit : IMappingMatrix
        {
            private readonly Matrix explicitBpb;

            private ScalingBooleanMatrixExplicit(Matrix explicitBpb)
            {
                this.explicitBpb = explicitBpb;
            }

            public int NumColumns => explicitBpb.NumColumns;

            public int NumRows => explicitBpb.NumRows;

            internal static Dictionary<int, IMappingMatrix> CreateBpbOfSubdomains(
                HeterogeneousStiffnessDistribution stiffnessDistribution, ILagrangeMultipliersEnumerator lagrangeEnumerator, 
                Dictionary<int, SignedBooleanMatrixColMajor> boundarySignedBooleanMatrices)
            {
                // According to Fragakis PhD (e.q. 3.28): 
                // Bpb = Dλ * Bb * inv(Db(s)), Dλ[λ,λ] = K(i)[b,b] * K(j)[b,b] / Sum(K(1)[b,b] + K(2)[b,b] + ...)
                // where K(s)[b,b] is the diagonal entry of (s) subdomain's stiffess matrix corresponding to the boundary dof b 
                // and (i, j) are the subdomains connected via the Lagrange multiplier λ. 

                Matrix Dlambda = stiffnessDistribution.BuildDlambda(lagrangeEnumerator).CopyToFullMatrix(); // Common for all subdomains
                var matricesBpb = new Dictionary<int, IMappingMatrix>();
                foreach (ISubdomain subdomain in stiffnessDistribution.subdomains)
                {
                    SignedBooleanMatrixColMajor Bb = boundarySignedBooleanMatrices[subdomain.ID];
                    Matrix invDb = stiffnessDistribution.InvertBoundaryDofStiffnesses(subdomain).CopyToFullMatrix(); //TODO: This can be reused from previous analysis steps
                    Matrix Bpb = Bb.MultiplyRight(invDb).MultiplyLeft(Dlambda);
                    matricesBpb[subdomain.ID] = new ScalingBooleanMatrixExplicit(Bpb);
                }
                return matricesBpb;
            }

            public Vector Multiply(Vector vector, bool transposeThis = false)
                => explicitBpb.Multiply(vector, transposeThis);

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
                => explicitBpb.MultiplyRight(other, transposeThis);
        }

        /// <summary>
        /// Stores the matrices Bb, Dλ and inv(Db). Matrix-vector and matrix-matrix multiplications with Bpb = Dλ * Bb * inv(Db) 
        /// are performed implicitly, e.g. Bpb * x = Dλ * (Bb * (inv(Db) * x)).
        /// </summary>
        private class ScalingBooleanMatrixImplicit : IMappingMatrix
        {
            /// <summary>
            /// Signed boolean matrix with only the boundary dofs of the subdomain as columns. 
            /// </summary>
            private readonly SignedBooleanMatrixColMajor Bb;

            /// <summary>
            /// Diagonal matrix that stores for each dof the product of the stiffnesses corresponding to that dof in each 
            /// subdomain, divided by their sum.
            /// </summary>
            private readonly DiagonalMatrix Dlambda;

            /// <summary>
            /// Inverse of the diagonal matrix that stores the multiplicity of each boundary dof of the subdomain.
            /// </summary>
            private readonly DiagonalMatrix invDb;

            private ScalingBooleanMatrixImplicit(DiagonalMatrix Dlambda, SignedBooleanMatrixColMajor Bb, DiagonalMatrix invMb)
            {
                this.Dlambda = Dlambda;
                this.Bb = Bb;
                this.invDb = invMb;
            }

            public int NumColumns => invDb.NumColumns;

            public int NumRows => Dlambda.NumRows;

            internal static Dictionary<int, IMappingMatrix> CreateBpbOfSubdomains(
                HeterogeneousStiffnessDistribution stiffnessDistribution, ILagrangeMultipliersEnumerator lagrangeEnumerator,
                Dictionary<int, SignedBooleanMatrixColMajor> boundarySignedBooleanMatrices)
            {
                // According to Fragakis PhD (e.q. 3.28): 
                // Bpb = Dλ * Bb * inv(Db(s)), Dλ[λ,λ] = K(i)[b,b] * K(j)[b,b] / Sum(K(1)[b,b] + K(2)[b,b] + ...)
                // where K(s)[b,b] is the diagonal entry of (s) subdomain's stiffess matrix corresponding to the boundary dof b 
                // and (i, j) are the subdomains connected via the Lagrange multiplier λ. 

                DiagonalMatrix Dlambda = stiffnessDistribution.BuildDlambda(lagrangeEnumerator); // Common for all subdomains
                var matricesBpb = new Dictionary<int, IMappingMatrix>();
                foreach (ISubdomain subdomain in stiffnessDistribution.subdomains)
                {
                    //if (!subdomain.MaterialsModified) continue;
                    SignedBooleanMatrixColMajor Bb = boundarySignedBooleanMatrices[subdomain.ID];
                    DiagonalMatrix invDb = stiffnessDistribution.InvertBoundaryDofStiffnesses(subdomain); //TODO: This can be reused from previous analysis steps
                    matricesBpb[subdomain.ID] = new ScalingBooleanMatrixImplicit(Dlambda, Bb, invDb);
                }
                return matricesBpb;
            }

            public Vector Multiply(Vector vector, bool transposeThis = false)
            {
                //TODO: Perhaps I can reuse the temporary vectors to reduce allocations/deallocations.
                if (transposeThis)
                {
                    // Bpb^T * x = (Dλ * Bb * inv(Db))^T * x = inv(Db)^T * Bb^T * Dλ^T * x = inv(Db) * (Bb^T * (Dλ * x));
                    return invDb * Bb.Multiply(Dlambda * vector, true);
                }
                else
                {
                    // Bpb * x = Dλ * (Bb * (inv(Db) * x))
                    return Dlambda * Bb.Multiply(invDb * vector);
                }
            }

            public Matrix MultiplyRight(Matrix other, bool transposeThis = false)
            {
                //TODO: Perhaps I can reuse the temporary matrices to reduce allocations/deallocations.
                if (transposeThis)
                {
                    // Bpb^T * X = (Dλ * Bb * inv(Db))^T * X = inv(Db)^T * Bb^T * Dλ^T * X = inv(Db) * (Bb^T * (Dλ * X));
                    return invDb * Bb.MultiplyRight(Dlambda * other, true);
                }
                else
                {
                    // Bpb * X = Dλ * (Bb * (inv(Db) * X))
                    return Dlambda * Bb.MultiplyRight(invDb * other);
                }
            }
        }
    }
}
