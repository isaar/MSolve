using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Preconditioning
{
    public class Feti1DiagonalDirichletPreconditioner : IFetiPreconditioner
    {
        private readonly Dictionary<int, Matrix> preconditioningBoundarySignedBooleanMatrices;
        private readonly Dictionary<int, Matrix> stiffnessesBoundaryBoundary;
        private readonly Dictionary<int, Matrix> stiffnessesBoundaryInternal;
        private readonly Dictionary<int, DiagonalMatrix> stiffnessesInternalInternalInverseDiagonal;
        private readonly int[] subdomainIDs;

        private Feti1DiagonalDirichletPreconditioner(int[] subdomainIDs, Dictionary<int, Matrix> stiffnessesBoundaryBoundary,
            Dictionary<int, Matrix> stiffnessesBoundaryInternal, 
            Dictionary<int, DiagonalMatrix> stiffnessesInternalInternalInverseDiagonal,
            Dictionary<int, Matrix> preconditioningBoundarySignedBooleanMatrices)
        {
            this.subdomainIDs = subdomainIDs;
            this.preconditioningBoundarySignedBooleanMatrices = preconditioningBoundarySignedBooleanMatrices;
            this.stiffnessesBoundaryBoundary = stiffnessesBoundaryBoundary;
            this.stiffnessesBoundaryInternal = stiffnessesBoundaryInternal;
            this.stiffnessesInternalInternalInverseDiagonal = stiffnessesInternalInternalInverseDiagonal;
        }

        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            lhs.Clear(); //TODO: this should be avoided
            foreach (int id in subdomainIDs)
            {
                Matrix Bpb = preconditioningBoundarySignedBooleanMatrices[id];
                Matrix Kbb = stiffnessesBoundaryBoundary[id];
                Matrix Kbi = stiffnessesBoundaryInternal[id];
                DiagonalMatrix invDii = stiffnessesInternalInternalInverseDiagonal[id];

                // inv(F) * y = Bpb * S * Bpb^T * y
                // S = Kbb - Kbi * inv(Dii) * Kib
                Vector By = Bpb.Multiply(rhs, true);
                Vector SBy = Kbb.Multiply(By) - Kbi.Multiply(invDii.Multiply(Kbi.Multiply(By, true)));
                Vector subdomainContribution = Bpb.Multiply(SBy);
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public void SolveLinearSystems(Matrix rhs, Matrix lhs)
        {
            lhs.Clear(); //TODO: this should be avoided
            foreach (int id in subdomainIDs)
            {
                Matrix Bpb = preconditioningBoundarySignedBooleanMatrices[id];
                Matrix Kbb = stiffnessesBoundaryBoundary[id];
                Matrix Kbi = stiffnessesBoundaryInternal[id];
                DiagonalMatrix invDii = stiffnessesInternalInternalInverseDiagonal[id];

                // inv(F) * Y = Bpb * S * Bpb^T * Y
                // S = Kbb - Kbi * inv(Dii) * Kib
                Matrix BY = Bpb.MultiplyRight(rhs, true);
                Matrix SBY = Kbb.MultiplyRight(BY) - Kbi.MultiplyRight(invDii.MultiplyRight(Kbi.MultiplyRight(BY, true)));
                Matrix subdomainContribution = Bpb.MultiplyRight(SBY);
                lhs.AddIntoThis(subdomainContribution);
            }
        }

        public class Factory : Feti1PreconditionerFactoryBase
        {
            public override IFetiPreconditioner CreatePreconditioner(IFeti1StiffnessDistribution stiffnessDistribution,
                Feti1DofSeparator dofSeparator, LagrangeMultipliersEnumerator lagrangeEnumerator,
                Dictionary<int, IMatrixView> stiffnessMatrices)
            {
                int[] subdomainIDs = dofSeparator.BoundaryDofIndices.Keys.ToArray();
                Dictionary<int, Matrix> boundaryBooleans = CalcBoundaryPreconditioningBooleanMatrices(stiffnessDistribution,
                    dofSeparator, lagrangeEnumerator);
                Dictionary<int, Matrix> stiffnessesBoundaryBoundary = 
                    ExtractStiffnessesBoundaryBoundary(dofSeparator, stiffnessMatrices);
                Dictionary<int, Matrix> stiffnessesBoundaryInternal = 
                    ExtractStiffnessBoundaryInternal(dofSeparator, stiffnessMatrices);
                Dictionary<int, DiagonalMatrix> stiffnessesInternalInternalInverseDiagonal = 
                    InvertStiffnessInternalInternalDiagonal(dofSeparator.InternalDofIndices, stiffnessMatrices);

                return new Feti1DiagonalDirichletPreconditioner(subdomainIDs, stiffnessesBoundaryBoundary, 
                    stiffnessesBoundaryInternal, stiffnessesInternalInternalInverseDiagonal, boundaryBooleans);
            }

            private Dictionary<int, DiagonalMatrix> InvertStiffnessInternalInternalDiagonal(Dictionary<int, int[]> internalDofs, 
                Dictionary<int, IMatrixView> stiffnessMatrices)
            {
                var stiffnessesInternalInternalInverse = new Dictionary<int, DiagonalMatrix>();
                foreach (int id in internalDofs.Keys)
                {
                    var diagonal = new double[internalDofs[id].Length];
                    for (int i = 0; i < diagonal.Length; ++i)
                    {
                        int idx = internalDofs[id][i];
                        diagonal[i] = stiffnessMatrices[id][idx, idx];
                    }
                    var matrix = DiagonalMatrix.CreateFromArray(diagonal, false);
                    matrix.Invert();
                    stiffnessesInternalInternalInverse.Add(id, matrix);
                }
                return stiffnessesInternalInternalInverse;
            }
        }
    }
}
