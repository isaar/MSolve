using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Feti1
{
    public class Feti1DiagonalDirichletPreconditioner : IFetiPreconditioner
    {
        private readonly Dictionary<int, Matrix> boundaryBooleanMatrices;
        private readonly Dictionary<int, Matrix> stiffnessesBoundaryBoundary;
        private readonly Dictionary<int, Matrix> stiffnessesBoundaryInternal;
        private readonly Dictionary<int, DiagonalMatrix> stiffnessesInternalInternalInverseDiagonal;
        private readonly int[] subdomainIDs;
        private readonly DiagonalMatrix weightMatrix;

        private Feti1DiagonalDirichletPreconditioner(int[] subdomainIDs, DiagonalMatrix weightMatrix,
            Dictionary<int, Matrix> boundaryBooleanMatrices, Dictionary<int, Matrix> stiffnessesBoundaryBoundary,
            Dictionary<int, Matrix> stiffnessesBoundaryInternal, 
            Dictionary<int, DiagonalMatrix> stiffnessesInternalInternalInverseDiagonal)
        {
            this.subdomainIDs = subdomainIDs;
            this.weightMatrix = weightMatrix;
            this.boundaryBooleanMatrices = boundaryBooleanMatrices;
            this.stiffnessesBoundaryBoundary = stiffnessesBoundaryBoundary;
            this.stiffnessesBoundaryInternal = stiffnessesBoundaryInternal;
            this.stiffnessesInternalInternalInverseDiagonal = stiffnessesInternalInternalInverseDiagonal;
        }

        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            lhs.Clear(); //TODO: this should be avoided
            Vector Wy = weightMatrix.Multiply(rhs);
            foreach (int id in subdomainIDs)
            {
                Matrix Bb = boundaryBooleanMatrices[id];
                Matrix Kbb = stiffnessesBoundaryBoundary[id];
                Matrix Kbi = stiffnessesBoundaryInternal[id];
                DiagonalMatrix invDii = stiffnessesInternalInternalInverseDiagonal[id];

                // inv(F) * y = W * Bb * S * Bb^T * W * y
                // S = Kbb - Kbi * inv(Kii) * Kib
                Vector BWy = Bb.Multiply(Wy, true);
                Vector SBWy = Kbb.Multiply(BWy) - Kbi.Multiply(invDii.Multiply(Kbi.Multiply(BWy, true)));
                Vector contribution = weightMatrix.Multiply(Bb.Multiply(SBWy));
                lhs.AddIntoThis(contribution);
            }
        }

        public class Factory : FetiPreconditionerFactoryBase
        {
            public override IFetiPreconditioner CreatePreconditioner(Dictionary<int, int[]> boundaryDofs,
                Dictionary<int, int[]> internalDofs, ContinuityEquationsCalculator continuityEquations,
                Dictionary<int, IMatrixView> stiffnessMatrices)
            {
                int[] subdomainIDs = boundaryDofs.Keys.ToArray();
                Dictionary<int, Matrix> boundaryBooleans = ExtractBoundaryBooleanMatrices(boundaryDofs, continuityEquations);
                Dictionary<int, Matrix> stiffnessesBoundaryBoundary = 
                    ExtractStiffnessesBoundaryBoundary(boundaryDofs, stiffnessMatrices);
                Dictionary<int, Matrix> stiffnessesBoundaryInternal = 
                    ExtractStiffnessBoundaryInternal(boundaryDofs, internalDofs, stiffnessMatrices);
                Dictionary<int, DiagonalMatrix> stiffnessesInternalInternalInverseDiagonal = 
                    InvertStiffnessInternalInternalDiagonal(internalDofs, stiffnessMatrices);

                return new Feti1DiagonalDirichletPreconditioner(subdomainIDs, continuityEquations.WeightMatrix, boundaryBooleans,
                    stiffnessesBoundaryBoundary, stiffnessesBoundaryInternal, stiffnessesInternalInternalInverseDiagonal);
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
