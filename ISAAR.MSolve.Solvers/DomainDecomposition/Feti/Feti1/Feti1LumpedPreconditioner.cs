using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Feti1
{
    public class Feti1LumpedPreconditioner : IFetiPreconditioner
    {
        private readonly Dictionary<int, Matrix> boundaryBooleanMatrices;
        private readonly Dictionary<int, Matrix> stiffnessesBoundaryBoundary;
        private readonly int[] subdomainIDs;
        private readonly DiagonalMatrix weightMatrix;

        private Feti1LumpedPreconditioner(int[] subdomainIDs, DiagonalMatrix weightMatrix,
            Dictionary<int, Matrix> boundaryBooleanMatrices, Dictionary<int, Matrix> stiffnessesBoundaryBoundary)
        {
            this.subdomainIDs = subdomainIDs;
            this.weightMatrix = weightMatrix;
            this.boundaryBooleanMatrices = boundaryBooleanMatrices;
            this.stiffnessesBoundaryBoundary = stiffnessesBoundaryBoundary;
        }

        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            // inv(F) * y = W * Bb * Kbb * Bb^T * W * y
            lhs.Clear(); //TODO: this should be avoided
            Vector Wy = weightMatrix.Multiply(rhs);
            foreach (int id in subdomainIDs)
            {
                Matrix Bb = boundaryBooleanMatrices[id];
                Matrix Kbb = stiffnessesBoundaryBoundary[id];
                Vector contribution = weightMatrix.Multiply(Bb.Multiply(Kbb.Multiply(Bb.Multiply(Wy, true))));
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
                return new Feti1LumpedPreconditioner(subdomainIDs, continuityEquations.WeightMatrix, boundaryBooleans,
                    stiffnessesBoundaryBoundary);
            }
        }
    }
}
