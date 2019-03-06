using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.DomainDecomposition.Feti;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti1
{
    internal class Feti1LumpedPreconditioner : IFetiPreconditioner
    {
        private readonly IReadOnlyList<Subdomain_v2> subdomains;
        private readonly ContinuityEquationsCalculator continuityEquations;
        private readonly Dictionary<int, int[]> boundaryDofs;

        private Dictionary<int, Matrix> boundaryBooleanMatrices;
        private Dictionary<int, Matrix> boundaryStiffnessMatrices;

        public Feti1LumpedPreconditioner(IReadOnlyList<Subdomain_v2> subdomains, ContinuityEquationsCalculator continuityEquations,
            Dictionary<int, int[]> boundaryDofs)
        {
            this.subdomains = subdomains;
            this.continuityEquations = continuityEquations;
            this.boundaryDofs = boundaryDofs;
        }

        public void CreatePreconditioner(Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            ExtractBoundaryBooleanMatrices(); //TODO: perhaps this should be done somewhere more centrally too.
            ExtractBoundaryStiffnessMatrices(stiffnessMatrices);
        }

        public void SolveLinearSystem(Vector rhs, Vector lhs)
        {
            lhs.Clear(); //TODO: this should be avoided
            DiagonalMatrix W = continuityEquations.WeightMatrix;
            Vector Wy = W.Multiply(rhs);
            foreach (Subdomain_v2 subdomain in subdomains)
            {
                int id = subdomain.ID;
                Matrix Bb = boundaryBooleanMatrices[id];
                Matrix Kbb = boundaryStiffnessMatrices[id];
                Vector contribution = W.Multiply(Bb.Multiply(Kbb.Multiply(Bb.Multiply(Wy, true))));
                lhs.AddIntoThis(contribution);
            }
        }

        private void ExtractBoundaryBooleanMatrices()
        {
            int numContinuityEquations = continuityEquations.NumContinuityEquations;
            int[] rowsToKeep = Enumerable.Range(0, numContinuityEquations).ToArray(); // Same for all subdomains
            boundaryBooleanMatrices = new Dictionary<int, Matrix>();
            foreach (Subdomain_v2 subdomain in subdomains)
            {
                int id = subdomain.ID;
                Matrix booleanMatrix = continuityEquations.BooleanMatrices[id].CopyToFullMatrix(false);
                Matrix boundaryBooleanMatrix = booleanMatrix.GetSubmatrix(rowsToKeep, boundaryDofs[id]);
                boundaryBooleanMatrices.Add(id, boundaryBooleanMatrix);
            }
        }

        private void ExtractBoundaryStiffnessMatrices(Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            boundaryStiffnessMatrices = new Dictionary<int, Matrix>();
            foreach (Subdomain_v2 subdomain in subdomains)
            {
                int id = subdomain.ID;
                IMatrixView stiffnessMatrix = stiffnessMatrices[subdomain.ID];
                Matrix boundaryStiffnessMatrix = stiffnessMatrix.GetSubmatrix(boundaryDofs[id], boundaryDofs[id]);
                boundaryStiffnessMatrices.Add(id, boundaryStiffnessMatrix);
            }
        }
    }
}
