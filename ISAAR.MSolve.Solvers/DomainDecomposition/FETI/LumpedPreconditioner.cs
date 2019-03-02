using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.FETI
{
    internal class LumpedPreconditioner : IFetiPreconditioner
    {
        private readonly IReadOnlyList<Subdomain_v2> subdomains;
        private readonly ContinuityEquationsCalculator continuityEquations;

        private Dictionary<int, int[]> boundaryDofs;
        private Dictionary<int, Matrix> boundaryBooleanMatrices;
        private Dictionary<int, Matrix> boundaryStiffnessMatrices;

        public LumpedPreconditioner(IReadOnlyList<Subdomain_v2> subdomains, ContinuityEquationsCalculator continuityEquations)
        {
            this.subdomains = subdomains;
            this.continuityEquations = continuityEquations;
        }

        public void CreatePreconditioner(Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            ExtractBoundaryDofs(); //TODO: this should be done somewhere more centrally.
            ExtractBoundaryBooleanMatrices(); //TODO: perhaps this too.
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

        private void ExtractBoundaryDofs()
        {
            boundaryDofs = new Dictionary<int, int[]>();
            foreach (Subdomain_v2 subdomain in subdomains)
            {
                var boundaryDofsOfSubdomain = new SortedSet<int>();
                foreach (Node_v2 node in subdomain.Nodes)
                {
                    int multiplicity = node.SubdomainsDictionary.Count;
                    if (multiplicity > 1) // boundary node
                    {
                        foreach (int dof in subdomain.FreeDofOrdering.FreeDofs.GetValuesOfRow(node))
                        {
                            boundaryDofsOfSubdomain.Add(dof);
                        }
                    }
                }
                boundaryDofs.Add(subdomain.ID, boundaryDofsOfSubdomain.ToArray());
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
