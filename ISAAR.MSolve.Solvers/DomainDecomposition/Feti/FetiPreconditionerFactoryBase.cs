using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: perhaps these helper methods should be somewhere more centrally.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public abstract class FetiPreconditionerFactoryBase : IFetiPreconditionerFactory
    {
        public abstract IFetiPreconditioner CreatePreconditioner(Dictionary<int, int[]> boundaryDofs,
            Dictionary<int, int[]> internalDofs, ContinuityEquationsCalculator continuityEquations,
            Dictionary<int, IMatrixView> stiffnessMatrices);

        protected Dictionary<int, Matrix> ExtractBoundaryBooleanMatrices(Dictionary<int, int[]> boundaryDofs, 
            ContinuityEquationsCalculator continuityEquations)
        {
            int numContinuityEquations = continuityEquations.NumContinuityEquations;
            int[] rowsToKeep = Enumerable.Range(0, numContinuityEquations).ToArray(); // Same for all subdomains
            var boundaryBooleanMatrices = new Dictionary<int, Matrix>();
            foreach (int id in boundaryDofs.Keys)
            {
                Matrix booleanMatrix = continuityEquations.BooleanMatrices[id].CopyToFullMatrix(false);
                Matrix boundaryBooleanMatrix = booleanMatrix.GetSubmatrix(rowsToKeep, boundaryDofs[id]);
                boundaryBooleanMatrices.Add(id, boundaryBooleanMatrix);
            }
            return boundaryBooleanMatrices;
        }

        protected Dictionary<int, Matrix> ExtractBoundaryStiffnessMatrices(Dictionary<int, int[]> boundaryDofs, 
            Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            var boundaryStiffnessMatrices = new Dictionary<int, Matrix>();
            foreach (int id in boundaryDofs.Keys)
            {
                IMatrixView stiffnessMatrix = stiffnessMatrices[id];
                Matrix boundaryStiffnessMatrix = stiffnessMatrix.GetSubmatrix(boundaryDofs[id], boundaryDofs[id]);
                boundaryStiffnessMatrices.Add(id, boundaryStiffnessMatrix);
            }
            return boundaryStiffnessMatrices;
        }
    }
}
