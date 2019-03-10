using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: perhaps these helper methods should be somewhere more centrally, which will also include extracting Kib, Kii
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

        protected Dictionary<int, Matrix> ExtractStiffnessesBoundaryBoundary(Dictionary<int, int[]> boundaryDofs, 
            Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            var stiffnessesBoundaryBoundary = new Dictionary<int, Matrix>();
            foreach (int id in boundaryDofs.Keys)
            {
                Matrix stiffnessBoundaryBoundary = stiffnessMatrices[id].GetSubmatrix(boundaryDofs[id], boundaryDofs[id]);
                stiffnessesBoundaryBoundary.Add(id, stiffnessBoundaryBoundary);
            }
            return stiffnessesBoundaryBoundary;
        }

        protected Dictionary<int, Matrix> ExtractStiffnessBoundaryInternal(Dictionary<int, int[]> boundaryDofs,
            Dictionary<int, int[]> internalDofs, Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            var stiffnessesBoundaryInternal = new Dictionary<int, Matrix>();
            foreach (int id in boundaryDofs.Keys)
            {
                Matrix stiffnessBoundaryInternal = stiffnessMatrices[id].GetSubmatrix(boundaryDofs[id], internalDofs[id]);
                stiffnessesBoundaryInternal.Add(id, stiffnessBoundaryInternal);
            }
            return stiffnessesBoundaryInternal;
        }
    }
}
