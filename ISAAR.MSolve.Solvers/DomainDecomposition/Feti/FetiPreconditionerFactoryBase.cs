using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: perhaps these helper methods should be somewhere more centrally, which will also include extracting Kib, Kii
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public abstract class FetiPreconditionerFactoryBase : IFetiPreconditionerFactory
    {
        public abstract IFetiPreconditioner CreatePreconditioner(IStiffnessDistribution stiffnessDistribution,
            DofSeparator dofSeparator, LagrangeMultipliersEnumerator lagrangeEnumerator,
            Dictionary<int, IMatrixView> stiffnessMatrices);

        protected Dictionary<int, Matrix> CalcBoundaryPreconditioningBooleanMatrices(IStiffnessDistribution stiffnessDistribution, 
            DofSeparator dofSeparator, LagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numContinuityEquations = lagrangeEnumerator.NumLagrangeMultipliers;
            int[] rowsToKeep = Enumerable.Range(0, numContinuityEquations).ToArray(); // Same for all subdomains
            var matricesBb = new Dictionary<int, Matrix>();
            foreach (int id in dofSeparator.BoundaryDofs.Keys)
            {
                Matrix B = lagrangeEnumerator.BooleanMatrices[id].CopyToFullMatrix(false);
                Matrix Bb = B.GetSubmatrix(rowsToKeep, dofSeparator.BoundaryDofs[id]);
                matricesBb[id] = Bb;
            }
            Dictionary<int, Matrix> matricesBpb = stiffnessDistribution.CalcBoundaryPreconditioningSignedBooleanMatrices(
                dofSeparator, lagrangeEnumerator, matricesBb);
            
            return matricesBpb;
        }

        protected Dictionary<int, Matrix> ExtractBoundaryBooleanMatrices(DofSeparator dofSeparator,
            LagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numContinuityEquations = lagrangeEnumerator.NumLagrangeMultipliers;
            int[] rowsToKeep = Enumerable.Range(0, numContinuityEquations).ToArray(); // Same for all subdomains
            var boundaryBooleanMatrices = new Dictionary<int, Matrix>();
            foreach (int id in dofSeparator.BoundaryDofs.Keys)
            {
                Matrix booleanMatrix = lagrangeEnumerator.BooleanMatrices[id].CopyToFullMatrix(false);
                Matrix boundaryBooleanMatrix = booleanMatrix.GetSubmatrix(rowsToKeep, dofSeparator.BoundaryDofs[id]);
                boundaryBooleanMatrices.Add(id, boundaryBooleanMatrix);
            }
            return boundaryBooleanMatrices;
        }

        protected Dictionary<int, Matrix> ExtractStiffnessesBoundaryBoundary(DofSeparator dofSeparator,
            Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            var stiffnessesBoundaryBoundary = new Dictionary<int, Matrix>();
            foreach (int id in dofSeparator.BoundaryDofs.Keys)
            {
                int[] boundaryDofs = dofSeparator.BoundaryDofs[id];
                Matrix stiffnessBoundaryBoundary = stiffnessMatrices[id].GetSubmatrix(boundaryDofs, boundaryDofs);
                stiffnessesBoundaryBoundary.Add(id, stiffnessBoundaryBoundary);
            }
            return stiffnessesBoundaryBoundary;
        }

        protected Dictionary<int, Matrix> ExtractStiffnessBoundaryInternal(DofSeparator dofSeparator, 
            Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            var stiffnessesBoundaryInternal = new Dictionary<int, Matrix>();
            foreach (int id in dofSeparator.BoundaryDofs.Keys)
            {
                int[] boundaryDofs = dofSeparator.BoundaryDofs[id];
                int[] internalDofs = dofSeparator.InternalDofs[id];
                Matrix stiffnessBoundaryInternal = stiffnessMatrices[id].GetSubmatrix(boundaryDofs, internalDofs);
                stiffnessesBoundaryInternal.Add(id, stiffnessBoundaryInternal);
            }
            return stiffnessesBoundaryInternal;
        }
    }
}
