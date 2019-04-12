using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg;

//TODO: perhaps these helper methods should be somewhere more centrally, which will also include extracting Kib, Kii
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Preconditioning
{
    public abstract class Feti1PreconditionerFactoryBase : IFetiPreconditionerFactory
    {
        public abstract IFetiPreconditioner CreatePreconditioner(IFeti1StiffnessDistribution stiffnessDistribution,
            Feti1DofSeparator dofSeparator, LagrangeMultipliersEnumerator lagrangeEnumerator,
            Dictionary<int, IMatrixView> stiffnessMatrices);

        protected Dictionary<int, Matrix> CalcBoundaryPreconditioningBooleanMatrices(IFeti1StiffnessDistribution stiffnessDistribution, 
            Feti1DofSeparator dofSeparator, LagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numContinuityEquations = lagrangeEnumerator.NumLagrangeMultipliers;
            int[] rowsToKeep = Enumerable.Range(0, numContinuityEquations).ToArray(); // Same for all subdomains
            var matricesBb = new Dictionary<int, Matrix>();
            foreach (int id in dofSeparator.BoundaryDofIndices.Keys)
            {
                Matrix B = lagrangeEnumerator.BooleanMatrices[id].CopyToFullMatrix(false);
                Matrix Bb = B.GetSubmatrix(rowsToKeep, dofSeparator.BoundaryDofIndices[id]);
                matricesBb[id] = Bb;
            }
            Dictionary<int, Matrix> matricesBpb = stiffnessDistribution.CalcBoundaryPreconditioningSignedBooleanMatrices(
                lagrangeEnumerator, matricesBb);
            
            return matricesBpb;
        }

        protected Dictionary<int, Matrix> ExtractBoundaryBooleanMatrices(Feti1DofSeparator dofSeparator,
            LagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numContinuityEquations = lagrangeEnumerator.NumLagrangeMultipliers;
            int[] rowsToKeep = Enumerable.Range(0, numContinuityEquations).ToArray(); // Same for all subdomains
            var boundaryBooleanMatrices = new Dictionary<int, Matrix>();
            foreach (int id in dofSeparator.BoundaryDofIndices.Keys)
            {
                Matrix booleanMatrix = lagrangeEnumerator.BooleanMatrices[id].CopyToFullMatrix(false);
                Matrix boundaryBooleanMatrix = booleanMatrix.GetSubmatrix(rowsToKeep, dofSeparator.BoundaryDofIndices[id]);
                boundaryBooleanMatrices.Add(id, boundaryBooleanMatrix);
            }
            return boundaryBooleanMatrices;
        }

        protected Dictionary<int, Matrix> ExtractStiffnessesBoundaryBoundary(Feti1DofSeparator dofSeparator,
            Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            var stiffnessesBoundaryBoundary = new Dictionary<int, Matrix>();
            foreach (int id in dofSeparator.BoundaryDofIndices.Keys)
            {
                int[] boundaryDofs = dofSeparator.BoundaryDofIndices[id];
                Matrix stiffnessBoundaryBoundary = stiffnessMatrices[id].GetSubmatrix(boundaryDofs, boundaryDofs);
                stiffnessesBoundaryBoundary.Add(id, stiffnessBoundaryBoundary);
            }
            return stiffnessesBoundaryBoundary;
        }

        protected Dictionary<int, Matrix> ExtractStiffnessBoundaryInternal(Feti1DofSeparator dofSeparator, 
            Dictionary<int, IMatrixView> stiffnessMatrices)
        {
            var stiffnessesBoundaryInternal = new Dictionary<int, Matrix>();
            foreach (int id in dofSeparator.BoundaryDofIndices.Keys)
            {
                int[] boundaryDofs = dofSeparator.BoundaryDofIndices[id];
                int[] internalDofs = dofSeparator.InternalDofIndices[id];
                Matrix stiffnessBoundaryInternal = stiffnessMatrices[id].GetSubmatrix(boundaryDofs, internalDofs);
                stiffnessesBoundaryInternal.Add(id, stiffnessBoundaryInternal);
            }
            return stiffnessesBoundaryInternal;
        }
    }
}
