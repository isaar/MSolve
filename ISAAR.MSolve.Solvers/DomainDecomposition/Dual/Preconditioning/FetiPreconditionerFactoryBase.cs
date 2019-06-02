using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;

//TODO: perhaps these helper methods should be somewhere more centrally, which will also include extracting Kib, Kii
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public abstract class FetiPreconditionerFactoryBase : IFetiPreconditionerFactory
    {
        public abstract IFetiPreconditioner CreatePreconditioner(IStiffnessDistribution stiffnessDistribution,
            IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator,
            Dictionary<int, IMatrixView> stiffnessMatrices);

        protected Dictionary<int, IMappingMatrix> CalcBoundaryPreconditioningBooleanMatrices(IStiffnessDistribution stiffnessDistribution,
            IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numContinuityEquations = lagrangeEnumerator.NumLagrangeMultipliers;
            var matricesBb = new Dictionary<int, SignedBooleanMatrixColMajor>();
            foreach (int id in dofSeparator.BoundaryDofIndices.Keys)
            {
                SignedBooleanMatrixColMajor B = lagrangeEnumerator.BooleanMatrices[id];
                SignedBooleanMatrixColMajor Bb = B.GetColumns(dofSeparator.BoundaryDofIndices[id], false);
                matricesBb[id] = Bb;
            }
            Dictionary<int, IMappingMatrix> matricesBpb = stiffnessDistribution.CalcBoundaryPreconditioningSignedBooleanMatrices(
                lagrangeEnumerator, matricesBb);
            
            return matricesBpb;
        }

        protected Dictionary<int, Matrix> ExtractStiffnessesBoundaryBoundary(IDofSeparator dofSeparator,
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

        protected Dictionary<int, Matrix> ExtractStiffnessBoundaryInternal(IDofSeparator dofSeparator, 
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
