using System;
using System.Collections.Generic;
using System.Linq;
//using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: perhaps these helper methods should be somewhere more centrally, which will also include extracting Kib, Kii
//TODO: Reanalysis: if the global lagranges have not changed, Bpb can be reused.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning
{
    public abstract class FetiPreconditionerFactoryBase : IFetiPreconditionerFactory
    {
        public abstract bool ReorderInternalDofsForFactorization { get; }

        public abstract IFetiPreconditioner CreatePreconditioner(IStructuralModel model, 
            IStiffnessDistribution stiffnessDistribution, IDofSeparator dofSeparator, 
            ILagrangeMultipliersEnumerator lagrangeEnumerator, Dictionary<int, IFetiSubdomainMatrixManager> matrixManagers);

        protected Dictionary<int, IMappingMatrix> CalcBoundaryPreconditioningBooleanMatrices(IStiffnessDistribution stiffnessDistribution,
            IDofSeparator dofSeparator, ILagrangeMultipliersEnumerator lagrangeEnumerator)
        {
            int numContinuityEquations = lagrangeEnumerator.NumLagrangeMultipliers;
            var matricesBb = new Dictionary<int, SignedBooleanMatrixColMajor>();
            foreach (int s in dofSeparator.BoundaryDofIndices.Keys)
            {
                SignedBooleanMatrixColMajor B = lagrangeEnumerator.BooleanMatrices[s];
                SignedBooleanMatrixColMajor Bb = B.GetColumns(dofSeparator.BoundaryDofIndices[s], false);
                matricesBb[s] = Bb;
            }
            Dictionary<int, IMappingMatrix> matricesBpb = stiffnessDistribution.CalcBoundaryPreconditioningSignedBooleanMatrices(
                lagrangeEnumerator, matricesBb);
            
            return matricesBpb;
        }
    }
}
