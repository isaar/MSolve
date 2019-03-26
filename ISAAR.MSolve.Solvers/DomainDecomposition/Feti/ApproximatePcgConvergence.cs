using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;

//TODO: According to Fragakis PhD this is valid only for Lumped preconditioner. For other preconditioners we need to isolate
//      sum(Bpb * Kbb * Bpb^T) * residual
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public class ApproximatePcgConvergence : IPcgResidualConvergence
    {
        private readonly double globalForcesNorm;

        public ApproximatePcgConvergence(double globalForcesNorm)
        {
            this.globalForcesNorm = globalForcesNorm;
        }

        public double EstimateResidualNormRatio(PcgAlgorithmBase pcg) 
            => pcg.PrecondResidual.Norm2() / globalForcesNorm;

        public void Initialize(PcgAlgorithmBase pcg) { } // Do nothing 
    }
}
