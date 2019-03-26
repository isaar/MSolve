using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: According to Fragakis PhD this is valid only for Lumped preconditioner. For other preconditioners we need to isolate
//      sum(Bpb * Kbb * Bpb^T) * residual
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti.Pcpg
{
    public class ApproximateResidualConvergence : IPcgResidualConvergence, IPcpgResidualConvergence
    {
        public ApproximateResidualConvergence(double globalForces)
        {
            this.GlobalForcesNorm = globalForces;
        }

        //TODO: this should be injected in when FETI initializes the object, thus the user needs another way to specify this strategy
        internal double GlobalForcesNorm { get; } = 0.0;

        public double EstimateResidualNormRatio(PcgAlgorithmBase pcg)
        {
            Debug.Assert(GlobalForcesNorm != 0.0, "norm2(globalForces) must be set first");
            return pcg.PrecondResidual.Norm2() / GlobalForcesNorm;
        }

        public double EstimateResidualNormRatio(IVectorView lagrangeMultipliers, IVectorView projectedPrecondResidual)
        {
            Debug.Assert(GlobalForcesNorm != 0.0, "norm2(globalForces) must be set first");
            return projectedPrecondResidual.Norm2() / GlobalForcesNorm;
        }

        public void Initialize(PcgAlgorithmBase pcg) { } // Do nothing
    }
}
