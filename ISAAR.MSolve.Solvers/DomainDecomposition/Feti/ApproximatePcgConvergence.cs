using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Iterative;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: According to Fragakis PhD this is valid only for Lumped preconditioner. For other preconditioners we need to isolate
//      sum(Bpb * Kbb * Bpb^T) * residual
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    public class ApproximatePcgConvergence : IPcgConvergenceStrategy
    {
        private readonly double globalForcesNorm;
        private double residualTolerance;

        public ApproximatePcgConvergence(double globalForcesNorm)
        {
            this.globalForcesNorm = globalForcesNorm;
        }

        public bool HasConverged(IVectorView solutionVector, IVectorView preconditionedResidualVector, double residualDotProduct)
            => preconditionedResidualVector.Norm2() / globalForcesNorm <= residualTolerance;

        public void Initialize(ILinearTransformation matrix, IVectorView rhsVector, double residualTolerance, 
            double initialResidualDotProduct)
        {
            this.residualTolerance = residualTolerance; //TODO: perhaps it should be injected into the constructor.
        }
    }
}
