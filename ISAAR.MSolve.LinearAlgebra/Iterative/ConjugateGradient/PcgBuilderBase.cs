using ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    public abstract class PcgBuilderBase
    {
        /// <summary>
        /// Specifies how to calculate the maximum iterations that the PCG algorithm will run for.
        /// </summary>
        public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);

        /// <summary>
        /// Specifies how the PCG algorithm will check that convergence has been reached.
        /// </summary>
        public IResidualConvergence ResidualConvergence { get; set; } = new SimpleConvergence();

        /// <summary>
        /// Specifies how often the residual vector will be corrected by an exact (but costly) calculation.
        /// </summary>
        public IResidualCorrection ResidualCorrection { get; set; } = new NoResidualCorrection();

        /// <summary>
        /// The PCG algorithm will converge when sqrt(r*inv(M)*r) / sqrt(r0*inv(M)*r0) &lt;= <paramref name="ResidualTolerance"/>,
        /// where M is the preconditioner, r = A*x is the current residual vector and r0 = A*x0 the initial residual vector.
        /// </summary>
        public double ResidualTolerance { get; set; } = 1E-10;
    }
}
