using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.ResidualUpdate;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Base abstract class for Preconditioned Conjugate Gradient implementations.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public abstract class PcgAlgorithmBase
    {
        protected readonly IMaxIterationsProvider maxIterationsProvider;
        protected readonly double residualTolerance;
        protected readonly IResidualConvergence residualConvergence;
        protected readonly IResidualCorrection residualCorrection;

        protected PcgAlgorithmBase(double residualTolerance, IMaxIterationsProvider maxIterationsProvider, 
            IResidualConvergence residualConvergence, IResidualCorrection residualCorrection)
        {
            this.residualTolerance = residualTolerance;
            this.maxIterationsProvider = maxIterationsProvider;
            this.residualConvergence = residualConvergence;
            this.residualCorrection = residualCorrection;
        }

        /// <summary>
        /// Solves the linear system A * x = b by solving the preconditioned system inv(P) * A * inv(P)^T * y = inv(P) * b, 
        /// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/>, x is the solution, y = P^T * x,
        /// P*P^T = <paramref name="preconditioner"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
        /// <param name="rhs">
        /// The right hand side vector b of the linear system A * x = b. Constraints:
        /// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="preconditioner">
        /// A preconditioner matrix that is also symmetric positive definite and has the same dimensions as A.
        /// </param>
        /// <param name="solution">
        /// The vector from which to start refining the solution vector x. Constraints:
        /// <paramref name="solution"/>.<see cref="IIndexable1D.Length"/>
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <param name="initialGuessIsZero">
        /// If <paramref name="solution"/> is 0, then set <paramref name="initialGuessIsZero"/> to true to avoid performing the
        /// operation b-A*0 before starting.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="rhs"/> or <paramref name="solution"/> violate the described constraints.
        /// </exception>
        public CGStatistics Solve(IMatrixView matrix, IPreconditioner preconditioner, IVectorView rhs, IVector solution,
            bool initialGuessIsZero, Func<IVector> zeroVectorInitializer) //TODO: find a better way to handle the case x0=0
        {
            return Solve(new ExplicitMatrixTransformation(matrix), preconditioner, rhs, solution, initialGuessIsZero, 
                zeroVectorInitializer);
        }

        /// <summary>
        /// Solves the linear system A * x = b by solving the preconditioned system inv(P) * A * inv(P)^T * y = inv(P) * b, 
        /// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/>, x is the solution, y = P^T * x,
        /// P*P^T = <paramref name="preconditioner"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">
        /// Represents the matrix A of the linear system A * x = b, which must be symmetric positive definite.
        /// </param>
        /// <param name="rhs">
        /// The right hand side vector b of the linear system A * x = b. Constraints:
        /// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <param name="preconditioner">
        /// A preconditioner matrix that is also symmetric positive definite and has the same dimensions as A.
        /// </param>
        /// <param name="solution">
        /// The vector from which to start refining the solution vector x. Constraints:
        /// <paramref name="solution"/>.<see cref="IIndexable1D.Length"/>
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <param name="initialGuessIsZero">
        /// If <paramref name="solution"/> is 0, then set <paramref name="initialGuessIsZero"/> to true to avoid performing the
        /// operation b-A*0 before starting.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="rhs"/> or <paramref name="solution"/> violate the described constraints.
        /// </exception>
        public CGStatistics Solve(ILinearTransformation matrix, IPreconditioner preconditioner, IVectorView rhs,
            IVector solution, bool initialGuessIsZero, Func<IVector> zeroVectorInitializer) //TODO: find a better way to handle the case x0=0
        {
            Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
            Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

            // r = b - A * x
            IVector res;
            if (initialGuessIsZero) res = rhs.Copy();
            else res = ExactResidual.Calculate(matrix, rhs, solution);
            return SolveInternal(matrix, preconditioner, rhs, solution, res, zeroVectorInitializer);
        }

        protected abstract CGStatistics SolveInternal(ILinearTransformation matrix, IPreconditioner preconditioner,
            IVectorView rhs, IVector solution, IVector residual, Func<IVector> zeroVectorInitializer);
    }
}
