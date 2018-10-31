using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Needs Builder pattern
//TODO: initialization should be done by a vector factory, instead of new Vector(..)
namespace ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG
{
    /// <summary>
    /// Implements the Conjugate Gradient algorithm with preconditioning for solving linear systems with symmetric positive 
    /// definite matrices.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PreconditionedConjugateGradient
    {
        private readonly int maxIterations;
        private readonly double residualTolerance;

        /// <summary>
        /// Initializes a new instance of <see cref="PreconditionedConjugateGradient"/> with the specified convergence criteria. 
        /// If any criterion is met the algorithm will temrinate.
        /// </summary>
        /// <param name="maxIterations">The maximum number of iterations before the algorithm terminates.</param>
        /// <param name="residualTolerance">If norm2(b-A*x) / norm2(b-A*x0) &lt;= <paramref name="residualTolerance"/>, the 
        ///     algorithm will terminate, where x is the current solution vector and x0 the initial guess.</param>
        public PreconditionedConjugateGradient(int maxIterations, double residualTolerance)
        {
            this.maxIterations = maxIterations;
            this.residualTolerance = residualTolerance;
        }

        /// <summary>
        /// Solves the linear system A * x = b by solving the preconditioned system inv(P) * A * inv(P)^T * y = inv(P) * b, 
        /// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/>, x is the solution, y = P^T * x,
        /// P*P^T = <paramref name="preconditioner"/>.
        /// Initially x = 0 and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
        /// <param name="rhsVector">The right hand side vector b of the linear system A * x = b. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="preconditioner">A preconditioner matrix that is also symmetric positive definite and has the same 
        ///     dimensions as A.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVector"/> violates the described 
        ///     constraint.</exception>
        public (IVector solution, CGStatistics stats) Solve(IMatrixView matrix, IVectorView rhs, IPreconditioner preconditioner)
        {
            CheckInput(matrix, rhs);
            Vector sol = Vector.CreateZero(matrix.NumColumns); // Start from x = 0
            IVector res = rhs.Copy(); // No need to do the multiplication A*0 = 0
            return SolveInternal(matrix, preconditioner, sol, res);
        }

        /// <summary>
        /// Solves the linear system A * x = b by solving the preconditioned system inv(P) * A * inv(P)^T * y = inv(P) * b, 
        /// where A = <paramref name="matrix"/>, b = <paramref name="rhsVector"/>, x is the solution, y = P^T * x,
        /// P*P^T = <paramref name="preconditioner"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
        /// <param name="rhsVector">The right hand side vector b of the linear system A * x = b. Constraints:
        ///     <paramref name="rhsVector"/>.<see cref="IIndexable1D.Length"/> == 
        ///     <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.</param>
        /// <param name="preconditioner">A preconditioner matrix that is also symmetric positive definite and has the same 
        ///     dimensions as A.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rhsVector"/> violates the described 
        ///     constraint.</exception>
        public (IVector solution, CGStatistics stats) Solve(IMatrixView matrix, IVectorView rhs, IPreconditioner preconditioner, 
            IVectorView initialGuess)
        {
            CheckInput(matrix, rhs, initialGuess);
            IVector sol = initialGuess.Copy(); // Should I copy this?
            IVector res = rhs.Subtract(matrix.MultiplyRight(sol));
            return SolveInternal(matrix, preconditioner, sol, res);
        }

        private (IVector solution, CGStatistics stats) SolveInternal(IMatrixView matrix, IPreconditioner preconditioner, 
            IVector sol, IVector res)
        {
            IVector z = preconditioner.SolveLinearSystem(res);
            IVector dir = z.Copy(); // TODO: Do I need to copy it?
            double zrDotCurrent = z.DotProduct(res);
            double resNormInit = res.Norm2(); // In basic CG, I could just take the sqrt(r*r), but here I have z*r.
            double resNormRatio = 1.0;
            CGStatistics statistics;

            for (int i = 0; i < maxIterations; ++i)
            {
                IVector matrixTimesDir = matrix.MultiplyRight(dir);
                double step = zrDotCurrent / (dir.DotProduct(matrixTimesDir));
                sol.AxpyIntoThis(dir, step);
                res.AxpyIntoThis(matrixTimesDir, - step);

                resNormRatio = Math.Sqrt(res.Norm2()) / resNormInit;
                if (resNormRatio < residualTolerance) // resNormRatio is non negative
                {
                    statistics = new CGStatistics
                    {
                        AlgorithmName = "Preconditioned Conjugate Gradient",
                        HasConverged = true, IterationsRequired = i + 1, NormRatio = resNormRatio
                    };
                    return (sol, statistics);
                }

                z = preconditioner.SolveLinearSystem(res); 
                double zrDotNext = z.DotProduct(res); //Fletcher-Reeves formula. TODO: For variable preconditioning use Polak-Ribiere
                double beta = zrDotNext / zrDotCurrent;
                dir = z.Axpy(dir, beta);
                zrDotCurrent = zrDotNext;
            }

            statistics = new CGStatistics
            {
                HasConverged = false, IterationsRequired = maxIterations, NormRatio = resNormRatio
            };
            return (sol, statistics);
        }

        private static void CheckInput(IMatrixView matrix, IVectorView rhs, IVectorView initialGuess)
        {
            if (matrix.NumColumns != initialGuess.Length) throw new NonMatchingDimensionsException(
                $"The matrix is {matrix.NumRows}-by-{matrix.NumColumns}),"
                + $" while the right hand side vector is {initialGuess.Length}-by-1");
            CheckInput(matrix, rhs);
        }

        private static void CheckInput(IMatrixView matrix, IVectorView rhs)
        {
            if (matrix.NumRows != matrix.NumColumns) throw new NonMatchingDimensionsException(
                "The matrix must be square, symmetric and positive definite.");
            Preconditions.CheckSystemSolutionDimensions(matrix, rhs);
        }
    }
}
