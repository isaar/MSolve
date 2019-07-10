using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Iterative.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Base abstract class for Preconditioned Conjugate Gradient implementations.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public abstract class PcgAlgorithmBase
    {
        protected readonly IPcgResidualConvergence convergence;
        protected readonly IMaxIterationsProvider maxIterationsProvider;
        protected readonly double residualTolerance;
        protected readonly IPcgResidualUpdater residualUpdater;

        protected IVector direction;
        protected int iteration;
        protected IVector matrixTimesDirection;
        protected double paramBeta;
        protected IVector precondResidual;
        protected double resDotPrecondRes;
        protected double resDotPrecondResOld;
        protected IVector residual;
        protected IVector solution;
        protected double stepSize;

        protected PcgAlgorithmBase(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
            IPcgResidualConvergence convergence, IPcgResidualUpdater residualUpdater)
        {
            this.residualTolerance = residualTolerance;
            this.maxIterationsProvider = maxIterationsProvider;
            this.convergence = convergence;
            this.residualUpdater = residualUpdater;
        }

        /// <summary>
        /// The direction vector d, used to update the solution vector: x = x + α * d
        /// </summary>
        public IVectorView Direction => direction;

        /// <summary>
        /// The current iteration of the algorithm. It belongs to the interval [0, maxIterations).
        /// </summary>
        public int Iteration => iteration;

        /// <summary>
        /// The matrix A of the linear system or another object that implements matrix-vector multiplications.
        /// </summary>
        public ILinearTransformation Matrix { get; private set; }

        /// <summary>
        /// The vector that results from <see cref="Matrix"/> * <see cref="Direction"/>.
        /// </summary>
        public IVectorView MatrixTimesDirection => matrixTimesDirection;

        /// <summary>
        /// The β parameter of Conjugate Gradient that ensures conjugacy between the direction vectors.
        /// </summary>
        public double ParamBeta => paramBeta;

        /// <summary>
        /// The preconditioner M, such that inv(M) ~= inv(A).
        /// </summary>
        public IPreconditioner Preconditioner { get; private set; }

        /// <summary>
        /// The vector s = inv(M) * r
        /// </summary>
        public IVectorView PrecondResidual => precondResidual;

        /// <summary>
        /// The dot product r(t) * (inv(M) * r(t)) of the current iteration t.
        /// </summary>
        public double ResDotPrecondRes => resDotPrecondRes;

        /// <summary>
        /// The dot product r(t-1) * (inv(M) * r(t-1)) of the previous iteration t-1.
        /// </summary>
        public double ResDotPrecondResOld => resDotPrecondResOld;

        /// <summary>
        /// The residual vector r = b - A * x.
        /// </summary>
        public IVectorView Residual => residual;

        /// <summary>
        /// The right hand side of the linear system b = A * x.
        /// </summary>
        public IVectorView Rhs { get; private set; }

        /// <summary>
        /// The current approximation to the solution of the linear system A * x = b
        /// </summary>
        public IVectorView Solution => solution;

        /// <summary>
        /// The step α taken along <see cref="Direction"/> to update the solution vector: x = x + α * d
        /// </summary>
        public double StepSize => stepSize;

        /// <summary>
        /// Releases references to the vectors and matrices used by this object and sets scalars to their default values.
        /// </summary>
        public virtual void Clear()
        {
            Matrix = null;
            Rhs = null;
            solution = null;
            residual = null;
            direction = null;
            matrixTimesDirection = null;
            Preconditioner = null;
            precondResidual = null;
            resDotPrecondRes = 0.0;
            resDotPrecondResOld = 0.0;
            stepSize = 0.0;
            paramBeta = 0.0;
            iteration = -1;
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
        public IterativeStatistics Solve(IMatrixView matrix, IPreconditioner preconditioner, IVectorView rhs, IVector solution,
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
        public IterativeStatistics Solve(ILinearTransformation matrix, IPreconditioner preconditioner, IVectorView rhs,
            IVector solution, bool initialGuessIsZero, Func<IVector> zeroVectorInitializer)
        {
            //TODO: find a better way to handle optimizations for the case x0=0, than using an initialGuessIsZero flag

            Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
            Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

            this.Matrix = matrix;
            this.Preconditioner = preconditioner;
            this.Rhs = rhs;
            this.solution = solution;

            // r = b - A * x
            if (initialGuessIsZero) residual = rhs.Copy();
            else residual = ExactResidual.Calculate(matrix, rhs, solution);
            return SolveInternal(maxIterationsProvider.GetMaxIterations(matrix.NumColumns), zeroVectorInitializer);
        }

        protected abstract IterativeStatistics SolveInternal(int maxIterations, Func<IVector> zeroVectorInitializer);
    }
}
