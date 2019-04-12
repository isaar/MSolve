using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Iterative.Termination;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: needs to throw exceptions or at least report indefinite, nonsymmetric and singular matrices.
//TODO: zero vector initialization should be done by a vector factory
//TODO: Exposing properties is more flexible than pushing data to the strategies, but how can I ensure that the properties are
//      initialized when the strategies will access them?
namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
{
    /// <summary>
    /// Implements the Conjugate Gradient algorithm for solving linear systems with a positive definite matrix.
    /// The implementation is based on the algorithm presented in section B2 of 
    /// "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", Jonathan Richard Shewchuk, 1994
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CGAlgorithm
    {
        private const string name = "Conjugate Gradient";
        private readonly IMaxIterationsProvider maxIterationsProvider;
        private readonly ICGResidualConvergence residualConvergence;
        private readonly ICGResidualUpdater residualUpdater;
        private readonly double residualTolerance;

        private IVector direction;
        private IVector matrixTimesDirection;
        private double resDotRes;
        private IVector residual;
        private IVector solution;

        private CGAlgorithm(double residualTolerance, IMaxIterationsProvider maxIterationsProvider,
            ICGResidualConvergence residualConvergence, ICGResidualUpdater residualUpdater)
        {
            this.residualTolerance = residualTolerance;
            this.maxIterationsProvider = maxIterationsProvider;
            this.residualConvergence = residualConvergence;
            this.residualUpdater = residualUpdater;
        }

        /// <summary>
        /// The direction vector d, used to update the solution vector: x = x + α * d
        /// </summary>
        public IVectorView Direction => direction;

        /// <summary>
        /// The current iteration of the algorithm. It belongs to the interval [0, maxIterations).
        /// </summary>
        public int Iteration { get; private set; }

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
        public double ParamBeta { get; private set; }

        /// <summary>
        /// The dot product <see cref="Residual"/> * <see cref="Residual"/>.
        /// </summary>
        public double ResDotRes => resDotRes;

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
        public double StepSize { get; private set; }

        /// <summary>
        /// Releases references to the vectors and matrices used by this object and sets scalars to their default values.
        /// </summary>
        public void Clear()
        {
            Matrix = null;
            Rhs = null;
            solution = null;
            residual = null;
            direction = null;
            matrixTimesDirection = null;
            resDotRes = 0.0;
            StepSize = 0.0;
            ParamBeta = 0.0;
            Iteration = -1;
        }

        /// <summary>
        /// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhs"/>.
        /// Initially x = <paramref name="initialGuess"/> and then it converges to the solution.
        /// </summary>
        /// <param name="matrix">The matrix A of the linear system A * x = b. It must be symmetric positive definite.</param>
        /// <param name="rhs">
        /// The right hand side vector b of the linear system A * x = b. Constraints:
        /// <paramref name="rhs"/>.<see cref="IIndexable1D.Length"/> 
        /// == <paramref name="matrix"/>.<see cref="IIndexable2D.NumRows"/>.
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
        public IterativeStatistics Solve(IMatrixView matrix, IVectorView rhs, IVector solution, bool initialGuessIsZero) //TODO: find a better way to handle the case x0=0
            => Solve(new ExplicitMatrixTransformation(matrix), rhs, solution, initialGuessIsZero);

        /// <summary>
        /// Solves the linear system A * x = b, where A = <paramref name="matrix"/> and b = <paramref name="rhs"/>.
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
        public IterativeStatistics Solve(ILinearTransformation matrix, IVectorView rhs, IVector solution, 
            bool initialGuessIsZero) //TODO: find a better way to handle the case x0=0
        {
            //TODO: these will also be checked by the matrix vector multiplication.
            Preconditions.CheckMultiplicationDimensions(matrix.NumColumns, solution.Length);
            Preconditions.CheckSystemSolutionDimensions(matrix.NumRows, rhs.Length);

            this.Matrix = matrix;
            this.Rhs = rhs;
            this.solution = solution;

            // r = b - A * x
            if (initialGuessIsZero) residual = rhs.Copy();
            else residual = ExactResidual.Calculate(matrix, rhs, solution);

            return SolveInternal(maxIterationsProvider.GetMaxIterations(matrix.NumColumns));
        }

        private IterativeStatistics SolveInternal(int maxIterations)
        {
            // δnew = δ0 = r * r
            resDotRes = residual.DotProduct(residual);

            // The convergence criterion must be initialized immediately after the first r and r*r are computed.
            residualConvergence.Initialize(this);

            // This is also used as output
            double residualNormRatio = double.NaN;

            // d = r
            direction = residual.Copy();

            // Allocate memory for other vectors, which will be reused during each iteration
            matrixTimesDirection = Rhs.CreateZeroVectorWithSameFormat();

            for (Iteration = 0; Iteration < maxIterations; ++Iteration)
            {
                // q = A * d
                Matrix.Multiply(direction, matrixTimesDirection);

                // α = δnew / (d * q)
                StepSize = ResDotRes / (direction.DotProduct(matrixTimesDirection));

                // x = x + α * d
                solution.AxpyIntoThis(direction, StepSize);

                // δold = δnew
                double resDotResOld = ResDotRes;

                // Normally the residual vector is updated as: r = r - α * q and δnew = r * r. 
                // However corrections might need to be applied.
                residualUpdater.UpdateResidual(this, residual, out resDotRes);

                // At this point we can check if CG has converged and exit, thus avoiding the uneccesary operations that follow.
                residualNormRatio = residualConvergence.EstimateResidualNormRatio(this);
                if (residualNormRatio <= residualTolerance)
                {
                    return new IterativeStatistics
                    {
                        AlgorithmName = name,
                        HasConverged = true,
                        NumIterationsRequired = Iteration + 1,
                        ResidualNormRatioEstimation = residualNormRatio
                    };
                }

                // β = δnew / δold
                ParamBeta = ResDotRes / resDotResOld;

                // d = r + β * d
                //TODO: benchmark the two options to find out which is faster
                //direction = residual.Axpy(direction, beta); //This allocates a new vector d, copies r and GCs the existing d.
                direction.LinearCombinationIntoThis(ParamBeta, residual, 1.0); //This performs additions instead of copying and needless multiplications.
            }

            // We reached the max iterations before CG converged
            return new IterativeStatistics
            {
                AlgorithmName = name,
                HasConverged = false,
                NumIterationsRequired = maxIterations,
                ResidualNormRatioEstimation = residualNormRatio
            };
        }

        /// <summary>
        /// Constructs <see cref="CGAlgorithm"/> instances, allows the user to specify some or all of the required parameters and 
        /// provides defaults for the rest.
        /// Author: Serafeim Bakalakos
        /// </summary>
        public class Builder
        {
            /// <summary>
            /// Specifies how to calculate the maximum iterations that the CG algorithm will run for.
            /// </summary>
            public IMaxIterationsProvider MaxIterationsProvider { get; set; } = new PercentageMaxIterationsProvider(1.0);

            /// <summary>
            /// Specifies how the CG algorithm will check that convergence has been reached.
            /// </summary>
            public ICGResidualConvergence ResidualConvergence { get; set; } = new RegularCGConvergence();

            /// <summary>
            /// Specifies how often the residual vector will be corrected by an exact (but costly) calculation.
            /// </summary>
            public ICGResidualUpdater ResidualUpdater { get; set; } = new RegularCGResidualUpdater();

            /// <summary>
            /// Normally the CG will converge when norm2(r) / norm2(r0) &lt;= <paramref name="ResidualTolerance"/>, 
            /// where r = A*x is the current residual vector and r0 = A*x0 the initial residual vector. Depending on 
            /// <see cref="ResidualConvergence"/>, some other criterion might be used.
            /// </summary>
            public double ResidualTolerance { get; set; } = 1E-10;

            /// <summary>
            /// Creates a new instance of <see cref="CGAlgorithm"/>.
            /// </summary>
            public CGAlgorithm Build()
                => new CGAlgorithm(ResidualTolerance, MaxIterationsProvider, ResidualConvergence, ResidualUpdater);
        }
    }
}
