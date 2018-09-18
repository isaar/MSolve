using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms.CG;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;

//TODO: perhaps the user should choose the PCG settings himself and pass it. In this case, this should be named IterativeSolver.
//TODO: the maxIterations of PCG should be able to use the order of the matrix as a default value.
//TODO: IIndexable2D is not a good choice if all solvers must cast it to the matrix types the operate on.
namespace ISAAR.MSolve.Solvers.PCG
{
    /// <summary>
    /// Iterative solver for models with only 1 subdomain. Uses the Proconditioned Conjugate Gradient algorithm.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgSolver
    {
        private const string name = "PcgSolver"; // for error messages
        private readonly PreconditionedConjugateGradient pcgAlgorithm;
        private readonly IPreconditionerBuilder preconditionerBuilder;
        private IMatrixView matrix;
        private IPreconditioner preconditioner;

        public PcgSolver(int maxIterations, double residualTolerance, IPreconditionerBuilder preconditionerBuilder)
        {
            pcgAlgorithm = new PreconditionedConjugateGradient(maxIterations, residualTolerance);
            this.preconditionerBuilder = preconditionerBuilder;
        }

        public void Initialize()
        {
        }

        /// <summary>
        /// Forces the solver to replace the previous linear system matrix for each subdomain with the new ones. Any processing 
        /// done on these matrices (e.g. factorization) will be repeated.
        /// </summary>
        /// <param name="subdomainMatrices"></param>
        public void SetLinearSystemMatrices(IReadOnlyList<IIndexable2D> subdomainMatrices)
        {
            if (subdomainMatrices.Count != 1) throw new InvalidSolverException(
                name + " only works when there is a single subdomain.");
            if (subdomainMatrices[0] is IMatrixView matrixMultipliable)
            {
                matrix = matrixMultipliable;
                preconditioner = preconditionerBuilder.BuildPreconditioner(matrix);
            }
            else throw new InvalidSolverException(name + " can only operate on matrices that can be multiplied with vectors.");
        }

        /// <summary>
        /// Solves the linear systems using the subdomain matrices stored for each subdomain.
        /// </summary>
        /// <param name="subdomainRhsVectors">The right hand side vectors for the linear systems of all subdomains. They must 
        ///     be provided in the same order as the subdomain matrices are provided in 
        ///     <see cref="SetLinearSystemMatrices(IReadOnlyList{IIndexable2D})"/>.</param>
        /// <returns></returns>
        public IReadOnlyList<Numerical.LinearAlgebra.Interfaces.IVector> Solve(IReadOnlyList<Vector> subdomainRhsVectors)
        {
            if (subdomainRhsVectors.Count != 1) throw new InvalidSolverException(
                name + " only works when there is a single subdomain.");

            (Vector solution, CGStatistics stats) = pcgAlgorithm.Solve(matrix, subdomainRhsVectors[0], preconditioner);
            return new Numerical.LinearAlgebra.Interfaces.IVector[] { solution.ToLegacyVector() };
        }
    }
}
