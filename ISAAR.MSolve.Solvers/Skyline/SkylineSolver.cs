using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;

//TODO: factorization should be done during Solve() to time correctly. Alternatively, an observer should record the durations.
//TODO: investigate if it is possible to avoid casting the matrix provided by the analyzer/assembler into skyline. Perhaps the
//      the matrix could be obtained directly from the assembler and the analyzer could provide delegates with the operations it 
//      wants done on the matrix, instead of doing them itself.
//TODO: try to abstract the subdomain logic from ther analyzers. I doubt it is possible though.
namespace ISAAR.MSolve.Solvers.Skyline
{
    /// <summary>
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on symmetric positive definite matrices 
    /// stored in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineSolver
    {
        private const string name = "SkylineSolver"; // for error messages
        private readonly double factorizationPivotTolerance;
        private CholeskySkyline factorizedMatrix;

        public SkylineSolver(double factorizationPivotTolerance = 1E-15)
        {
            this.factorizationPivotTolerance = factorizationPivotTolerance;
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
            if (subdomainMatrices[0] is SkylineMatrix skyline)
            {
                factorizedMatrix = skyline.FactorCholesky(true, factorizationPivotTolerance);
            }
            else throw new InvalidSolverException(name + " can only operate on matrices stored in Skyline format.");
        }

        /// <summary>
        /// Solves the linear systems using the subdomain matrices stored for each subdomain.
        /// </summary>
        /// <param name="subdomainRhsVectors">The right hand side vectors for the linear systems of all subdomains. They must 
        ///     be provided in the same order as the subdomain matrices are provided in 
        ///     <see cref="SetLinearSystemMatrices(IReadOnlyList{IIndexable2D})"/>.</param>
        /// <returns></returns>
        public IReadOnlyList<Vector> Solve(IReadOnlyList<Vector> subdomainRhsVectors)
        {
            if (subdomainRhsVectors.Count != 1) throw new InvalidSolverException(
                name + " only works when there is a single subdomain.");
            Vector solution = factorizedMatrix.SolveLinearSystem(subdomainRhsVectors[0]);
            return new Vector[] { solution };
        }
    }
}
