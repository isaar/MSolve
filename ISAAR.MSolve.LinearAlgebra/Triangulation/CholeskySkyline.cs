using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: reduce indexing the skyline arrays by incrementing/decrementing the offsets of previous iterations as much as possible
//TODO: use BLAS for the dot product / axpy during the factorization and the back/forward solve. Perhaps some cutoffs are needed.
//TODO: I think that top to bottom skyline format would have better memory access patterns, since the algorithm moves top to
//      bottom. It is also used by MKL so it will be easier to abstract the providers.
//TODO: The implementations of the factorization, forward and back substitution belong to SparseBLAS, SparseLAPACK providers.
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// Cholesky factorization of a symmetric positive definite matrix: A = L * transpose(L) = transpose(U) * U. The matrix is  
    /// stored in skyline format. Only the active columns of the upper triangle part of the matrix is stored and factorized. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CholeskySkyline : SkylineFactorizationBase
    {
        /// <summary>
        /// The default value under which a diagonal entry (pivot) is considered to be 0 during Cholesky factorization.
        /// </summary>
        public const double PivotTolerance = 1e-15;

        private CholeskySkyline(int order, double[] values, int[] diagOffsets) : base(order, values, diagOffsets)
        { }

        /// <summary>
        /// Calculates the Cholesky factorization of a symmetric positive definite matrix, such that A = transpose(U) * U. 
        /// Does not need any extra memory.
        /// </summary>
        /// <param name="order">The number of rows/ columns of the square matrix.</param>
        /// <param name="skyValues">
        /// The non-zero entries of the original <see cref="SkylineMatrix"/>. This array will be overwritten during the 
        /// factorization.
        /// </param>
        /// <param name="skyDiagOffsets">
        /// The indexes of the diagonal entries into <paramref name="skyValues"/>. The new <see cref="CholeskySkyline"/> 
        /// instance will hold a reference to <paramref name="skyDiagOffsets"/>. However they do not need copying, since they 
        /// will not be altered during or after the factorization.
        /// </param>
        /// <param name="pivotTolerance">
        /// If a diagonal entry is &lt;= <paramref name="pivotTolerance"/> it means that the original matrix is not symmetric 
        /// positive definite and an <see cref="IndefiniteMatrixException"/> will be thrown.
        /// </param>
        /// <exception cref="IndefiniteMatrixException">
        /// Thrown if the original skyline matrix turns out to not be symmetric positive definite.
        /// </exception>
        public static CholeskySkyline Factorize(int order, double[] skyValues, int[] skyDiagOffsets,
            double pivotTolerance = CholeskySkyline.PivotTolerance)
        {
            FactorizeInternal(order, skyValues, skyDiagOffsets, pivotTolerance);
            return new CholeskySkyline(order, skyValues, skyDiagOffsets);
        }

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        public override double CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Explicitly creates the upper triangular matrix U that resulted from the Cholesky factorization: A = transpose(U) * U,
        /// where A and U are n-by-n. 
        /// This method is safe to use as the factorization data are copied (if necessary). However, it is inefficient if the 
        /// generated matrix is only used once.
        /// </summary>
        public TriangularUpper GetFactorU()
        {
            // The factorization A = transpose(u) * D * u, u = unit upper triangular is stored. Thus U = sqrt(D) * u.
            // Since D is diagonal, we need to scale each column j of u by sqrt(D[j,j]).
            var upper = TriangularUpper.CreateZero(Order);
            for (int j = 0; j < Order; ++j)
            {
                int colheight = diagOffsets[j+1] - diagOffsets[j] - 1;
                for (int t = 0; t < colheight + 1; ++t)
                {
                    int i = j - t;
                    upper[i, j] = values[diagOffsets[j] + t];
                }
            }
            return upper;
        }

        /// <summary>
        /// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
        /// </summary>
        public override void SolveLinearSystem(Vector rhs, Vector solution)
        {
            Preconditions.CheckSystemSolutionDimensions(this, rhs);
            Preconditions.CheckMultiplicationDimensions(Order, solution.Length);

            SubstituteForward(Order, values, diagOffsets, rhs.RawData, solution.RawData);
            SubstituteBack(Order, values, diagOffsets, solution.RawData);
        }

        internal static void FactorizeInternal(int order, double[] values, int[] diagOffsets, double pivotTolerance)
        {
            // Process column j
            for (int j = 0; j < order; ++j)
            {
                int offsetAjj = diagOffsets[j];

                // The number of non-zero entries in column j, above the diagonal and excluding it
                int heightColJ = diagOffsets[j + 1] - offsetAjj - 1; //TODO: reuse the diagOffset form previous iteration.

                // The row index above which col j has only zeros
                int topColJ = j - heightColJ;

                // Update each A[i,j] by traversing the column j from the top until the diagonal. 
                // The top now is the min row with a non-zero entry instead of 0.
                for (int i = topColJ; i < j; ++i)
                {
                    int offsetAii = diagOffsets[i];         // TODO: increment/decrement the offset from previous iteration
                    int offsetAij = offsetAjj + j - i;      // TODO: increment/decrement the offset from previous iteration

                    // The number of non-zero entries in column i, above the diagonal and excluding it
                    int heightColI = diagOffsets[i + 1] - offsetAii - 1; //TODO: reuse the diagOffset form previous iteration.

                    // The row index above which col j has only zeros
                    int topColI = i - heightColI;

                    // The row index above which either col j or col i has only zeros: max(topColJ, topColI)
                    int topCommon = (topColI > topColJ) ? topColI : topColJ;

                    // Dot product of the parts of columns j, i (i < j) between: [the common top non-zero row, row i)
                    // for (int k = max(topRowOfColJ, topRowOfColI; k < i; ++k) dotColsIJ += A[k,i] * A[k,j]
                    int numDotEntries = i - topCommon;
                    double dotColsIJ = 0.0;
                    for (int t = 1; t <= numDotEntries; ++t) //TODO: BLAS dot
                    {
                        dotColsIJ += values[offsetAii + t] * values[offsetAij + t];
                    }

                    // A[i,j] = (A[i,j] - dotIJ) / A[i,i]
                    values[offsetAij] = (values[offsetAij] - dotColsIJ) / values[offsetAii];
                }

                // Update the diagonal term
                // Dot product with itself of the part of column j between: [the top non-zero row, row j).
                // for (int k = topRowOfColJ; k < j; ++k) dotColsJJ += A[k,j]^2
                double dotColsJJ = 0.0;
                double valueAkj;
                for (int t = 1; t <= heightColJ; ++t) //TODO: BLAS dot. //TODO: the managed version would be better off using the offsetAkj as a loop variable
                { 
                    valueAkj = values[offsetAjj + t];
                    dotColsJJ += valueAkj * valueAkj;
                }

                // if A[j,j] = sqrt(A[j,j]-dotColsJJ), but if the subroot is <= 0, then the matrix is not positive definite
                double subroot = values[offsetAjj] - dotColsJJ;
                if (subroot <= pivotTolerance)
                {
                    throw new IndefiniteMatrixException($"The leading minor of order {j} (and therefore the matrix itself)"
                        + " is not positive-definite, and the factorization could not be completed.");
                }
                values[offsetAjj] = Math.Sqrt(subroot);
            }
        }

        internal static void SubstituteBack(int order, double[] values, int[] diagOffsets, double[] sol)
        {
            // Column / vector version of the algorithm: process U column-by-column 
            for (int j = order - 1; j >= 0; --j)
            {
                int offsetUjj = diagOffsets[j];

                // The number of non-zero entries in column j, above the diagonal and excluding it
                int heightColJ = diagOffsets[j + 1] - offsetUjj - 1;

                // Update the diagonal term first
                sol[j] /= values[offsetUjj];

                // Scale and subtract the column i above the diagonal from x
                // for (int i = 0; i < j; ++i) x[i] += diagValue * U[i, j];
                double diagValue = -sol[j];
                for (int t = 1; t <= heightColJ; ++t) //TODO: BLAS axpy(), but needs top-to-bottom skyline or Array.Reverse(sol)
                {
                    sol[j - t] += diagValue * values[offsetUjj + t];
                }
            }
        }

        internal static void SubstituteForward(int order, double[] values, int[] diagOffsets, double[] rhs, double[] sol)
        {
            // Dot product version of the algorithm: process L column-by-column <=> process U row-by-row
            for (int i = 0; i < order; ++i)
            {
                int offsetUii = diagOffsets[i];

                // The number of non-zero entries in column i, above the diagonal and excluding it
                int heightColI = diagOffsets[i + 1] - offsetUii - 1;

                // for (int j = 0; j < i; ++j) dot += L[i,j] * x[j], where L[i,j] = U[j,i]
                double dot = 0;
                for (int t = 1; t <= heightColI; ++t) //TODO: BLAS dot(), but needs top-to-bottom skyline or Array.Reverse(sol)
                {
                    dot += values[offsetUii + t] * sol[i - t];
                }

                // x[i] = (b[i] - dot) / L[i,i]
                sol[i] = (rhs[i] - dot) / values[offsetUii];
            }
        }
    }
}
