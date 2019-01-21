using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Triangulation.SampleImplementations;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: reduce indexing the skyline arrays by incrementing/decrementing the offsets of previous iterations as much as possible
//TODO: use BLAS for the dot products during the factorization and the back/forward solve. Perhaps some cutoffs are needed.
//TODO: I think that top to bottom skyline format would have better memory access patterns, since the algorithm moves top to
//      bottom. It is also used by MKL so it will be easier to abstract the providers.
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// Cholesky factorization of a symmetric positive definite matrix: A = L * transpose(L) = transpose(U) * U. The matrix is  
    /// stored in skyline format. Only the active columns of the upper triangle part of the matrix is stored and factorized. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CholeskySkyline : IIndexable2D, ISparseMatrix, ITriangulation
    {
        /// <summary>
        /// The default value under which a diagonal entry (pivot) is considered to be 0 during Cholesky factorization.
        /// </summary>
        public const double PivotTolerance = 1e-15;

        private readonly double[] values;
        private readonly int[] diagOffsets;

        private CholeskySkyline(int order, double[] values, int[] diagOffsets)
        {
            this.NumColumns = order;
            this.values = values;
            this.diagOffsets = diagOffsets;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get { return NumColumns; } }

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>.
        /// </summary>
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                SkylineMatrix.ProcessIndices(ref rowIdx, ref colIdx);
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffsets[colIdx] - 1; // excluding diagonal
                if (entryHeight > maxColumnHeight) return 0.0; // outside stored non zero pattern
                else return values[diagOffsets[colIdx] + entryHeight];
            }
        }

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
            double pivotTolerance = LdlSkyline.PivotTolerance)
        {
            // Process column j
            for (int j = 0; j < order; ++j)
            {
                int offsetAjj = skyDiagOffsets[j];

                // The number of non-zero entries in column j, above the diagonal and excluding it
                int heightColJ = skyDiagOffsets[j + 1] - offsetAjj - 1; //TODO: reuse the diagOffset form previous iteration.

                // The row index above which col j has only zeros
                int topColJ = j - heightColJ;

                // Update each A[i,j] by traversing the column j from the top until the diagonal. 
                // The top now is the min row with a non-zero entry instead of 0.
                for (int i = topColJ; i < j; ++i)
                {
                    int offsetAii = skyDiagOffsets[i];         // TODO: increment/decrement the offset from previous iteration
                    int offsetAij = offsetAjj + j - i;      // TODO: increment/decrement the offset from previous iteration

                    // The number of non-zero entries in column i, above the diagonal and excluding it
                    int heightColI = skyDiagOffsets[i + 1] - offsetAii - 1; //TODO: reuse the diagOffset form previous iteration.

                    // The row index above which col j has only zeros
                    int topColI = i - heightColI;

                    // The row index above which either col j or col i has only zeros: max(topColJ, topColI)
                    int topCommon = (topColI > topColJ) ? topColI : topColJ;

                    // Dot product of the parts of columns j, i (i < j) between: [the common top non-zero row, row i)
                    // for (int k = max(topRowOfColJ, topRowOfColI; k < i; ++k) dotColsIJ += A[k,i] * A[k,j]
                    int numDotEntries = i - topCommon;
                    int startColI = offsetAii + 1;
                    int startColJ = offsetAij + 1;
                    double dotColsIJ = 0.0;
                    for (int t = 0; t < numDotEntries; ++t)
                    {
                        dotColsIJ += skyValues[startColI + t] * skyValues[startColJ + t];
                    }

                    // A[i,j] = (A[i,j] - dotIJ) / A[i,i]
                    skyValues[offsetAij] = (skyValues[offsetAij] - dotColsIJ) / skyValues[offsetAii];
                }

                // Update the diagonal term
                // Dot product with itself of the part of column j between: [the top non-zero row, row j).
                // for (int k = topRowOfColJ; k < j; ++k) dotColsJJ += A[k,j]^2
                double dotColsJJ = 0.0;
                int start = offsetAjj + 1;
                double valueAkj;
                for (int t = 0; t < heightColJ; ++t)
                {
                    valueAkj = skyValues[start + t];
                    dotColsJJ += valueAkj * valueAkj;
                }

                // if A[j,j] = sqrt(A[j,j]-dotColsJJ), but if the subroot is <= 0, then the matrix is not positive definite
                double subroot = skyValues[offsetAjj] - dotColsJJ;
                if (subroot < pivotTolerance)
                {
                    throw new IndefiniteMatrixException($"The leading minor of order {j} (and therefore the matrix itself)"
                        + " is not positive-definite, and the factorization could not be completed.");
                }
                skyValues[offsetAjj] = Math.Sqrt(subroot);
            }
            return new CholeskySkyline(order, skyValues, skyDiagOffsets);
        }

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        public double CalcDeterminant()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>.
        /// </summary>
        public int CountNonZeros() => values.Length;

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
            => SkylineMatrix.CreateFromArrays(NumColumns, values, diagOffsets, false, false).EnumerateNonZeros();

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
            => SkylineMatrix.CreateFromArrays(NumColumns, values, diagOffsets, false, false).Equals(other, tolerance);

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
            var upper = TriangularUpper.CreateZero(NumColumns);
            for (int j = 0; j < NumColumns; ++j)
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
        /// See <see cref="ISparseMatrix.GetSparseFormat"/>.
        /// </summary>
        public SparseFormat GetSparseFormat()
            => SkylineMatrix.CreateFromArrays(NumColumns, values, diagOffsets, false, false).GetSparseFormat();

        /// <summary>
        /// See <see cref="ITriangulation.SolveLinearSystem(Vector, Vector)"/>.
        /// </summary>
        public void SolveLinearSystem(Vector rhs, Vector solution)
        {
            Preconditions.CheckSystemSolutionDimensions(this, rhs);
            Preconditions.CheckMultiplicationDimensions(NumColumns, solution.Length);

            throw new NotImplementedException();
        }
    }
}
