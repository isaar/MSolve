using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: reduce indexing the skyline arrays by incrementing/decrementing the offsets of previous iterations as much as possible
//TODO: use BLAS for the dot product / axpy during the factorization and the back/forward solve. Perhaps some cutoffs are needed.
//TODO: I think that top to bottom skyline format would have better memory access patterns, since the algorithm moves top to
//      bottom. It is also used by MKL so it will be easier to abstract the providers.
//TODO: The implementations of the factorization, forward and back substitution belong to SparseBLAS, SparseLAPACK providers.
//TODO: The nullspace should be a different class and the user should choose if he wants to calculate it.
//TODO: Reduce code duplication between this and CholeskySkyline
//TODO: In the examples tested so far, the dependent columns turn out to be the last 3 or 6 columns. Does this always happen?
//      Does this factorization work even if a dependent column is found before an independent one?
namespace ISAAR.MSolve.LinearAlgebra.Triangulation
{
    /// <summary>
    /// Cholesky-like factorization of a symmetric positive semi-definite matrix: A = L * transpose(L) = transpose(U) * U.
    /// This factorization will apply the Cholesky algorithm on the independent columns of the matrix, set the dependent ones 
    /// equal to columns of the identity matrix and return the nullspace of the matrix.
    /// The matrix is stored in skyline format. Only the active columns of the upper triangle part of the matrix is stored and 
    /// factorized. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SemidefiniteCholeskySkyline
    {
        /// <summary>
        /// The default value under which a diagonal entry (pivot) is considered to be 0 during Cholesky factorization.
        /// For the default value see "The FETI Level 1 Method: Theory and Implementation, Kamath Chandrika, 2000"
        /// </summary>
        public const double PivotTolerance = 1e-7; //

        private readonly double[] values;
        private readonly int[] diagOffsets;
        private readonly List<double[]> nullSpaceBasis; //TODO: should these be sparse vectors?

        private SemidefiniteCholeskySkyline(int order, double[] values, int[] diagOffsets, 
            List<int> dependentColumns, List<double[]> nullSpaceBasis)
        {
            this.Order = order;
            this.values = values;
            this.diagOffsets = diagOffsets;
            this.DependentColumns = dependentColumns;
            this.nullSpaceBasis = nullSpaceBasis;
        }

        public IReadOnlyList<int> DependentColumns { get; }

        public IReadOnlyList<double[]> NullSpaceBasis => nullSpaceBasis; //TODO: return IVectorView

        /// <summary>
        /// The number of rows/columns of the original square matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// Applies the Cholesky factorization to the independent columns of a symmetric positive semi-definite matrix,
        /// sets the dependent ones equal to columns of the identity matrix and return the nullspace of the matrix. Requires 
        /// extra memory for the basis vectors of the nullspace.
        /// </summary>
        /// <param name="order">The number of rows/ columns of the square matrix.</param>
        /// <param name="skyValues">
        /// The non-zero entries of the original <see cref="SkylineMatrix"/>. This array will be overwritten during the 
        /// factorization.
        /// </param>
        /// <param name="skyDiagOffsets">
        /// The indexes of the diagonal entries into <paramref name="skyValues"/>. The new 
        /// <see cref="SemidefiniteCholeskySkyline"/> instance will hold a reference to <paramref name="skyDiagOffsets"/>. 
        /// However they do not need copying, since they will not be altered during or after the factorization.
        /// </param>
        /// <param name="pivotTolerance">
        /// If a diagonal entry is &lt;= <paramref name="pivotTolerance"/> it means that the corresponding column is dependent 
        /// on the rest. The Cholesky factorization only applies to independent column, while dependent ones are used to compute
        /// the nullspace. Therefore it is important to select a tolerance that will identify small pivots that result from 
        /// singularity, but not from ill-conditioning.
        /// </param>
        public static SemidefiniteCholeskySkyline Factorize(int order, double[] skyValues, int[] skyDiagOffsets,
            double pivotTolerance = SemidefiniteCholeskySkyline.PivotTolerance)
        {
            (List<int> dependentColumns, List<double[]> nullSpaceBasis) = 
                FactorizeInternal(order, skyValues, skyDiagOffsets, pivotTolerance);
            Debug.Assert(dependentColumns.Count == nullSpaceBasis.Count);
            return new SemidefiniteCholeskySkyline(order, skyValues, skyDiagOffsets, dependentColumns, nullSpaceBasis);
        }

        /// <summary>
        /// See <see cref="ITriangulation.CalcDeterminant"/>.
        /// </summary>
        public double CalcDeterminant()
        {
            if (DependentColumns.Count > 0) return 0.0;
            else throw new NotImplementedException(); //TODO: call the code of the regular CholeskySkyline
        }

        /// <summary>
        /// Performs the operation: <paramref name="result"/> = generalized_inverse(A) * <paramref name="vector"/>. 
        /// The resulting vector overwrites <paramref name="result"/>.
        /// </summary>
        /// <param name="vector">The vector that will be multiplied. Its <see cref="IIndexable1D.Length"/> must be equal to 
        /// <see cref="IIndexable2D.NumRows"/> of the original matrix A.
        /// </param>
        /// <param name="result">
        /// Output vector that will be overwritten with the solution of the linear system. Its <see cref="IIndexable1D.Length"/>  
        /// must be equal to <see cref="IIndexable2D.NumColumns"/> of the original matrix A.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="vector"/> or <paramref name="result"/> violate the described constraints.
        /// </exception>
        public void MultiplyGeneralizedInverseMatrixTimesVector(Vector vector, Vector result)
        {
            Preconditions.CheckSystemSolutionDimensions(Order, vector.Length);
            Preconditions.CheckMultiplicationDimensions(Order, result.Length);

            // A^+ * b = [ Aii^-1 * bi; 0], where i are the independent rows/columns
            // TODO: Is this correct?
            CholeskySkyline.SubstituteForward(Order, values, diagOffsets, vector.RawData, result.RawData);
            CholeskySkyline.SubstituteBack(Order, values, diagOffsets, result.RawData);
            foreach (int row in DependentColumns) result[row] = 0.0;
        }

        /// <summary>
        /// Performs the operation: result = generalized_inverse(A) * <paramref name="vector"/>. The resul is written to a new
        /// vector and returned.
        /// </summary>
        /// <param name="vector">The vector that will be multiplied. Its <see cref="IIndexable1D.Length"/> must be equal to 
        /// <see cref="IIndexable2D.NumRows"/> of the original matrix A.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if <paramref name="vector"/> violates the described constraint.
        /// </exception>
        public Vector MultiplyGeneralizedInverseMatrixTimesVector(Vector vector)
        {
            var result = Vector.CreateZero(Order);
            MultiplyGeneralizedInverseMatrixTimesVector(vector, result);
            return result;
        }

        internal static (List<int> dependentColumns, List<double[]> nullSpaceBasis) FactorizeInternal(int order, double[] values, 
            int[] diagOffsets, double pivotTolerance)
        {
            var dependentColumns = new List<int>();
            var nullSpaceBasis = new List<double[]>();

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
                if (subroot > pivotTolerance) values[offsetAjj] = Math.Sqrt(subroot); // positive definite
                else if (subroot < -pivotTolerance) // not positive semidefinite
                {
                    throw new IndefiniteMatrixException($"The leading minor of order {j} (and therefore the matrix itself)"
                        + " is negative, and the factorization could not be completed.");
                }
                else // positive semidefinite with rank deficiency
                {
                    // Set the dependent column to identity and partially calculate a new basis vector of the null space.
                    dependentColumns.Add(j);
                    var basisVector = new double[order];
                    nullSpaceBasis.Add(basisVector);

                    // The diagonal entries of the skyline column and the null space basis vector are set to 1. 
                    basisVector[j] = 1.0;
                    values[offsetAjj] = 1.0;

                    // Superdiagonal entries are set to 0 after copying them to the nullspace (and negating them). 
                    for (int t = 1; t <= heightColJ; ++t) //TODO: indexing could be written more efficiently
                    {
                        basisVector[j - t] = -values[offsetAjj + t];
                        values[offsetAjj + t] = 0;
                    }

                    // Subdiagonal entries are also set to 0, but they are stored in subsequent columns i of the skyline format.
                    for (int i = j + 1; i < order; ++i)
                    {
                        int offsetAij = diagOffsets[i] + i - j; // TODO: reuse the indexing from the previous iteration
                        if (offsetAij < diagOffsets[i + 1]) values[offsetAij] = 0; //TODO: originally this was <= ? Which is the correct one?
                    }
                }
            }

            // The independent columns have been factorized successfully. Apply back substitution to finish the calculation
            // of each vector in the null space basis.
            foreach (var basisVector in nullSpaceBasis)
            {
                CholeskySkyline.SubstituteBack(order, values, diagOffsets, basisVector);
            }
            return (dependentColumns, nullSpaceBasis);
        }
    }
}
