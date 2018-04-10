using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public class SkylineMatrix: IMatrixView, ISparseMatrix
    {
        /// <summary>
        /// Contains the non zero superdiagonal entries of the matrix in column major order, starting from the diagonal and going
        /// upwards. Its length is nnz.
        /// </summary>
        private double[] values;

        /// <summary>
        /// Contains the indices into values of the diagonal entries of the matrix. Its length = order + 1, with the last entry
        /// being equal to nnz.
        /// </summary>
        private int[] diagOffsets;

        private SkylineMatrix(int order, double[] values, int[] diagOffsets)
        {
            this.values = values;
            this.diagOffsets = diagOffsets;
            this.NumColumns = order;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get { return NumColumns; } }

        // TODO: should I add index bound checking?
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffsets[colIdx] - 1; // excluding diagonal
                if (entryHeight > maxColumnHeight) return 0.0; // outside stored non zero pattern
                else return values[diagOffsets[colIdx] + entryHeight];
            }
            set
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffsets[colIdx] - 1; // excluding diagonal
                if (entryHeight > maxColumnHeight)
                {
                    throw new SparsityPatternModifiedException($"In column {colIdx} only rows [{maxColumnHeight}, {colIdx}]"
                        + $" can be changed, but you are trying to set entry ({rowIdx}, {colIdx})");
                }
                else values[diagOffsets[colIdx] + entryHeight] = value;
            }
        }

        public static SkylineMatrix CopyFromMatrix(SkylineMatrix originalMatrix)
        {
            return CreateFromArrays(originalMatrix.NumColumns, originalMatrix.values, originalMatrix.diagOffsets, true, false);
        }

        public static SkylineMatrix CreateFromArrays(int order, double[] values, int[] diagOffsets, 
            bool copyArrays, bool checkInput)
        {
            if (checkInput)
            {
                if (diagOffsets.Length != order + 1)
                {
                    throw new ArgumentException("The length of the Skyline diagonal offsets array must be equal to the number of"
                        + " rows/columns + 1, but was " + diagOffsets.Length);
                }
                if (diagOffsets[diagOffsets.Length - 1] != values.Length)
                {
                    throw new ArgumentException("The last entry of the Skyline diagonal offsets array must be equal to the number"
                        + " of non zero entries, but was " + diagOffsets[diagOffsets.Length - 1]);
                }
            }
            if (copyArrays)
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                int[] diagOffsetsCopy = new int[diagOffsets.Length];
                Array.Copy(diagOffsets, diagOffsetsCopy, diagOffsets.Length);
                return new SkylineMatrix(order, valuesCopy, diagOffsetsCopy);
            }
            else return new SkylineMatrix(order, values, diagOffsets);
        }

        public static SkylineMatrix CreateZeroWithPattern(int order, int[] diagOffsets, bool checkInput)
        {
            if (checkInput)
            {
                if (diagOffsets.Length != order + 1)
                {
                    throw new ArgumentException("The length of the Skyline diagonal offsets array must be equal to the number of"
                        + " rows/columns + 1, but was " + diagOffsets.Length);
                }
            }
            int nnz = diagOffsets[diagOffsets.Length] - 1;
            return new SkylineMatrix(order, new double[nnz], diagOffsets);
        }

        public double[,] CopyToArray2D()
        {
            double[,] array2D = new double[NumColumns, NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
                array2D[j, j] = values[colOffset]; // diagonal entry
                for (int i = columnTop; i < j; ++i) // non zero entries stored above diagonal
                {
                    double value = values[colOffset + j - i];
                    array2D[j, i] = value;
                    array2D[i, j] = value;
                }
            }
            return array2D;
        }

        public Matrix CopyToFullMatrix()
        {
            Matrix fullMatrix = Matrix.CreateZero(this.NumColumns, this.NumColumns);
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
                fullMatrix[j, j] = values[colOffset]; // diagonal entry
                for (int i = columnTop; i < j; ++i) // non zero entries stored above diagonal
                {
                    double value = values[colOffset + j - i];
                    fullMatrix[j, i] = value;
                    fullMatrix[i, j] = value; 
                }
            }
            return fullMatrix;
        }

        public int CountNonZeros()
        {
            return values.Length;
        }

        public IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is SkylineMatrix otherSKY) // In case both matrices have the exact same index arrays
            {
                if (this.diagOffsets == otherSKY.diagOffsets)
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherSKY.values[i]);
                    }
                    return new SkylineMatrix(NumColumns, resultValues, diagOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        public void DoEntrywiseIntoThis(SkylineMatrix other, Func<double, double, double> binaryOperation)
        {
            //Preconditions.CheckSameMatrixDimensions(this, other); // no need if the indexing arrays are the same
            if (this.diagOffsets != other.diagOffsets)
            {
                throw new SparsityPatternModifiedException("Only allowed if the index arrays are the same");
            }
            for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], other.values[i]);
        }

        IMatrixView IMatrixView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            // Only apply the operation on non zero entries
            double[] newValues = new double[values.Length];
            for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Copy the index arrays. TODO: See if we can use the same index arrays (e.g. if this class does not change them (it shouldn't))
                int[] diagOffsetsCopy = new int[diagOffsets.Length];
                Array.Copy(diagOffsets, diagOffsetsCopy, diagOffsets.Length);
                return new SkylineMatrix(NumColumns, newValues, diagOffsetsCopy);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new SkylineMatrix(NumColumns, newValues, diagOffsets).CopyToFullMatrix();
            }
        }

        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0))
            {
                for (int i = 0; i < values.Length; ++i) values[i] = unaryOperation(values[i]);
            }
            else
            {
                throw new SparsityPatternModifiedException("This operation will change the sparsity pattern");
            }
        }

        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
                yield return (j, j, values[colOffset]); // diagonal entry
                for (int i = columnTop; i < j; ++i)
                {
                    double value = values[colOffset + j - i];
                    yield return (i, j, value);
                    yield return (j, i, value);
                }
            }
        }

        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
            var comparer = new ValueComparer(1e-13);
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j+1] + colOffset + 1;
                for (int i = 0; i < columnTop; ++i) // zero entries above stored column
                {
                    if (!( comparer.AreEqual(0.0, other[i, j]) && comparer.AreEqual(0.0, other[j, i]) )) return false;
                }
                for (int i = columnTop; i < j; ++i) // non zero entries of column, excluding diafonal
                {
                    double value = values[colOffset + j - i];
                    if (!(comparer.AreEqual(value, other[i, j]) && comparer.AreEqual(value, other[j, i]))) return false;
                }
                if (!comparer.AreEqual(values[colOffset], other[j, j])) return false; // non zero diagonal entry
            }
            return true; // At this point all entries have been checked and are equal
        }

        /// <summary>
        /// Calculate the cholesky factorization. The matrix must be positive definite, otherwise an
        /// <see cref="IndefiniteMatrixException"/> will be thrown. If inPlace is set to true, this object must not be used 
        /// again, otherwise a <see cref="NullReferenceException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">False, to copy the internal non zero entries before factorization. True, to overwrite them with
        ///     the factorized data, thus saving memory and time. However, that will make this object unusable, so you MUST NOT 
        ///     call any other members afterwards.</param>
        /// <param name="tolerance">If a diagonal entry is closer to zero than this tolerance, an 
        ///     <see cref="IndefiniteMatrixException"/> exception will be thrown.</param>
        /// <returns></returns>
        public SkylineCholesky FactorCholesky(bool inPlace, double tolerance = SkylineCholesky.PivotTolerance)
        {
            if (inPlace)
            {
                var factor = SkylineCholesky.CalcFactorization(NumColumns, values, diagOffsets);
                // Set the skyline arrays to null to force NullReferenceException if they are accessed again.
                // TODO: perhaps there is a better way to handle this.
                values = null;
                diagOffsets = null;
                return factor;
            }
            else
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                return SkylineCholesky.CalcFactorization(NumColumns, values, diagOffsets);
            }
        }

        public SparseFormat GetSparseFormat()
        {
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Diagonal offsets", diagOffsets);
            return format;
        }

        /// <summary>
        /// Not optimized yet.
        /// </summary>
        /// <param name="other"></param>
        /// <param name="transposeThis">Does not matter since <see cref="SkylineMatrix"/> is symmetric.</param>
        /// <param name="transposeOther"></param>
        /// <returns></returns>
        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(other, this, transposeOther, transposeThis);
        }

        /// <summary>
        /// Not optimized yet.
        /// </summary>
        /// <param name="other"></param>
        /// <param name="transposeThis">Does not matter since <see cref="SkylineMatrix"/> is symmetric.</param>
        /// <param name="transposeOther"></param>
        /// <returns></returns>
        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(this, other, transposeThis, transposeOther);
        }

        public VectorMKL MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            if (vector is VectorMKL) return MultiplyRight((VectorMKL)vector, transposeThis);
            else throw new NotImplementedException();
        }

        /// <summary>
        /// Performance considerations: Only accesses the matrix entries once, but the result vector entries will be accessed 
        /// multiple times. 
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="transposeThis">Does not matter since <see cref="SkylineMatrix"/> is symmetric.</param>
        /// <returns></returns>
        public VectorMKL MultiplyRight(VectorMKL vector, bool transposeThis = false)
        {
            int n = vector.Length;
            Preconditions.CheckMultiplicationDimensions(NumColumns, n);
            double[] result = new double[n];
            // A*x = (L+D)*x + U*x
            // (L+D)*x is easy, since the non zero entries of row i left of the diagonal are stored contiguously in column i and
            // we can easily take its dot product with the vector.
            // U*x is trickier, since we cannot access contiguously the non zero entries of row i. Instead think of it as
            // U*x = linear combination of columns of U (accessed contiguously) with the entries of vector as coefficients. Then 
            // we can deal with them while we process the next columns (i, n-1]. This way the matrix is only indexed once, but 
            // not the result vector entry result[i].
            for (int j = 0; j < NumColumns; ++j)
            {
                int diagOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + diagOffset + 1;
                double linearCombinationCoeff = vector[j];
                // Dot product of the (L+D) part of the row * vector
                double dotLower = values[diagOffset] * linearCombinationCoeff; // Contribution of diagonal entry: A[j,j] * x[j]
                for (int i = j - 1; i >= columnTop; --i) // Process the rest of the non zero entries of the column
                {
                    double aij = values[diagOffset + j - i]; // Thus the matrix is only indexed once

                    // Contribution of the L part of the row, which is identical to the stored column j.
                    // Thus A[j,i]=A[i,j] and sum(A[i,j]*x[j]) = sum(A[i,j]*x[i])
                    dotLower += aij * vector[i];

                    // Contribution of the U part of the column: result += coefficient * column j of U. This will update all rows
                    // [columnTop, j) of the result vector need to be updated to account for the current column j. 
                    result[i] += aij * linearCombinationCoeff;
                }
                // Column j alters rows [0,j) of the result vector, thus this should be the 1st time result[j] is written.
                Debug.Assert(result[j] == 0);
                result[j] = dotLower; // contribution of the (L+D) part of the row. 
            }
            return VectorMKL.CreateFromArray(result, false);
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int nnz = values.Length;
            for (int i = 0; i < nnz; ++i) aggregator = processEntry(values[i], aggregator);
            aggregator = processZeros(NumColumns * NumColumns - nnz, aggregator);
            return finalize(aggregator);
        }

        public IMatrixView Transpose()
        {
            return Transpose(true);
        }

        public SkylineMatrix Transpose(bool copyInternalArrays)
        {
            return CreateFromArrays(NumColumns, values, diagOffsets, copyInternalArrays, false);
        }

        /// <summary>
        /// Perhaps this should be manually inlined. Testing needed.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static void ProcessIndices(ref int rowIdx, ref int colIdx)
        {
            if (rowIdx > colIdx)
            {
                int swap = rowIdx;
                rowIdx = colIdx;
                colIdx = swap;
            }
        }
    }
}
