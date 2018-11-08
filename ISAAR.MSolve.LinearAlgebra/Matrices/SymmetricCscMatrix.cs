using System;
using System.Runtime.CompilerServices;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.SuiteSparse;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Symmetric sparse matrix in Compressed Sparse Columns format, with only the non-zero entries of the upper triangle being 
    /// explicitly stored. Do not use this, since it is an experimantal class, which will probably be removed.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SymmetricCscMatrix: IIndexable2D
    {
        /// <summary>
        /// values.Length = number of non zeros in upper triangle.
        /// </summary>
        private readonly double[] values;

        /// <summary>
        /// rowIndices.Length = number of non zeros in upper triangle.
        /// </summary>
        private readonly int[] rowIndices;

        /// <summary>
        /// colOffsets.Length = number of rows/columns +1. 
        /// colOffsets[colOffsets.Length-1] = number of non zeros in upper triangle.
        /// </summary>
        private readonly int[] colOffsets;

        /// <summary>
        /// Does not copy anything
        /// </summary>
        /// <param name="values"></param>
        /// <param name="rowIndices"></param>
        /// <param name="colOffsets"></param>
        /// <param name="checkInput">Set to true to perform basic dimension checks.</param>
        public SymmetricCscMatrix(double[] values, int[] rowIndices, int[] colOffsets, bool checkInput)
        {
            int nnz = colOffsets[colOffsets.Length - 1];
            if (checkInput)
            {
                if ((nnz != values.Length) || (nnz != rowIndices.Length))
                {
                    throw new ArgumentException("Mismatch in dimensions of the CSC arrays. Check that colOffsets.Length ="
                        + " number of rows/columns + 1 and that colOffsets[colOffsets.Length-1] = values.Length ="
                        + " rowIndices.Length = number of non zeros in upper triangle");
                }
            }
            this.values = values;
            this.rowIndices = rowIndices;
            this.colOffsets = colOffsets;
            this.NumColumns = colOffsets.Length - 1;
            this.NumRows = colOffsets.Length - 1;
            this.NumNonZerosUpper = nnz;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// Only structural non zeros
        /// </summary>
        /// <returns></returns>
        public int NumNonZeros { get { throw new NotImplementedException(); } }

        public int NumNonZerosUpper { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        public double this[int rowIdx, int colIdx]
        {
            get
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                int colStart = colOffsets[colIdx];
                int colEnd = colOffsets[colIdx+1];
                for (int i = colStart; i < colEnd; ++i) //Only scan the nnz entries for the given column
                {
                    if (rowIndices[i] == rowIdx) return values[i];
                }
                return 0.0; //The requensted entry is not stored.
            }
        }

        public IIndexable2D DoPointwise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is SymmetricCscMatrix) // In case both matrices have the exact same index arrays
            {
                SymmetricCscMatrix otherCSC = (SymmetricCscMatrix)other;
                if (HaveSameIndexArrays(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherCSC.values[i]);
                    }
                    return new SymmetricCscMatrix(resultValues, rowIndices, colOffsets, false);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            Matrix result = Matrix.CreateZero(this.NumRows, this.NumColumns);
            for (int j = 0; j < this.NumColumns; ++j)
            {
                for (int i = 0; i < this.NumRows; ++i)
                {
                    result[i, j] = binaryOperation(this[i, j], other[i, j]);
                }
            }
            return result;
        }

        public void DoPointwiseIntoThis(SymmetricCscMatrix other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            if (!HaveSameIndexArrays(other))
            {
                throw new SparsityPatternModifiedException("Only allowed if the index arrays are the same");
            }
            for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], other.values[i]);
        }

        public IIndexable2D DoToAllEntries(Func<double, double> unaryOperation)
        {
            // Only apply the operation on non zero entries
            double[] newValues = new double[values.Length];
            for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Copy the index arrays. TODO: See if we can use the same index arrays (e.g. if this class does not change them (it shouldn't))
                int[] cloneRowIndices = new int[rowIndices.Length];
                Array.Copy(rowIndices, cloneRowIndices, rowIndices.Length);
                int[] cloneColOffsets = new int[colOffsets.Length];
                Array.Copy(colOffsets, cloneColOffsets, colOffsets.Length);
                return new SymmetricCscMatrix(newValues, cloneRowIndices, cloneColOffsets, false);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new SymmetricCscMatrix(newValues, rowIndices, colOffsets, false).ToFullMatrix();
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

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
        }

        public CholeskySuiteSparse FactorCholesky(SuiteSparseOrdering ordering)
        {
            return CholeskySuiteSparse.Factorize(NumColumns, NumNonZerosUpper, values, rowIndices, colOffsets, false, ordering);
        }

        /// <summary>
        /// 2 or more sparse matrices may have the same (2 references to 1 object) index arrays. 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool HaveSameIndexArrays(SymmetricCscMatrix other)
        {
            return (this.rowIndices == other.rowIndices) && (this.colOffsets == other.colOffsets);
        }

        /// <summary>
        /// Matrix-vector multiplication, with the vector on the right: matrix * vector or transpose(matrix) * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
        /// <param name="transposeThis">Set to true to transpose this (the left matrix). Unless the transpose matrix is used in 
        ///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
        /// <returns></returns>
        public Vector MultiplyRight(Vector vector, bool transposeThis = false)
        {
            throw new NotImplementedException();
        }

        public Matrix ToFullMatrix()
        {
            Matrix fullMatrix = Matrix.CreateZero(this.NumRows, this.NumColumns);
            for (int j = 0; j < this.NumColumns; ++j) //Column major order
            {
                int colCurrent = colOffsets[j];
                int colNext = colOffsets[j + 1]; // TODO: 1 of the two accesses can be removed
                for (int i = colCurrent; i < colNext; ++i)
                {
                    fullMatrix[rowIndices[i], colCurrent] = values[i];
                }
            }
            return fullMatrix;
        }

        /// <summary>
        /// Perhaps this should be manually inlined. Testing needed.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ProcessIndices(ref int rowIdx, ref int colIdx)
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
