using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// Compressed Sparse Columns format.
    /// </summary>
    public class SymmetricCSC: IIndexable2D
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
        public SymmetricCSC(double[] values, int[] rowIndices, int[] colOffsets, bool checkInput)
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
            if (other is SymmetricCSC) // In case both matrices have the exact same index arrays
            {
                SymmetricCSC otherCSC = (SymmetricCSC)other;
                if (HaveSameIndexArrays(otherCSC))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherCSC.values[i]);
                    }
                    return new SymmetricCSC(resultValues, rowIndices, colOffsets, false);
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

        public void DoPointwiseIntoThis(SymmetricCSC other, Func<double, double, double> binaryOperation)
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
                return new SymmetricCSC(newValues, cloneRowIndices, cloneColOffsets, false);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new SymmetricCSC(newValues, rowIndices, colOffsets, false).ToFullMatrix();
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

        public bool Equals(IIndexable2D other, ValueComparer comparer = null)
        {
            if (!Preconditions.AreSameMatrixDimensions(this, other)) return false;
            if (comparer == null) comparer = new ValueComparer(1e-13);
            // All entries must be checked. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            for (int j = 0; j < this.NumColumns; ++j)
            {
                for (int i = 0; i < this.NumRows; ++i)
                {
                    if (!comparer.AreEqual(this[i, j], other[i, j])) return false;
                }
            }
            return true;
        }

        public SuiteSparseCholesky FactorCholesky()
        {
            return SuiteSparseCholesky.CalcFactorization(NumColumns, NumNonZerosUpper, values, rowIndices, colOffsets);
        }

        /// <summary>
        /// 2 or more sparse matrices may have the same (2 references to 1 object) index arrays. 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool HaveSameIndexArrays(SymmetricCSC other)
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
        public VectorMKL MultiplyRight(VectorMKL vector, bool transposeThis = false)
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

        // To write the full matrix use the extension IIndexable.WriteToConsole()
        public void WriteToConsole()
        {
            Console.WriteLine("Values:");
            for (int i = 0; i < values.Length; ++i)
            {
                Console.Write(values[i]);
                Console.Write(" ");
            }
            Console.WriteLine("Row indices:");
            for (int i = 0; i < rowIndices.Length; ++i)
            {
                Console.Write(rowIndices[i]);
                Console.Write(" ");
            }
            Console.WriteLine("Column offsets:");
            for (int i = 0; i < colOffsets.Length; ++i)
            {
                Console.Write(colOffsets[i]);
                Console.Write(" ");
            }
        }

        public void WriteToFile(string path, bool append = false, Array2DFormatting format = null)
        {
            using (var writer = new StreamWriter(path, append))
            {
                writer.WriteLine("Values:");
                for (int i = 0; i < values.Length; ++i)
                {
                    writer.Write(values[i]);
                    writer.Write(" ");
                }
                writer.WriteLine("Row indices:");
                for (int i = 0; i < rowIndices.Length; ++i)
                {
                    writer.Write(rowIndices[i]);
                    writer.Write(" ");
                }
                writer.WriteLine("Column offsets:");
                for (int i = 0; i < colOffsets.Length; ++i)
                {
                    writer.Write(colOffsets[i]);
                    writer.Write(" ");
                }

#if DEBUG
                writer.Flush(); // If the user inspects the file while debugging, make sure the contentss are written.
#endif
            }

            
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
