using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

//TODO: align data using mkl_malloc
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// General matrix. Dense (full) storage. Uses MKL. Stored as 1D column major array.
    /// </summary>
    public class Matrix: IMatrixView, ISliceable2D
    {
        protected readonly double[] data;

        protected Matrix(double[] data, int numRows, int numColumns)
        {
            this.data = data;
            this.NumRows = numRows;
            this.NumColumns = numColumns;
        }

        public bool IsSquare { get { return NumRows == NumColumns; } }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// Only structural non zeros
        /// </summary>
        /// <returns></returns>
        public int NumNonZeros { get { return NumRows * NumColumns; } }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// The entry with row index = i and column index = j. 
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= i &lt; <see cref="NumRows"/></param>
        /// <param name="colIdx">The column index: 0 &lt;= j &lt; <see cref="NumColumns"/></param>
        /// <returns>The entry with indices i, j</returns>
        public double this[int rowIdx, int colIdx] //TODO: Should I add bound checking?
        {
            get { return data[colIdx * NumRows + rowIdx]; }
            set { data[colIdx * NumRows + rowIdx] = value; }
        }

        /// <summary>
        /// Create a new <see cref="Matrix"/> from a provided array. The array will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional array containing the elements of the matrix</param>
        /// <returns></returns>
        public static Matrix CreateFromArray(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numCols = array2D.GetLength(1); 
            return new Matrix(Conversions.Array2DToFullColMajor(array2D), numRows, numCols);
        }


        /// <summary>
        /// Create a new <see cref="Matrix"/> from a provided array. The array will can be copied for extra safety or not for 
        /// extra performance.
        /// </summary>
        /// <param name="array1D">A 1-dimensional array containing the elements of the matrix in column major order. Its length 
        /// must be equal to <see cref="numRows"/> + <see cref="NumColumns"/>. It will not be checked.</param>
        /// <param name="numRows">The number of rows of the matrix</param>
        /// <param name="numColumns">The number of columns of the matrix</param>
        /// <param name="copyArray">True to make a deep copy of <see cref="array1D"/>. 
        /// False (default) to use <see cref="array1D"/> as its internal storage.</param>
        /// <returns></returns>
        public static Matrix CreateFromArray(double[] array1D, int numRows, int numColumns, bool copyArray = false)
        {
            if (copyArray)
            {
                var clone = new double[array1D.Length];
                Array.Copy(array1D, clone, clone.Length);
                return new Matrix(clone, numRows, numColumns);
            }
            else
            {
                return new Matrix(array1D, numRows, numColumns);
            }
        }

        /// <summary>
        /// The original matrix will be copied.
        /// </summary>
        /// <param name="original"></param>
        /// <returns></returns>
        public static Matrix CreateFromMatrix(Matrix original) 
        {
            //TODO: Perhaps this should use BLAS. 
            //TODO: Perhaps it should be an instance method CopyToMatrix(). Or the instance method would return an interface.
            double[] data = original.data;
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return new Matrix(clone, original.NumRows, original.NumColumns);
        }

        public static Matrix CreateWithValue(int numRows, int numColumns, double value)
        {
            double[] data = new double[numRows * numColumns];
            for (int i = 0; i < data.Length; ++i) data[i] = value;
            return new Matrix(data, numRows, numColumns);
        }

        /// <summary>
        /// Create a new <see cref="Matrix"/> with the specified dimensions and all entries equal to 0.
        /// </summary> 
        /// <param name="numRows">The number of rows of the matrix.</param>
        /// <param name="numColumns">The number of rows of the matrix.</param>
        /// <returns></returns>
        public static Matrix CreateZero(int numRows, int numColumns)
        {
            double[] data = new double[numRows * numColumns];
            return new Matrix(data, numRows, numColumns);
        }

        #region operators (use extension operators when they become available)
        public static Matrix operator +(Matrix matrix1, Matrix matrix2) => matrix1.Axpy(1.0, matrix2);
        public static Matrix operator -(Matrix matrix1, Matrix matrix2) => matrix1.Axpy(-1.0, matrix2);
        public static Matrix operator *(double scalar, Matrix matrix) => matrix.Scale(scalar);
        public static Matrix operator *(Matrix matrix, double scalar)=> matrix.Scale(scalar);
        public static Matrix operator *(Matrix matrixLeft, Matrix matrixRight)
            => matrixLeft.MultiplyRight(matrixRight, false, false);
        public static VectorMKL operator *(Matrix matrixLeft, VectorMKL vectorRight)
            => matrixLeft.MultiplyRight(vectorRight, false);
        public static VectorMKL operator *(VectorMKL vectorRight, Matrix matrixLeft)
            => matrixLeft.MultiplyRight(vectorRight, true);
        #endregion

        /// <summary>
        /// result = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        /// <returns></returns>
        public Matrix Axpy(double scalar, Matrix other)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[this.data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpy(this.data.Length, scalar, ref other.data[0], 1, ref result[0], 1);
            return new Matrix(result, NumRows, NumColumns);
        }

        /// <summary>
        /// this = this + scalar * other
        /// </summary>
        /// <param name="other"></param>
        /// <param name="scalar"></param>
        public void AxpyIntoThis(double scalar, Matrix other)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            CBlas.Daxpy(this.data.Length, scalar, ref other.data[0], 1, ref this.data[0], 1);
        }

        public double CalcDeterminant()
        {
            if ((NumRows == 2) && (NumColumns == 2))
            {
                return AnalyticFormulas.Matrix2x2ColMajorDeterminant(data);
            }
            else if ((NumRows == 3) && (NumColumns == 3))
            {
                return AnalyticFormulas.Matrix3x3ColMajorDeterminant(data);
            }
            else return FactorLU().CalcDeterminant();
        }

        /// <summary>
        /// Copy the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
        /// and length(1) = <see cref="NumColumns"/>. 
        /// </summary>
        /// <returns>A new <see cref="double"/>[<see cref="NumRows"/>, <see cref="NumRows"/>] array 
        /// with the entries of the matrix</returns>
        public double[,] CopyToArray2D()
        {
            return Conversions.FullColMajorToArray2D(data, NumRows, NumColumns);
        }

        public IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is Matrix) return DoEntrywise((Matrix)other, binaryOperation);
            else return other.DoEntrywise(this, binaryOperation); // To avoid accessing zero entries
        }

        public Matrix DoEntrywise(Matrix other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            var result = new double[NumRows * NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                for (int i = 0; i < NumRows; ++i)
                {
                    int idx = j * NumRows + NumColumns;
                    result[idx] = binaryOperation(this.data[idx], other.data[idx]);
                }
            }
            return new Matrix(result, NumRows, NumColumns);
        }

        public void DoEntrywiseIntoThis(Matrix other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            for (int j = 0; j < NumColumns; ++j)
            {
                for (int i = 0; i < NumRows; ++i)
                {
                    int idx = j * NumRows + NumColumns;
                    this.data[idx] = binaryOperation(this.data[idx], other.data[idx]);
                }
            }
        }

        IMatrixView IMatrixView.DoToAllEntries(Func<double, double> unaryOperation)
        {
            return DoToAllEntries(unaryOperation);
        }

        public Matrix DoToAllEntries(Func<double, double> unaryOperation)
        {
            var result = new double[NumRows * NumColumns];
            for (int i = 0; i < NumRows * NumColumns; ++i)
            {
                result[i] = unaryOperation(data[i]);
            }
            return new Matrix(result, NumRows, NumColumns);
        }

        // Ok for a DenseMatrix, but for sparse formats some operation (e.g scale) maintain the sparsity pattern,
        // while others don't
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            for (int i = 0; i < NumRows * NumColumns; ++i)
            {
                data[i] = unaryOperation(data[i]);
            }
        }

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            if (other is Matrix)
            {
                //Check each dimension, rather than the lengths of the internal buffers
                if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
                double[] otherData = ((Matrix)other).data;
                var comparer = new ValueComparer(1e-13);
                for (int i = 0; i < this.data.Length; ++i)
                {
                    if (!comparer.AreEqual(this.data[i], otherData[i])) return false;
                }
                return true;
            }
            else return other.Equals(this, tolerance); // To avoid accessing zero entries
        }

        public LUFactorization FactorLU()
        {
            if (IsSquare) return LUFactorization.CalcFactorization(NumRows, data);
            else throw new NonMatchingDimensionsException($"The matrix must be square, but was {NumRows}x{NumColumns}");
        }

        public Matrix Invert()
        {
            if ((NumRows == 2) && (NumColumns == 2))
            {
                (double[] inverse, double det) = AnalyticFormulas.Matrix2x2ColMajorInvert(data);
                return new Matrix(inverse, 2, 2);
            }
            else if ((NumRows == 3) && (NumColumns == 3))
            {
                (double[] inverse, double det) = AnalyticFormulas.Matrix3x3ColMajorInvert(data);
                return new Matrix(inverse, 3, 3);
            }
            else return FactorLU().InvertInPlace();
        }

        public (Matrix inverse, double determinant) InvertAndDetermninat()
        {
            if ((NumRows == 2) && (NumColumns == 2))
            {
                (double[] inverse, double det) = AnalyticFormulas.Matrix2x2ColMajorInvert(data);
                return (new Matrix(inverse, 2, 2), det);
            }
            else if ((NumRows == 3) && (NumColumns == 3))
            {
                (double[] inverse, double det) = AnalyticFormulas.Matrix3x3ColMajorInvert(data);
                return (new Matrix(inverse, 3, 3), det);
            }
            else
            {
                LUFactorization factor = FactorLU();
                return (factor.InvertInPlace(), factor.CalcDeterminant());
            }
                
        }

        /// <summary>
        /// result = thisScalar * this + otherScalar * otherMatrix
        /// </summary>
        /// <returns></returns>
        public Matrix LinearCombination(double thisScalar, double otherScalar, Matrix otherMatrix)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[this.data.Length];
            Array.Copy(this.data, result, this.data.Length);
            CBlas.Daxpby(this.data.Length, otherScalar, ref otherMatrix.data[0], 1, thisScalar, ref result[0], 1);
            return new Matrix(result, this.NumRows, this.NumColumns);
        }

        /// <summary>
        /// this = this + scalar * otherMatrix
        /// </summary>
        /// <returns></returns>
        public void LinearCombinationIntoThis(double thisScalar, double otherScalar, Matrix otherMatrix)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpby(this.data.Length, otherScalar, ref otherMatrix.data[0], 1, thisScalar, ref this.data[0], 1);
        }

        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return other.MultiplyRight(this, transposeOther, transposeThis);
        }

        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            if (other is Matrix) return MultiplyRight((Matrix)other, transposeThis);
            else return other.MultiplyLeft(this, transposeOther, transposeThis);
        }

        /// <summary>
        /// Matrix-matrix multiplication, with the other matrix on the right: this [m x k] * other [k x n] 
        /// or transpose(this [k x m]) * other [k x n].
        /// </summary>
        /// <param name="other">A matrix with as many rows as the column of this matrix.</param>
        /// <param name="transposeThis">Set to true to transpose this (the left matrix). Unless the transpose matrix is used in 
        ///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
        /// <param name="transposeOther">Set to true to transpose other (the right matrix). Unless the transpose matrix is used in 
        ///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
        /// <returns>A matrix with dimensions (m x n)</returns>
        public Matrix MultiplyRight(Matrix other, bool transposeThis = false, bool transposeOther = false)
        {
            int leftRows, leftCols, rightRows, rightCols;
            CBLAS_TRANSPOSE transposeLeft, transposeRight;
            if (transposeThis)
            {
                transposeLeft = CBLAS_TRANSPOSE.CblasTrans;
                leftRows = this.NumColumns;
                leftCols = this.NumRows;
            }
            else
            {
                transposeLeft = CBLAS_TRANSPOSE.CblasNoTrans;
                leftRows = this.NumRows;
                leftCols = this.NumColumns;
            }
            if (transposeOther)
            {
                transposeRight = CBLAS_TRANSPOSE.CblasTrans;
                rightRows = other.NumColumns;
                rightCols = other.NumRows;
            }
            else
            {
                transposeRight = CBLAS_TRANSPOSE.CblasNoTrans;
                rightRows = other.NumRows;
                rightCols = other.NumColumns;
            }

            Preconditions.CheckMultiplicationDimensions(leftCols, rightRows);
            double[] result = new double[leftRows * rightCols];
            CBlas.Dgemm(CBLAS_LAYOUT.CblasColMajor, transposeLeft, transposeRight,
                leftRows, rightCols, leftCols,
                1.0, ref this.data[0], this.NumRows, 
                ref other.data[0], other.NumRows,
                1.0, ref result[0], leftRows);
            return new Matrix(result, leftRows, rightCols);
        }

        public VectorMKL MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            if (vector is VectorMKL) return MultiplyRight((VectorMKL)vector, transposeThis);
            else throw new NotImplementedException();
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
            int leftRows, leftCols;
            CBLAS_TRANSPOSE transpose;
            if (transposeThis)
            {
                transpose = CBLAS_TRANSPOSE.CblasTrans;
                leftRows = NumColumns;
                leftCols = NumRows;
            }
            else
            {
                transpose = CBLAS_TRANSPOSE.CblasNoTrans;
                leftRows = NumRows;
                leftCols = NumColumns;
            }

            Preconditions.CheckMultiplicationDimensions(leftCols, vector.Length);
            double[] result = new double[leftRows];
            CBlas.Dgemv(CBLAS_LAYOUT.CblasColMajor, transpose, NumRows, NumColumns,
                1.0, ref data[0], NumRows,
                ref vector.InternalData[0], 1,
                0.0, ref result[0], 1);
            return VectorMKL.CreateFromArray(result, false);
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            for (int i = 0; i < data.Length; ++i) aggregator = processEntry(data[i], aggregator);
            // no zeros implied
            return finalize(aggregator);
        }

        /// <summary>
        /// result = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public Matrix Scale(double scalar)
        {
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Dscal(data.Length, scalar, ref result[0], 1);
            return new Matrix(result, NumRows, NumColumns);
        }

        /// <summary>
        /// this = scalar * this
        /// </summary>
        /// <param name="scalar"></param>
        public void ScaleIntoThis(double scalar)
        {
            CBlas.Dscal(data.Length, scalar, ref data[0], 1);
        }

        public void SetAll(double value)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = value;
        }

        /// <summary>
        /// Returns a subvector containing only the entries at the provided row and column indices
        /// </summary>
        /// <param name="rowIndices">Rows of the entries to be returned. They must be 0 &lt; = i &lt; 
        ///     <see cref="NumRows"/>.</param>
        /// <param name="colIndices">Columns of the entries to be returned. They must be 0 &lt; = i &lt; 
        ///     <see cref="NumRows"/>.</param>
        /// <returns></returns>
        public Matrix Slice(int[] rowIndices, int[] colIndices)
        {
            double[] submatrix = new double[colIndices.Length * rowIndices.Length];
            int idxCounter = -1;
            foreach (var j in colIndices)
            {
                foreach (var i in rowIndices)
                {
                    submatrix[++idxCounter] = data[j * NumRows + i];
                }
            }
            return new Matrix(submatrix, rowIndices.Length, colIndices.Length);
        }

        /// <summary>
        /// Returns a subvector containing the entries at the indices between the provided start (inclusive) and end (exclusive).
        /// </summary>
        /// <param name="rowStartInclusive">The first row from which to copy entries.</param>
        /// <param name="rowEndExclusive">The row after the last one until which to copy entries.</param>
        /// <param name="colStartInclusive">The first column from which to copy entries.</param>
        /// <param name="colEndExclusive">The column after the last one until which to copy entries.</param>
        /// <returns></returns>
        public Matrix Slice(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
        {
            int newNumRows = rowEndExclusive - rowStartInclusive;
            int newNumCols = colEndExclusive - colStartInclusive;
            double[] submatrix = new double[newNumCols * newNumRows];
            int idxCounter = -1;
            for (int j = colStartInclusive; j < colEndExclusive; ++j)
            {
                for (int i = rowStartInclusive; i < rowEndExclusive; ++i)
                {
                    submatrix[++idxCounter] = data[j * NumRows + i];
                }
            }
            return new Matrix(submatrix, newNumRows, newNumCols);
        }

        public VectorMKL SliceColumn(int colIndex)
        {
            double[] result = new double[NumRows];
            Array.Copy(data, colIndex * NumRows, result, 0, NumRows);
            return VectorMKL.CreateFromArray(data, false);
        }

        public VectorMKL SliceRow(int rowIndex)
        {
            double[] result = new double[NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                result[j] = data[j * NumRows + rowIndex];
            }
            return VectorMKL.CreateFromArray(result, false);
        }

        /// <summary>
        /// Doesn't copy anything. Remove this once the design is cleaned. 
        /// </summary>
        /// <returns></returns>
        public IMatrix2D ToLegacyMatrix()
        {
            return new Matrix2D(CopyToArray2D());
        }

        IMatrixView IMatrixView.Transpose()
        {
            return Transpose();
        }

        public Matrix Transpose()
        {
            //TODO: The wrapper library does not include MKL's blas-like extensions yet. Create my own wrapper or 
            // piggyback on another BLAS function.
            double[] transpose = Conversions.ColumnMajorToRowMajor(data, NumRows, NumColumns);
            return new Matrix(transpose, NumColumns, NumRows);
        }

        public void TransposeIntoThis()
        {
            throw new NotImplementedException("Use mkl_dimatcopy");
        }
    }
}
