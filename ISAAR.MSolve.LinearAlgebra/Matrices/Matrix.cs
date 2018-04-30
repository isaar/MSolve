using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.ArrayManipulations;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: align data using mkl_malloc
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// General matrix. Dense (full) storage. Uses MKL. Stored as 1D column major array.
    /// </summary>
    public class Matrix: IMatrix, ISliceable2D
    {
        private readonly double[] data;

        private Matrix(double[] data, int numRows, int numColumns)
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
        /// TODO: make this package-private. It should only be used for passing raw arrays to linear algebra libraries.
        /// </summary>
        internal double[] InternalData { get { return data; } }

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

        public static Matrix CreateIdentity(int order)
        {
            double[] data = new double[order * order];
            for (int j = 0; j < order; ++j) data[j * order + j] = 1.0;
            return new Matrix(data, order, order);
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
        public static Matrix operator +(Matrix matrix1, Matrix matrix2) => matrix1.Axpy(matrix2, 1.0);
        public static Matrix operator -(Matrix matrix1, Matrix matrix2) => matrix1.Axpy(matrix2, -1.0);
        public static Matrix operator *(double scalar, Matrix matrix) => matrix.Scale(scalar);
        public static Matrix operator *(Matrix matrix, double scalar)=> matrix.Scale(scalar);
        public static Matrix operator *(Matrix matrixLeft, Matrix matrixRight)
            => matrixLeft.MultiplyRight(matrixRight, false, false);
        public static Vector operator *(Matrix matrixLeft, Vector vectorRight)
            => matrixLeft.MultiplyRight(vectorRight, false);
        public static Vector operator *(Vector vectorLeft, Matrix matrixRight)
            => matrixRight.MultiplyRight(vectorLeft, true);
        #endregion

        public Matrix AppendBottom(Matrix other)
        {
            Preconditions.CheckSameColDimension(this, other);
            double[] result = ArrayColMajor.JoinVertically(this.NumRows, this.NumColumns, this.data,
                other.NumRows, other.NumColumns, other.data);
            return new Matrix(result, this.NumRows + other.NumColumns, NumColumns);
        }

        public Matrix AppendRight(Matrix other)
        {
            Preconditions.CheckSameRowDimension(this, other);
            double[] result = ArrayColMajor.JoinHorizontally(this.NumRows, this.NumColumns, this.data,
                other.NumRows, other.NumColumns, other.data);
            return new Matrix(result, NumRows, this.NumColumns + other.NumColumns);
        }

        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix casted) return Axpy(casted, otherCoefficient);
            else return otherMatrix.LinearCombination(otherCoefficient, this, 1.0); // To avoid accessing zero entries
        }

        public Matrix Axpy(Matrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref result[0], 1);
            return new Matrix(result, NumRows, NumColumns);
        }

        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix casted) AxpyIntoThis(casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                for (int j = 0; j < NumColumns; ++j)
                {
                    for (int i = 0; i < NumRows; ++i)
                    {
                        this.data[j * NumRows + i] += otherCoefficient * otherMatrix[i, j];
                    }
                }
            }
        }

        public void AxpyIntoThis(Matrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref this.data[0], 1);
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

        public Matrix Copy()
        {
            //TODO: Perhaps this should use BLAS. 
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return new Matrix(clone, NumRows, NumColumns);
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
            if (other is Matrix casted) return DoEntrywise(casted, binaryOperation);
            else return other.DoEntrywise(this, binaryOperation); // To avoid accessing zero entries
        }

        public Matrix DoEntrywise(Matrix other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            var result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(this.data[i], other.data[i]);
            return new Matrix(result, NumRows, NumColumns);
        }

        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is Matrix casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, other);
                for (int j = 0; j < NumColumns; ++j)
                {
                    for (int i = 0; i < NumRows; ++i)
                    {
                        int index1D = j * NumRows + i;
                        this.data[index1D] = binaryOperation(this.data[index1D], other[i, j]);
                    }
                }
            }
        }

        public void DoEntrywiseIntoThis(Matrix other, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, other);
            for (int i = 0; i < data.Length; ++i) this.data[i] = binaryOperation(this.data[i], other.data[i]);
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

        void IMatrix.DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            DoToAllEntriesIntoThis(unaryOperation);
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
            if (other is Matrix casted)
            {
                //Check each dimension, rather than the lengths of the internal buffers
                if (!Preconditions.AreSameMatrixDimensions(this, casted)) return false;
                double[] otherData = casted.data;
                var comparer = new ValueComparer(1e-13);
                for (int i = 0; i < this.data.Length; ++i)
                {
                    if (!comparer.AreEqual(this.data[i], otherData[i])) return false;
                }
                return true;
            }
            else return other.Equals(this, tolerance); // To avoid accessing zero entries
        }

        public CholeskyFull FactorCholesky()
        {
            Preconditions.CheckSquare(this);
            // Copy matrix. This may exceed available memory and needs an extra O(n^2) space. 
            // To avoid these, set "inPlace=true".
            return CholeskyFull.Factorize(NumColumns, CopyInternalData());
        }

        /// <summary>
        /// Calculates the LQ factorization of a matrix with <see cref="NumRows"/> &lt;= <see cref="NumColumns"/>, such that
        /// A = L*Q. Requires an extra 
        /// <see cref="NumRows"/>*<see cref="NumColumns"/> + min(<see cref="NumRows"/>,<see cref="NumColumns"/>) available
        /// memory. 
        /// </summary>
        /// <returns></returns>
        public LQFactorization FactorLQ()
        {
            return LQFactorization.Factorize(NumRows, NumColumns, CopyInternalData());
        }

        public LUFactorization FactorLU()
        {
            Preconditions.CheckSquare(this);
            // Copy matrix. This may exceed available memory and needs an extra O(n^2) space. 
            // To avoid these, set "inPlace=true".
            return LUFactorization.Factorize(NumColumns, CopyInternalData());
        }

        /// <summary>
        /// Calculates the QR factorization of a matrix with <see cref="NumRows"/> &gt;= <see cref="NumColumns"/>, such that
        /// A = Q*R. Requires an extra 
        /// <see cref="NumRows"/>*<see cref="NumColumns"/> + min(<see cref="NumRows"/>,<see cref="NumColumns"/>) available 
        /// memory. 
        /// </summary>
        /// <returns></returns>
        public QRFactorization FactorQR()
        {
            return QRFactorization.Factorize(NumRows, NumColumns, CopyInternalData());
        }

        public Vector GetDiagonal()
        {
            return Vector.CreateFromArray(GetDiagonalAsArray(), false);
        }

        public double[] GetDiagonalAsArray()
        {
            Preconditions.CheckSquare(this);
            double[] diag = new double[NumRows];
            for (int i = 0; i < NumRows; ++i) diag[i] = data[i * NumRows + i];
            return diag;
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
            else return FactorLU().Invert(true);
        }

        public (Matrix inverse, double determinant) InvertAndDetermninant()
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
                return (factor.Invert(true), factor.CalcDeterminant());
            }
        }

        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else return otherMatrix.LinearCombination(otherCoefficient, this, thisCoefficient); // To avoid accessing zero entries
        }

        public Matrix LinearCombination(double thisCoefficient, Matrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref result[0], 1);
            return new Matrix(result, NumRows, NumColumns);
        }

        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                for (int j = 0; j < NumColumns; ++j)
                {
                    for (int i = 0; i < NumRows; ++i)
                    {
                        int index1D = j * NumRows + i;
                        this.data[index1D] = thisCoefficient * this.data[index1D] + otherCoefficient * otherMatrix[i, j];
                    }
                }
            }
        }

        public void LinearCombinationIntoThis(double thisCoefficient, Matrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref this.data[0], 1);
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

        public Vector MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            if (vector is Vector casted) return MultiplyRight(casted, transposeThis);
            else throw new NotImplementedException();
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
            return Vector.CreateFromArray(result, false);
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            for (int i = 0; i < data.Length; ++i) aggregator = processEntry(data[i], aggregator);
            // no zeros implied
            return finalize(aggregator);
        }

        //TODO: perhaps I should transfer this to a permutation matrix (implemented as a vector).
        public Matrix Reorder(IReadOnlyList<int> permutation, bool oldToNew)
        {
            Preconditions.CheckSquare(this);
            if (permutation.Count != NumRows) throw new NonMatchingDimensionsException(
                $"This matrix has order = {NumRows}, while the permutation vector has {permutation.Count} entries.");
            if (oldToNew) return new Matrix(ArrayColMajor.ReorderOldToNew(NumRows, data, permutation), NumRows, NumRows);
            else return new Matrix(ArrayColMajor.ReorderNewToOld(NumRows, data, permutation), NumRows, NumRows);
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

        public void SetColumn(int colIdx, Vector colValues)
        {
            Preconditions.CheckIndexCol(this, colIdx);
            Preconditions.CheckSameRowDimension(this, colValues);
            ArrayColMajor.SetCol(NumRows, NumColumns, data, colIdx, colValues.InternalData);
        }

        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            data[colIdx * NumRows + rowIdx] = value;
        }

        public void SetRow(int rowIdx, Vector rowValues)
        {
            Preconditions.CheckIndexRow(this, rowIdx);
            Preconditions.CheckSameColDimension(this, rowValues);
            ArrayColMajor.SetRow(NumRows, NumColumns, data, rowIdx, rowValues.InternalData);
        }

        public void SetSubmatrix(int rowStart, int colStart, Matrix submatrix)
        {
            Preconditions.CheckIndices(this, rowStart, colStart);
            if ((rowStart + submatrix.NumRows > this.NumRows) || (colStart + submatrix.NumColumns > this.NumColumns))
            {
                throw new NonMatchingDimensionsException("The submatrix doesn't fit inside this matrix, at least when starting"
                    + " from the specified entry.");
            }
            ArrayColMajor.SetSubmatrix(this.NumRows, this.NumColumns, this.data, rowStart, colStart, 
                submatrix.NumRows, submatrix.NumColumns, submatrix.data);
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

        public Vector SliceColumn(int colIndex)
        {
            double[] result = new double[NumRows];
            Array.Copy(data, colIndex * NumRows, result, 0, NumRows);
            return Vector.CreateFromArray(data, false);
        }

        public Vector SliceRow(int rowIndex)
        {
            double[] result = new double[NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                result[j] = data[j * NumRows + rowIndex];
            }
            return Vector.CreateFromArray(result, false);
        }

        public void SVD(double[] w, double[,] v)
        {
            DenseStrategies.SVD(this, w, v);
        }

        /// <summary>
        /// Doesn't copy anything. Remove this once the design is cleaned. 
        /// </summary>
        /// <returns></returns>
        public Numerical.LinearAlgebra.Interfaces.IMatrix2D ToLegacyMatrix()
        {
            return new Numerical.LinearAlgebra.Matrix2D(CopyToArray2D());
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

        private double[] CopyInternalData()
        {
            double[] dataCopy = new double[data.Length];
            Array.Copy(data, dataCopy, data.Length);
            return dataCopy;
        }
    }
}
