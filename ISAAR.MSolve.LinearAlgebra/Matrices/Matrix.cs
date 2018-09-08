using System;
using System.Collections.Generic;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Factorizations;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: align data using mkl_malloc
//TODO: add inplace option for factorizations and leave all subsequent operations (determinant, system solution, etc.) to them
//TODO: remove legacy matrix conversions
//TODO: SetSubrow, SetSubcolumn, SetSubmatrix only need to check the stricter upper bounds.
//TODO: Se https://software.intel.com/en-us/mkl-developer-reference-c-lapmr, 
//      https://software.intel.com/en-us/mkl-developer-reference-c-laswp, 
//      https://software.intel.com/en-us/mkl-developer-reference-c-syswapr for reordering
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// General purpose matrix class. All entries are stored in an 1D column major array. Uses MKL for most operations. 
    /// Authors: Serafeim Bakalakos
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

        /// <summary>
        /// Returns true if <see cref="NumRows"/> == <see cref="NumColumns"/>.
        /// </summary>
        public bool IsSquare { get { return NumRows == NumColumns; } }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of non-zero and explicitly stored zero entries, which is the number of all entries in this 
        /// <see cref="Matrix"/>.
        /// </summary>
        public int NumNonZeros { get { return NumRows * NumColumns; } }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// The internal array that stores the entries of the matrix. It should only be used for passing the raw array to linear 
        /// algebra libraries.
        /// </summary>
        internal double[] InternalData { get { return data; } }

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>.
        /// </summary>
        /// <remarks>
        /// Also note that it may be possible to pass in <paramref name="rowIdx"/> &gt;= <see cref="IIndexable2D.NumRows"/> or
        /// <paramref name="rowIdx"/> &lt; 0, without throwing <see cref="IndexOutOfRangeException"/>, since the indices are not  
        /// checked explicitly. The constraints on <paramref name="colIdx"/> described in the interfaces will correctly throw
        /// <see cref="IndexOutOfRangeException"/> if violated.
        /// </remarks>
        public double this[int rowIdx, int colIdx] //TODO: Should I add bound checking?
        {
            get { return data[colIdx * NumRows + rowIdx]; }
            set { data[colIdx * NumRows + rowIdx] = value; }
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix"/> by copying the entries of <paramref name="array2D"/>.
        /// </summary>
        /// <param name="array2D">A 2-dimensional array containing the entries of the matrix. It will be copied.</param>
        public static Matrix CreateFromArray(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numCols = array2D.GetLength(1); 
            return new Matrix(Conversions.Array2DToFullColMajor(array2D), numRows, numCols);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix"/> with <paramref name="array1D"/> or a clone as its internal array.
        /// </summary>
        /// <param name="array1D">A 1-dimensional array containing the elements of the matrix in column major order. Its length 
        ///     must be equal to <see cref="numRows"/> * <see cref="NumColumns"/>. It will not be checked.</param>
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numColumns">The number of columns of the new matrix.</param>
        /// <param name="copyArray">If true, <paramref name="array1D"/> will be copied and the new <see cref="Matrix"/> instance 
        ///     will have a reference to the copy, which is safer. If false, the new matrix will have a reference to 
        ///     <paramref name="array1D"/> itself, which is faster.</param>
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
        /// Initializes a new instance of <see cref="Matrix"/> by copying the entries of <paramref name="matrix"/>.
        /// </summary>
        /// <param name="matrix">The entries of the legacy matrix instance <see cref="Numerical.LinearAlgebra.Matrix2D"/>.</param>
        public static Matrix CreateFromLegacyMatrix(Numerical.LinearAlgebra.Matrix2D matrix)
        {
            // The other matrix might be transposed internally.
            bool isTransposed = (matrix.Data[0, 1] != matrix[0, 1]);
            double[] data = isTransposed ? 
                Conversions.Array2DToFullRowMajor(matrix.Data) : Conversions.Array2DToFullColMajor(matrix.Data);
            return new Matrix(data, matrix.Rows, matrix.Columns);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix"/> that is equal to the identity matrix, namely a square matrix with 
        /// non-diagonal entries being equal to 0 and diagonal entries being equal to 1.
        /// </summary>
        /// <param name="order">The number of rows/columns of the identity matrix.</param>
        public static Matrix CreateIdentity(int order)
        {
            double[] data = new double[order * order];
            for (int j = 0; j < order; ++j) data[j * order + j] = 1.0;
            return new Matrix(data, order, order);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix"/> with all entries being equal to <paramref name="value"/>.
        /// </summary>
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numColumns">The number of columns of the new matrix.</param>
        /// <param name="value">The value that all entries of the new matrix will be initialized to.</param>
        public static Matrix CreateWithValue(int numRows, int numColumns, double value)
        {
            double[] data = new double[numRows * numColumns];
            for (int i = 0; i < data.Length; ++i) data[i] = value;
            return new Matrix(data, numRows, numColumns);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="Matrix"/> with all entries being equal to 0.
        /// </summary> 
        /// <param name="numRows">The number of rows of the new matrix.</param>
        /// <param name="numColumns">The number of rows of the new matrix.</param>
        /// <returns></returns>
        public static Matrix CreateZero(int numRows, int numColumns)
        {
            double[] data = new double[numRows * numColumns];
            return new Matrix(data, numRows, numColumns);
        }

        #region operators (use extension operators when they become available)
        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="matrix1"/>[i, j] + <paramref name="matrix2"/>[i, j], 
        /// for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>.
        /// The resulting entries are written to a new <see cref="Matrix"/> instance.
        /// </summary>
        /// <param name="matrix1">The first <see cref="Matrix"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix2"/>.</param>
        /// <param name="matrix2">The second <see cref="Matrix"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix1"/> and <paramref name="matrix2"/>
        ///     have a different number of <see cref="NumRows"/> or <see cref="NumColumns"/>.</exception>
        public static Matrix operator +(Matrix matrix1, Matrix matrix2) => matrix1.Axpy(matrix2, 1.0);

        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="matrix1"/>[i, j] - <paramref name="matrix2"/>[i, j], 
        /// for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>.
        /// The resulting entries are written to a new <see cref="Matrix"/> instance.
        /// </summary>
        /// <param name="matrix1">The first <see cref="Matrix"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix2"/>.</param>
        /// <param name="matrix2">The second <see cref="Matrix"/> operand. It must have as many rows and columns as 
        ///     <paramref name="matrix1"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix1"/> and <paramref name="matrix2"/>
        ///     have a different number of <see cref="NumRows"/> or <see cref="NumColumns"/>.</exception>
        public static Matrix operator -(Matrix matrix1, Matrix matrix2) => matrix1.Axpy(matrix2, -1.0);

        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="scalar"/> * <paramref name="matrix1"/>[i, j],
        /// for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>.
        /// The resulting entries are written to a new <see cref="Matrix"/> instance.
        /// </summary>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        /// <param name="matrix">The matrix to multiply.</param>
        public static Matrix operator *(double scalar, Matrix matrix) => matrix.Scale(scalar);

        /// <summary>
        /// Performs the operation: result[i, j] = <paramref name="scalar"/> * <paramref name="matrix1"/>[i, j],
        /// for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>.
        /// The resulting entries are written to a new <see cref="Matrix"/> instance.
        /// </summary>
        /// <param name="matrix">The matrix to multiply.</param>
        /// <param name="scalar">The scalar value that will be multiplied with all vector entries.</param>
        public static Matrix operator *(Matrix matrix, double scalar)=> matrix.Scale(scalar);

        /// <summary>
        /// Performs the matrix-matrix multiplication: result = <paramref name="matrixLeft"/> * <paramref name="matrixRight"/>.
        /// If <paramref name="matrixLeft"/> is m1-by-n1 and <paramref name="matrixRight"/> is m2-by-n2, then n1 must be equal to
        /// m2. The result will be an m1-by-n2 matrix, written to a new <see cref="Matrix"/> instance.
        /// </summary>
        /// <param name="matrixLeft">The <see cref="Matrix"/> operand on the left.</param>
        /// <param name="matrixRight">The <see cref="Matrix"/> operand on the right.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixLeft"/>.<see cref="NumColumns"/> is 
        ///     different than <paramref name="matrixRight"/>.<see cref="NumRows"/>.</exception>
        public static Matrix operator *(Matrix matrixLeft, Matrix matrixRight)
            => matrixLeft.MultiplyRight(matrixRight, false, false);

        /// <summary>
        /// Performs the matrix-vector multiplication: result = <paramref name="matrixLeft"/> * <paramref name="vectorRight"/>.
        /// If <paramref name="matrixLeft"/> is m1-by-n1 and <paramref name="vectorRight"/> has length = n2, then n1 must be 
        /// equal to n2. The result will be a vector with length = m1, written to a new <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="matrixLeft">The <see cref="Matrix"/> operand on the left.</param>
        /// <param name="vectorRight">The <see cref="Vector"/> operand on the right. It can be considered as a column 
        ///     vector.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixLeft"/>.<see cref="NumColumns"/> is 
        ///     different than <paramref name="vectorRight"/>.<see cref="Vector.Length"/>.</exception>
        public static Vector operator *(Matrix matrixLeft, Vector vectorRight)
            => matrixLeft.MultiplyRight(vectorRight, false);

        /// <summary>
        /// Performs the matrix-vector multiplication: result = <paramref name="vectorLeft"/> * <paramref name="matrixRight"/>.
        /// If <paramref name="matrixRight"/> is m1-by-n1 and <paramref name="vectorLeft"/> has length = n2, then m1 must be 
        /// equal to n2. The result will be a vector with length = n1, written to a new <see cref="Vector"/> instance.
        /// </summary>
        /// <param name="vectorLeft">The <see cref="Vector"/> operand on the left. It can be considered as a row vector.</param>
        /// <param name="matrixRight">The <see cref="Matrix"/> operand on the right.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrixRight"/>.<see cref="NumRows"/> is 
        ///     different than <paramref name="vectorLeft"/>.<see cref="Vector.Length"/>.</exception>
        public static Vector operator *(Vector vectorLeft, Matrix matrixRight)
            => matrixRight.MultiplyRight(vectorLeft, true);
        #endregion

        /// <summary>
        /// Creates a new <see cref="Matrix"/> that contains all rows of this <see cref="Matrix"/> instance, followed by all rows 
        /// of <paramref name="matrix"/>. If this is m1-by-n1 and <paramref name="matrix"/> is m2-by-n2, then n2 must be equal 
        /// to n1 and the resulting matrix will be (m1+m2)-by-n1.
        /// </summary>
        /// <param name="matrix">The matrix whose rows will be appended after all rows of this <see cref="Matrix"/> 
        ///     instance.</param>
        public Matrix AppendBottom(Matrix matrix)
        {
            Preconditions.CheckSameColDimension(this, matrix);
            double[] result = ArrayColMajor.JoinVertically(this.NumRows, this.NumColumns, this.data,
                matrix.NumRows, matrix.NumColumns, matrix.data);
            return new Matrix(result, this.NumRows + matrix.NumColumns, NumColumns);
        }

        /// <summary>
        /// Creates a new <see cref="Matrix"/> that contains all columns of this <see cref="Matrix"/> instance, followed by all 
        /// columns of <paramref name="matrix"/>. If this is m1-by-n1 and <paramref name="matrix"/> is m2-by-n2, then m2 must be 
        /// equal to mn1 and the resulting matrix will be m1-by-(n1+n2).
        /// </summary>
        /// <param name="matrix">The matrix whose columns will be appended after all columns of this <see cref="Matrix"/> 
        ///     instance.</param>
        public Matrix AppendRight(Matrix matrix)
        {
            Preconditions.CheckSameRowDimension(this, matrix);
            double[] result = ArrayColMajor.JoinHorizontally(this.NumRows, this.NumColumns, this.data,
                matrix.NumRows, matrix.NumColumns, matrix.data);
            return new Matrix(result, NumRows, this.NumColumns + matrix.NumColumns);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
        /// </summary>
        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix casted) return Axpy(casted, otherCoefficient);
            else return otherMatrix.LinearCombination(otherCoefficient, this, 1.0); // To avoid accessing zero entries
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix is written to a new <see cref="Matrix"/> and then returned.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
        ///     <see cref="Matrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
        public Matrix Axpy(Matrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref result[0], 1);
            return new Matrix(result, NumRows, NumColumns);
        }

        /// <summary>
        /// See <see cref="IMatrix.AxpyIntoThis(IMatrixView, double)"/>.
        /// </summary>
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

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="Matrix"/> instance.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
        ///     <see cref="Matrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
        public void AxpyIntoThis(Matrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref this.data[0], 1);
        }

        /// <summary>
        /// Calculates the determinant of this matrix, which must be square. If the inverse matrix is also needed, use
        /// <see cref="InvertAndDetermninant"/> instead.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if this matrix is not square.</exception>
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
        /// Initializes a new instance of <see cref="Matrix"/> by copying the entries of this instance.
        /// </summary>
        public Matrix Copy()
        {
            //TODO: Perhaps this should use BLAS. 
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return new Matrix(clone, NumRows, NumColumns);
        }

        /// <summary>
        /// Copies the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
        /// and length(1) = <see cref="NumColumns"/>. 
        /// </summary>
        public double[,] CopyToArray2D()
        {
            return Conversions.FullColMajorToArray2D(data, NumRows, NumColumns);
        }

        /// <summary>
        /// See <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public IMatrixView DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation)
        {
            if (matrix is Matrix casted) return DoEntrywise(casted, binaryOperation);
            else return matrix.DoEntrywise(this, binaryOperation); // To avoid accessing zero entries
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// result[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i, j]). 
        /// The resulting matrix is written to a new <see cref="Matrix"/> and then returned.
        /// </summary>
        /// <param name="matrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
        ///     <see cref="Matrix"/> instance.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix"/> has different 
        ///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
        public Matrix DoEntrywise(Matrix matrix, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, matrix);
            var result = new double[data.Length];
            for (int i = 0; i < data.Length; ++i) result[i] = binaryOperation(this.data[i], matrix.data[i]);
            return new Matrix(result, NumRows, NumColumns);
        }

        /// <summary>
        /// See <see cref="IMatrix.DoEntrywiseIntoThis(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IMatrixView matrix, Func<double, double, double> binaryOperation)
        {
            if (matrix is Matrix casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, matrix);
                for (int j = 0; j < NumColumns; ++j)
                {
                    for (int i = 0; i < NumRows; ++i)
                    {
                        int index1D = j * NumRows + i;
                        this.data[index1D] = binaryOperation(this.data[index1D], matrix[i, j]);
                    }
                }
            }
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i, j]). 
        /// The resulting matrix overwrites the entries of this <see cref="Matrix"/> instance.
        /// </summary>
        /// <param name="matrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
        ///     <see cref="Matrix"/> instance.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix"/> has different 
        ///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
        public void DoEntrywiseIntoThis(Matrix matrix, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, matrix);
            for (int i = 0; i < data.Length; ++i) this.data[i] = binaryOperation(this.data[i], matrix.data[i]);
        }

        /// <summary>
        /// See <see cref="IMatrixView.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        IMatrixView IMatrixView.DoToAllEntries(Func<double, double> unaryOperation) => DoToAllEntries(unaryOperation);

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// result[i, j] = <paramref name="unaryOperation"/>(this[i,j]). 
        /// The resulting matrix is written to a new <see cref="Matrix"/> and then returned.
        /// </summary>
        /// <param name="unaryOperation">A method that takes 1 argument and returns 1 result.</param>
        public Matrix DoToAllEntries(Func<double, double> unaryOperation)
        {
            var result = new double[NumRows * NumColumns];
            for (int i = 0; i < NumRows * NumColumns; ++i) result[i] = unaryOperation(data[i]);
            return new Matrix(result, NumRows, NumColumns);
        }

        /// <summary>
        /// See <see cref="IMatrix.DoToAllEntriesIntoThis(Func{double, double})"/>.
        /// </summary>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            for (int i = 0; i < NumRows * NumColumns; ++i) data[i] = unaryOperation(data[i]);
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            if (other is Matrix casted)
            {
                //Check each dimension, rather than the lengths of the internal buffers
                if (!Preconditions.AreSameMatrixDimensions(this, casted)) return false;
                double[] otherData = casted.data;
                var comparer = new ValueComparer(tolerance);
                for (int i = 0; i < this.data.Length; ++i)
                {
                    if (!comparer.AreEqual(this.data[i], otherData[i])) return false;
                }
                return true;
            }
            else return other.Equals(this, tolerance); // To avoid accessing zero entries
        }

        /// <summary>
        /// Calculates the Cholesky factorization of a symmetric positive definite matrix with n = <see cref="NumRows"/> = 
        /// <see cref="NumColumns"/>, such that A = L^T * L. L is a lower triangular n-by-n matrix. This only works if the matrix
        /// is symmetric positive definite. Requires extra available memory n^2 entries. 
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        /// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not symmetric positive definite.</exception>
        /// <exception cref="MklException">Thrown if the call to Intel MKL fails due to invalid input.</exception>
        public CholeskyFull FactorCholesky()
        {
            Preconditions.CheckSquare(this);
            // Copy matrix. This may exceed available memory and needs an extra O(n^2) space. 
            // To avoid these, set "inPlace=true".
            return CholeskyFull.Factorize(NumColumns, CopyInternalData());
        }

        /// <summary>
        /// Calculates the LQ factorization of a matrix with m = <see cref="NumRows"/> &lt;= <see cref="NumColumns"/> = n, such 
        /// that A = L * Q. Q is an orthogonal n-by-n matrix and L is a lower trapezoidal m-by-n matrix. Requires extra available  
        /// memory form * n + min(m, n) entries.
        /// </summary>
        /// <exception cref="MklException">Thrown if the call to Intel MKL fails due to invalid input.</exception>
        public LQFactorization FactorLQ()
        {
            return LQFactorization.Factorize(NumRows, NumColumns, CopyInternalData());
        }

        /// <summary>
        /// Calculates the LUP factorization of a square matrix with n = <see cref="NumRows"/> = <see cref="NumColumns"/>, such 
        /// that A = P * L * U. L is a lower triangular n-by-n matrix. U is an upper triangular n-by-n matrix. P is an n-by-n
        /// permutation matrix. Requires extra available memory n^2 + n entries. 
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        /// <exception cref="MklException">Thrown if the call to Intel MKL fails due to invalid input.</exception>
        public LUFactorization FactorLU()
        {
            Preconditions.CheckSquare(this);
            // Copy matrix. This may exceed available memory and needs an extra O(n^2) space. 
            // To avoid these, set "inPlace=true".
            return LUFactorization.Factorize(NumColumns, CopyInternalData());
        }

        /// <summary>
        /// Calculates the QR factorization of a matrix with m = <see cref="NumRows"/> &gt;= <see cref="NumColumns"/> = n, such 
        /// that A = Q * R. Q is an orthogonal m-by-m matrix and R is an upper trapezoidal m-by-n matrix. Requires extra 
        /// available memory for m * n + min(m, n) entries. 
        /// </summary>
        /// <exception cref="MklException">Thrown if the call to Intel MKL fails due to invalid input.</exception>
        public QRFactorization FactorQR()
        {
            return QRFactorization.Factorize(NumRows, NumColumns, CopyInternalData());
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetColumn(int)"/>.
        /// </summary>
        public Vector GetColumn(int colIndex)
        {
            double[] result = new double[NumRows];
            Array.Copy(data, colIndex * NumRows, result, 0, NumRows);
            return Vector.CreateFromArray(result, false);
        }

        /// <summary>
        /// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal. The matrix must be square.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        public Vector GetDiagonal() => Vector.CreateFromArray(GetDiagonalAsArray(), false);

        /// <summary>
        /// Returns an array with the entries of the matrix's main diagonal. The matrix must be square.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        public double[] GetDiagonalAsArray()
        {
            Preconditions.CheckSquare(this);
            double[] diag = new double[NumRows];
            for (int i = 0; i < NumRows; ++i) diag[i] = data[i * NumRows + i];
            return diag;
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetRow(int)"/>.
        /// </summary>
        public Vector GetRow(int rowIndex)
        {
            double[] result = new double[NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                result[j] = data[j * NumRows + rowIndex];
            }
            return Vector.CreateFromArray(result, false);
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetSubmatrix(int[], int[])"/>.
        /// </summary>
        public Matrix GetSubmatrix(int[] rowIndices, int[] colIndices)
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
        /// See <see cref="ISliceable2D.GetSubmatrix(int, int, int, int)"/>.
        /// </summary>
        public Matrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
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

        /// <summary>
        /// Calculates the inverse matrix and returns it in a new <see cref="Matrix"/> instance. This only works if this 
        /// <see cref="Matrix"/> is square and invertible. If the determinant matrix is also needed, use 
        /// <see cref="InvertAndDetermninant"/> instead.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        /// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
        /// <exception cref="MklException">Thrown if the call to Intel MKL fails due to invalid input.</exception>
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

        /// <summary>
        /// Calculates the determinant and the inverse matrix and returns the latter in a new <see cref="Matrix"/> instance. 
        /// This only works if this <see cref="Matrix"/> is square and invertible.
        /// </summary>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
        /// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
        /// <exception cref="MklException">Thrown if the call to Intel MKL fails due to invalid input.</exception>
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

        /// <summary>
        /// Returns true if: this[i, j] &lt;= <paramref name="tolerance"/>, for 0 &lt;= i &lt; <see cref="NumRows"/>, 
        /// 0 &lt;= j &lt; <see cref="NumColumns"/>. Otherwise false is returned.
        /// </summary>
        /// <param name="tolerance">The tolerance under which a matrix entry is considered to be 0. It can be set to 0, to check 
        ///     if the entries are exactly 0.</param>
        public bool IsZero(double tolerance) => DenseStrategies.IsZero(data, tolerance);

        /// <summary>
        /// See <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>.
        /// </summary>
        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is Matrix casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else return otherMatrix.LinearCombination(otherCoefficient, this, thisCoefficient); // To avoid accessing zero entries
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// result[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
        /// The resulting matrix is written to a new <see cref="Matrix"/> and then returned.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix"/>.</param>
        /// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
        ///     <see cref="Matrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
        public Matrix LinearCombination(double thisCoefficient, Matrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref result[0], 1);
            return new Matrix(result, NumRows, NumColumns);
        }

        /// <summary>
        /// See <see cref="IMatrix.LinearCombinationIntoThis(double, IMatrixView, double)"/>.
        /// </summary>
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

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// this[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="Matrix"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix"/>.</param>
        /// <param name="otherMatrix">A matrix with the same <see cref="NumRows"/> and <see cref="NumColumns"/> as this 
        ///     <see cref="Matrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="NumRows"/> or <see cref="NumColumns"/> than this instance.</exception>
        public void LinearCombinationIntoThis(double thisCoefficient, Matrix otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref this.data[0], 1);
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return other.MultiplyRight(this, transposeOther, transposeThis);
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            if (other is Matrix) return MultiplyRight((Matrix)other, transposeThis);
            else return other.MultiplyLeft(this, transposeOther, transposeThis);
        }

        /// <summary>
        /// Performs the matrix-matrix multiplication: oper(this) * oper(<paramref name="other"/>).
        /// </summary>
        /// <param name="other">A matrix such that the <see cref="NumRows"/> of oper(<paramref name="other"/>) 
        ///     are equal to the <see cref="NumColumns"/> of oper(this).</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <param name="transposeOther">If true, oper(<paramref name="other"/>) = transpose(<paramref name="other"/>). 
        ///     Otherwise oper(<paramref name="other"/>) = <paramref name="other"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if oper(<paramref name="otherMatrix"/>) has 
        ///     different <see cref="NumRows"/> than the <see cref="NumColumns"/> of oper(this).</exception>
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

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IVectorView, bool)"/>.
        /// </summary>
        public Vector MultiplyRight(IVectorView vector, bool transposeThis = false)
        {
            if (vector is Vector casted) return MultiplyRight(casted, transposeThis);
            else throw new NotImplementedException();
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: oper(this) * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// </summary>
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to the 
        ///     <see cref="IIndexable2D.NumColumns"/> of oper(this).</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="vector"/> is different than the <see cref="NumColumns"/> of oper(this).</exception>
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

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            for (int i = 0; i < data.Length; ++i) aggregator = processEntry(data[i], aggregator);
            // no zeros implied
            return finalize(aggregator);
        }

        /// <summary>
        /// Creates a new <see cref="Matrix"/> that contains the entries of this <see cref="Matrix"/> with a different order,
        /// which is specified by the provided <paramref name="permutation"/> and <paramref name="oldToNew"/>.
        /// </summary>
        /// <param name="permutation">An array that contains the row/column indices of this <see cref="Matrix"/> in a 
        ///     different order.</param>
        /// <param name="oldToNew">If true, 
        ///     reordered[<paramref name="permutation"/>[i], <paramref name="permutation"/>[j]] =  original[i, j]. If false, 
        ///     reordered[i, j] = original[<paramref name="permutation"/>[i], <paramref name="permutation"/>[j]].</param>
        public Matrix Reorder(IReadOnlyList<int> permutation, bool oldToNew)
        {
            Preconditions.CheckSquare(this);
            if (permutation.Count != NumRows) throw new NonMatchingDimensionsException(
                $"This matrix has order = {NumRows}, while the permutation vector has {permutation.Count} entries.");
            if (oldToNew) return new Matrix(ArrayColMajor.ReorderOldToNew(NumRows, data, permutation), NumRows, NumRows);
            else return new Matrix(ArrayColMajor.ReorderNewToOld(NumRows, data, permutation), NumRows, NumRows);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Scale(double)"/>.
        /// </summary>
        IMatrixView IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="NumRows"/>, 0 &lt;= j &lt; <see cref="NumColumns"/>:
        /// result[i, j] = <paramref name="scalar"/> * this[i, j].
        /// The resulting matrix is written to a new <see cref="Matrix"/> and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
        public Matrix Scale(double scalar)
        {
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(data, result, data.Length);
            CBlas.Dscal(data.Length, scalar, ref result[0], 1);
            return new Matrix(result, NumRows, NumColumns);
        }

        /// <summary>
        /// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar) => CBlas.Dscal(data.Length, scalar, ref data[0], 1);

        /// <summary>
        /// Sets all entries of this matrix to be equal to <paramref name="value"/>.
        /// </summary>
        /// <param name="value">The value that all entries of the this matrix will be equal to.</param>
        public void SetAll(double value)
        {
            for (int i = 0; i < data.Length; ++i) data[i] = value;
        }

        /// <summary>
        /// Sets some consecutive entries of the column with index = <paramref name="colIdx"/> to be equal to 
        /// <paramref name="colValues"/>, starting from the entry with row index = <paramref name="rowStart"/>.
        /// </summary>
        /// <param name="colIdx">The index of the column to set. Constraints: 
        ///     0 &lt;= <paramref name="colIdx"/> &lt; <see cref="NumColumns"/>.</param>
        /// <param name="rowStart">The first entry of column <paramref name="colIdx"/> to be modified. Constraints: 
        ///     1) 0 &lt;= <paramref name="rowStart"/> &lt; <see cref="NumRows"/>, 
        ///     2) <paramref name="rowStart"/> + <paramref name="colValues"/>.<see cref="IIndexable1D.Length"/> &lt;= 
        ///        <see cref="NumRows"/>.</param>
        /// <param name="colValues">The new values of the column entries. Constraints: <paramref name="rowStart"/>
        ///     + <paramref name="colValues"/>.<see cref="IIndexable1D.Length"/> &lt;= <see cref="NumRows"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="colIdx"/> or <paramref name="rowStart"/> 
        ///     violate the described constraints.</exception>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rowStart"/>
        ///     + <paramref name="colValues"/>.<see cref="IIndexable1D.Length"/> &gt; <see cref="NumRows"/>.</exception>
        public void SetSubcolumn(int colIdx, Vector colValues, int rowStart = 0)
        {
            Preconditions.CheckIndexCol(this, colIdx);
            if (rowStart + colValues.Length > this.NumRows) throw new NonMatchingDimensionsException(
                "The entries to set exceed this matrix's number of rows");
            ArrayColMajor.SetCol(NumRows, NumColumns, data, colIdx, rowStart, colValues.InternalData);
        }

        /// <summary>
        /// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            data[colIdx * NumRows + rowIdx] = value;
        }

        /// <summary>
        /// Sets some consecutive entries of this matrix to be equal to the entries of <paramref name="submatrix"/>, starting 
        /// from the entry at (<paramref name="rowStart"/>, <paramref name="colStart"/>).
        /// </summary>
        /// <param name="rowStart">The index of the first row to be modified. Constraints: 
        ///     1) 0 &lt;= <paramref name="rowStart"/> &lt; <see cref="NumRows"/>, 
        ///     2) <paramref name="rowStart"/> + <paramref name="submatrix"/>.<see cref="NumRows"/> &lt;= 
        ///        this.<see cref="NumRows"/>.</param>
        /// <param name="colStart">The index of the first column to be modified. Constraints: 
        ///     1) 0 &lt;= <paramref name="colStart"/> &lt; <see cref="NumColumns"/>, 
        ///     2) <paramref name="colStart"/> + <paramref name="submatrix"/>.<see cref="NumColumns"/> &lt;= 
        ///        this.<see cref="NumColumns"/>.</param>
        /// <param name="submatrix">The new values of this matrix's entries to be modified. Constraints:
        ///     1) <paramref name="rowStart"/> + <paramref name="submatrix"/>.<see cref="NumRows"/> &lt;= 
        ///        this.<see cref="NumRows"/>,
        ///     2) <paramref name="colStart"/> + <paramref name="submatrix"/>.<see cref="NumColumns"/> &lt;= 
        ///        this.<see cref="NumColumns"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowStart"/> or <paramref name="colStart"/> 
        ///     violate the described constraints.</exception>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="rowStart"/>, <paramref name="colStart"/>
        ///     or <paramref name="submatrix"/> exceed the desrcibed constraints.</exception>
        public void SetSubmatrix(int rowStart, int colStart, Matrix submatrix)
        {
            // TODO: create Preconditions.CheckOverflow1D() and 2D for such setters.
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
        /// Sets some consecutive entries of the row with index = <paramref name="rowIdx"/> to be equal to 
        /// <paramref name="rowValues"/>, starting from the entry with column index = <paramref name="colStart"/>.
        /// </summary>
        /// <param name="rowIdx">The index of the row to set. Constraints: 
        ///     0 &lt;= <paramref name="rowIdx"/> &lt; <see cref="NumRows"/>.</param>
        /// <param name="colStart">The first entry of row <paramref name="rowIdx"/> to be modified. Constraints: 
        ///     1) 0 &lt;= <paramref name="colStart"/> &lt; <see cref="NumColumns"/>, 
        ///     2) <paramref name="colStart"/> + <paramref name="rowValues"/>.<see cref="IIndexable1D.Length"/> &lt;= 
        ///        <see cref="NumColumns"/>.</param>
        /// <param name="rowValues">The new values of the row entries. Constraints: <paramref name="colStart"/>
        ///     + <paramref name="rowValues"/>.<see cref="IIndexable1D.Length"/> &lt;= <see cref="NumColumns"/>.</param>
        /// <exception cref="IndexOutOfRangeException">Thrown if <paramref name="rowIdx"/> or <paramref name="colStart"/> 
        ///     violate the described constraints.</exception>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="colStart"/>
        ///     + <paramref name="rowValues"/>.<see cref="IIndexable1D.Length"/> &gt; <see cref="NumColumns"/>.</exception>
        public void SetSubrow(int rowIdx, Vector rowValues, int colStart = 0)
        {
            Preconditions.CheckIndexRow(this, rowIdx);
            if (colStart + rowValues.Length > this.NumRows) throw new NonMatchingDimensionsException(
                "The entries to set exceed this matrix's number of columns");
            ArrayColMajor.SetRow(NumRows, NumColumns, data, rowIdx, colStart, rowValues.InternalData);
        }

        /// <summary>
        /// Calculates the Singular Value Decomposition of a matrix.
        /// </summary>
        /// <param name="w"></param>
        /// <param name="v"></param>
        public void SVD(double[] w, double[,] v)
        {
            DenseStrategies.SVD(this, w, v);
        }

        /// <summary>
        /// Creates a new instance of the legacy matrix class <see cref="Numerical.LinearAlgebra.Matrix"/>, by copying the 
        /// entries of this <see cref="Matrix"/> instance. 
        /// </summary>
        public Numerical.LinearAlgebra.Interfaces.IMatrix2D ToLegacyMatrix()
        {
            return new Numerical.LinearAlgebra.Matrix2D(CopyToArray2D());
        }

        /// <summary>
        /// See <see cref="IMatrixView.Transpose"/>.
        /// </summary>
        IMatrixView IMatrixView.Transpose() => Transpose();

        /// <summary>
        /// Initializes a new <see cref="Matrix"/> instance, that is transpose to this: result[i, j] = this[j, i]. The entries 
        /// will be explicitly copied.
        /// </summary>
        public Matrix Transpose()
        {
            //TODO: The wrapper library does not include MKL's blas-like extensions yet. Create my own wrapper or 
            // piggyback on another BLAS function.
            double[] transpose = Conversions.ColumnMajorToRowMajor(data, NumRows, NumColumns);
            return new Matrix(transpose, NumColumns, NumRows);
        }

        /// <summary>
        /// Transposes the matrix by modifying the entries of this <see cref="Matrix instance"/>: this[i, j] = this[j, i].
        /// </summary>
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
