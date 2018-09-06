using System;
using System.Runtime.CompilerServices;
using IntelMKL.LP64;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: DoEntryWise may change the sparsity pattern if binaryOperation(0, 0)!=0. Fix this in all triangular and sparse matrices.
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Lower triangular square matrix in row major Packed storage format (only stores the n*(n+1)/2 non zeros). Uses Intel MKL. 
    /// For the more information about the layout, see 
    /// <see cref="https://software.intel.com/en-us/mkl-developer-reference-c-matrix-storage-schemes-for-lapack-routines."/>
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class TriangularLower: IMatrix
    {
        /// <summary>
        /// Packed storage, row major order: L[i, j] = data[j + (i+1)*i/2] for 0 &lt;= j &lt;= i &lt; n.
        /// </summary>
        private readonly double[] data;

        private TriangularLower(double[] data, int order)
        {
            this.data = data;
            this.Order = order;
        }

        /// <summary>
        /// The number of columns of the square matrix.
        /// </summary>
        public int NumColumns { get { return Order; } }

        /// <summary>
        /// The number of rows of the square matrix.
        /// </summary>
        public int NumRows { get { return Order; } }

        /// <summary>
        /// The number of rows/columns of the square matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>.
        /// </summary>
        /// <remarks> This property is not that efficient, due to the necessary bound checking.</remarks>
        /// <exception cref="SparsityPatternModifiedException">Thrown if you try to set an entry outside the non-zero
        ///     lower triangle, that is if <paramref name="colIdx"/> &gt; <paramref name="rowIdx"/>.</exception>
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                Preconditions.CheckIndices(this, rowIdx, colIdx);
                if (rowIdx >= colIdx) return data[FindIndex1D(rowIdx, colIdx)];
                else return 0.0;
            }
            set
            {
                Preconditions.CheckIndices(this, rowIdx, colIdx);
                if (rowIdx >= colIdx) data[FindIndex1D(rowIdx, colIdx)] = value;
                else throw new SparsityPatternModifiedException(
                    $"Cannot change the superdiagonal entry A[{rowIdx}, {colIdx}] = 0.");
            }
        }

        /// <summary>
        /// Initializes a new instance of <see cref="TriangularLower"/> by copying the lower triangle of the provided 2D array.
        /// </summary>
        /// <param name="array2D">A 2-dimensional array containing the elements of the matrix. Constraints: 
        ///     <paramref name="array2D"/>.GetLength(0) == <paramref name="array2D"/>.GetLength(1).</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="array2D"/>.GetLength(0) != 
        ///     <paramref name="array2D"/>.GetLength(1).</exception>
        public static TriangularLower CreateFromArray(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numCols = array2D.GetLength(1);
            if (numRows != numCols)
            {
                string msg = $"Provided array must have the same dimensions, but was ({numRows}x{numCols})";
                throw new NonMatchingDimensionsException(msg);
            }
            return new TriangularLower(Conversions.Array2DToPackedLowerRowMajor(array2D), numRows);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="TriangularLower"/> with <paramref name="array1D"/> or a clone as its 
        /// internal array.
        /// </summary>
        /// <param name="order">The number of rows/columns of the new square matrix.</param>
        /// <param name="array1D">An 1-dimensional array containing the elements of the lower triangle of the matrix in row 
        ///     major order.</param>
        /// <param name="copyArray">If true, <paramref name="array1D"/> will be copied and the new <see cref="TriangularLower"/>  
        ///     instance will have a reference to the copy, which is safer. If false, the new matrix will have a reference to 
        ///     <paramref name="array1D"/> itself, which is faster.</param>
        public static TriangularLower CreateFromArray(int order, double[] array1D, bool copyArray = false)
        {
            if (copyArray)
            {
                var clone = new double[array1D.Length];
                Array.Copy(array1D, clone, array1D.Length);
                return new TriangularLower(clone, order);
            }
            else return new TriangularLower(array1D, order);
        }

        /// <summary>
        /// Initializes a new instance of <see cref="TriangularLower"/> with all entries being equal to 0. Only the lower 
        /// triangle zero entries are explictily stored though.
        /// </summary>
        /// <param name="order">The number of rows/columns of the new square matrix.</param>
        public static TriangularLower CreateZero(int order)
        {
            var data = new double[(order * (order + 1)) / 2];
            return new TriangularLower(data, order);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
        /// </summary>
        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is TriangularLower casted) return Axpy(casted, otherCoefficient);
            else return DoEntrywise(otherMatrix, (x1, x2) => x1 + otherCoefficient * x2); //TODO: optimize this
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="Order"/>, 0 &lt;= j &lt;= i:
        /// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix is written to a new <see cref="TriangularLower"/> and then returned.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same <see cref="Order"/> as this <see cref="TriangularLower"/> 
        ///     instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="Order"/> than this instance.</exception>
        public TriangularLower Axpy(TriangularLower otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref result[0], 1);
            return new TriangularLower(result, NumRows);
        }

        /// <summary>
        /// See <see cref="IMatrix.AxpyIntoThis(IMatrixView, double)"/>.
        /// </summary>
        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is TriangularLower casted) AxpyIntoThis(casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix is also lower triangular.");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="Order"/>, 0 &lt;= j &lt;= i:
        /// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="TriangularLower"/> instance.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same <see cref="Order"/> as this <see cref="TriangularLower"/> 
        ///     instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="Order"/> than this instance.</exception>
        public void AxpyIntoThis(TriangularLower otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref this.data[0], 1);
        }

        /// <summary>
        /// Calculates the determinant of this matrix.
        /// </summary>
        public double CalcDeterminant()
        {
            // TODO: Find more effienct formulas for the diagonal accesses.
            double det = 1.0;
            for (int i = 0; i < Order; ++i) det *= data[FindIndex1D(i, i)];
            return det;
        }

        /// <summary>
        /// Copies the entries of the matrix into a 2-dimensional array. The returned array has length(0) = length(1) = 
        /// <see cref="Order"/>. 
        /// </summary>
        public double[,] CopyToArray2D() => Conversions.PackedLowerRowMajorToArray2D(data);

        /// <summary>
        /// Initializes a new <see cref="Matrix"/> instance by copying the entries of this <see cref="TriangularLower"/> into
        /// the lower triangle of the new matrix.
        /// </summary>
        public Matrix CopyToFullMatrix()
        {
            // TODO: This won't work if the implementation of Matrix changes
            double[] fullArray = Conversions.PackedLowerRowMajorToFullColMajor(data, Order); 
            return Matrix.CreateFromArray(fullArray, Order, Order, false);
        }

        /// <summary>
        /// See <see cref="IMatrixView.DoEntrywise(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public IMatrixView DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation)
        {
            return DenseStrategies.DoEntrywise(this, matrix, binaryOperation);
        }

        /// <summary>
        /// See <see cref="IMatrix.DoEntrywiseIntoThis(IMatrixView, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IMatrixView matrix, Func<double, double, double> binaryOperation)
        {
            if (matrix is TriangularLower casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix is also lower triangular.");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="Order"/>, 0 &lt;= j &lt;= i:
        /// this[i, j] = <paramref name="binaryOperation"/>(this[i,j], <paramref name="matrix"/>[i, j]) 
        /// The resulting matrix overwrites the entries of this <see cref="TriangularLower"/> instance.
        /// </summary>
        /// <param name="matrix">A matrix with the same <see cref="Order"/> as this <see cref="TriangularLower"/> 
        ///     instance.</param>
        /// <param name="binaryOperation">A method that takes 2 arguments and returns 1 result.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="matrix"/> has different 
        ///     <see cref="Order"/> than this instance.</exception>
        public void DoEntrywiseIntoThis(TriangularLower matrix, Func<double, double, double> binaryOperation)
        {
            Preconditions.CheckSameMatrixDimensions(this, matrix);
            for (int i = 0; i < data.Length; ++i) this.data[i] = binaryOperation(this.data[i], matrix.data[i]);
        }

        /// <summary>
        /// See <see cref="IMatrixView.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        public IMatrixView DoToAllEntries(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Only apply the operation on non zero entries
                double[] newData = new double[data.Length];
                for (int i = 0; i < data.Length; ++i) newData[i] = unaryOperation(data[i]);
                {
                    return new TriangularLower(newData, Order);
                }
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                Matrix fullMatrix = Matrix.CreateZero(Order, Order);
                for (int i = 0; i < Order; ++i) //Row major order
                {
                    for (int j = 0; j <= i; ++j) fullMatrix[i, j] = unaryOperation(data[FindIndex1D(i, j)]);
                }
                return fullMatrix;
            }
        }
        
        /// <summary>
        /// See <see cref="IMatrix.DoToAllEntriesIntoThis(Func{double, double})"/>.
        /// </summary>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0))
            {
                for (int i = 0; i < data.Length; ++i) data[i] = unaryOperation(data[i]);
            }
            else throw new SparsityPatternModifiedException("This operation will change the sparsity pattern");
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1e-13) => DenseStrategies.AreEqual(this, other, tolerance);

        /// <summary>
        /// See <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>.
        /// </summary>
        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is TriangularLower casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else return DoEntrywise(otherMatrix, (x1, x2) => thisCoefficient * x1 + otherCoefficient * x2); //TODO: optimize this
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="Order"/>, 0 &lt;= j &lt;= i:
        /// result[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
        /// The resulting matrix is written to a new <see cref="TriangularLower"/> and then returned.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="TriangularLower"/>.</param>
        /// <param name="otherMatrix">A matrix with the same <see cref="Order"/> as this <see cref="TriangularLower"/> 
        ///     instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="Order"/> than this instance.</exception>
        public TriangularLower LinearCombination(double thisCoefficient, TriangularLower otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref result[0], 1);
            return new TriangularLower(result, NumRows);
        }

        /// <summary>
        /// See <see cref="IMatrix.LinearCombinationIntoThis(double, IMatrixView, double)"/>.
        /// </summary>
        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is TriangularLower casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix is also lower triangular.");
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="Order"/>, 0 &lt;= j &lt;= i:
        /// this[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="TriangularLower"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="TriangularLower"/>.</param>
        /// <param name="otherMatrix">A matrix with the same <see cref="Order"/> as this <see cref="TriangularLower"/> 
        ///     instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     <see cref="Order"/> than this instance.</exception>
        public void LinearCombinationIntoThis(double thisCoefficient, TriangularLower otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref this.data[0], 1);
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyLeft(IMatrixView matrix, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(matrix, this, transposeOther, transposeThis);
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IMatrixView, bool, bool)"/>.
        /// </summary>
        public Matrix MultiplyRight(IMatrixView matrix, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(this, matrix, transposeThis, transposeOther);
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
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to <see cref="Order"/> of this 
        ///     matrix.</param>
        /// <param name="transposeThis">If true, oper(this) = transpose(this). Otherwise oper(this) = this.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="vector"/> is different than <see cref="Order"/> of this matrix.</exception>
        public Vector MultiplyRight(Vector vector, bool transposeThis = false)
        {
            CBLAS_TRANSPOSE transpose = transposeThis ? CBLAS_TRANSPOSE.CblasTrans : CBLAS_TRANSPOSE.CblasNoTrans;
            Preconditions.CheckMultiplicationDimensions(Order, vector.Length);
            double[] result = vector.CopyToArray();
            CBlas.Dtpmv(CBLAS_LAYOUT.CblasRowMajor, CBLAS_UPLO.CblasLower, transpose, CBLAS_DIAG.CblasNonUnit, Order,
                ref data[0], ref result[0], 1);
            return Vector.CreateFromArray(result, false);
        }

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int nnz = data.Length;
            for (int i = 0; i < nnz; ++i) aggregator = processEntry(data[i], aggregator);
            aggregator = processZeros(Order * Order - nnz, aggregator);
            return finalize(aggregator);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Scale(double)"/>.
        /// </summary>
        IMatrixView IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the following operation for 0 &lt;= i &lt; <see cref="Order"/>, 0 &lt;= j &lt;= i:
        /// result[i, j] = <paramref name="scalar"/> * this[i, j].
        /// The resulting matrix is written to a new <see cref="TriangularLower"/> and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
        public TriangularLower Scale(double scalar)
        {
            int nnz = this.data.Length;
            double[] result = new double[nnz];
            Array.Copy(this.data, result, nnz);
            CBlas.Dscal(nnz, scalar, ref result[0], 1);
            return new TriangularLower(result, this.Order);
        }

        /// <summary>
        /// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar) => CBlas.Dscal(data.Length, scalar, ref data[0], 1);

        /// <summary>
        /// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            if (colIdx > rowIdx) throw new SparsityPatternModifiedException("Cannot modify upper triangle entries.");
            this[rowIdx, colIdx] = value;
        }

        /// <summary>
        /// Solves the linear equations system: this * result = <paramref name="rhs"/> by forward substitution. WARNING: this
        /// matrix must be invertible. No exception will be thrown if the matrix is singular.
        /// </summary>
        /// <param name="rhs">The right hand side vector of the linear system. Constraints: 
        ///     <paramref name="rhs"/>.<see cref="Vector.Length"/> == this.<see cref="Order"/>.</param>
        public Vector SolveLinearSystem(Vector rhs)
        {
            Preconditions.CheckSystemSolutionDimensions(this, rhs);
            double[] result = rhs.CopyToArray();
            CBlas.Dtpsv(CBLAS_LAYOUT.CblasRowMajor, CBLAS_UPLO.CblasLower, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit,
                Order, ref data[0], ref result[0], 1);
            return Vector.CreateFromArray(result, false);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Transpose"/>.
        /// </summary>
        public IMatrixView Transpose() => Transpose(true);

        /// <summary>
        /// Creates a new <see cref="TriangularUpper"/> matrix, that is transpose to this: result[i, j] = this[j, i]. The 
        /// internal array can be copied or shared with this <see cref="TriangularLower"/> matrix.
        /// </summary>
        /// <param name="copyInternalArray">If true, the internal array that stores the entries of this 
        ///     <see cref="TriangularLower"/> instance will be copied and the new <see cref="TriangularUpper"/> instance 
        ///     will have a reference to the copy, which is safer. If false, both the new matrix and this one will have  
        ///     a reference to the same internal array, which is faster.</param>
        public TriangularUpper Transpose(bool copyInternalArray)
            => TriangularUpper.CreateFromArray(Order, data, copyInternalArray); // trans(lower row major) = upper col major

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int FindIndex1D(int i, int j)
        {
            return j + ((i + 1) * i) / 2;
        }
    }
}
