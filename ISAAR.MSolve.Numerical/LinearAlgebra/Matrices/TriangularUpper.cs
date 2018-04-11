using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Factorizations;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

//TODO: Perhaps I should use row major for lower triangular, upper triangular or both.
//TODO: Perhaps I should have an abstract class that handles everything except the lower/upper specific stuff and concrete
//  private classes Lower, Upper. The indexer would be faster.
//TODO: align data using mkl_malloc
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    ///  Upper triangular matrix. Packed storage (only stores the n*(n+1)/2 non zeros) in column major order. Uses MKL. For the
    ///  layout see 
    ///  <see cref="https://software.intel.com/en-us/mkl-developer-reference-c-matrix-storage-schemes-for-lapack-routines."/>
    /// </summary>
    public class TriangularUpper: IMatrix
    {
        /// <summary>
        /// Packed storage, column major order: U[i, j] = data[i + j*(2*n-j-1)/2] for 0 &lt;= j &lt;= i &lt; n.
        /// Although not used here, lower triangular would be: L[i, j] = data[i + j*(j+1)/2] for 0 &lt;= i &lt;= j &lt; n.
        /// </summary>
        private readonly double[] data;

        private TriangularUpper(double[] data, int order)
        {
            this.data = data;
            this.Order = order;
        }

        /// <summary>
        /// The number of columns of the matrix.
        /// </summary>
        public int NumColumns { get { return Order; } }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get { return Order; } }

        /// <summary>
        /// The number of rows or columns of the square matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// The entry with row index = rowIdx and column index = colIdx. Warning: if you try to set an entry outside the non-zero
        /// upper triangle, that is if rowIdx &gt; colIdx, a <see cref="SparsityPatternModifiedException"/> will be thrown. This 
        /// property is not that efficient, due to the necessary bound checking.
        /// </summary>
        /// <param name="rowIdx">The row index: 0 &lt;= i &lt; <see cref="Order"/></param>
        /// <param name="colIdx">The column index: 0 &lt;= j &lt; <see cref="Order"/></param>
        /// <returns>The entry with indices i, j</returns>
        public double this[int rowIdx, int colIdx]
        {
            get
            {
                Preconditions.CheckIndices(this, rowIdx, colIdx);
                if (colIdx >= rowIdx) return data[FindIndex1D(rowIdx, colIdx)];
                else return 0.0;
            }
            set
            {
                Preconditions.CheckIndices(this, rowIdx, colIdx);
                if (colIdx >= rowIdx) data[FindIndex1D(rowIdx, colIdx)] = value;
                else throw new SparsityPatternModifiedException(
                    $"Cannot change the subdiagonal entry A[{rowIdx}, {colIdx}] = 0.");
            }
        }

        /// <summary>
        /// Create a new <see cref="TriangularUpper"/> from the upper (superdiagonal) portion of the provided array. The array
        /// entries will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional containing the elements of the matrix. 
        ///     Its lengths in both dimensions must be the same.</param>
        /// <returns></returns>
        public static TriangularUpper CreateFromArray(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numCols = array2D.GetLength(1);
            if (numRows != numCols)
            {
                string msg = $"Provided array must have the same dimensions, but was ({numRows}x{numCols})";
                throw new NonMatchingDimensionsException(msg);
            }
            return new TriangularUpper(Conversions.Array2DToPackedUpperColMajor(array2D), numRows);
        }
        
        /// <summary>
        /// Create a new <see cref="TriangularUpper"/> from a provided array. The array can be copied (for extra safety)
        /// or not (for extra performance).
        /// </summary>
        /// <param name="array1D">1-dimensional array containing the elements of the upper triangle of the matrix in column 
        ///     major order.</param>
        /// <param name="copyArray">True to make a deep copy of <see cref="array1D"/>. 
        ///     False (default) to use <see cref="array1D"/> as its internal storage.</param>
        /// <returns></returns>
        public static TriangularUpper CreateFromArray(double[] array1D, bool copyArray = false)
        {
            int order = Conversions.PackedLengthToOrder(array1D.Length);
            if (copyArray)
            {
                var clone = new double[array1D.Length];
                Array.Copy(array1D, clone, array1D.Length);
                return new TriangularUpper(clone, order);
            }
            else return new TriangularUpper(array1D, order);
        }

        public static TriangularUpper CreateZero(int order)
        {
            var data = new double[(order * (order + 1)) / 2];
            return new TriangularUpper(data, order);
        }

        public IMatrixView Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is TriangularUpper casted) return Axpy(casted, otherCoefficient);
            else return DoEntrywise(otherMatrix, (x1, x2) => x1 + otherCoefficient * x2); //TODO: optimize this
        }

        public TriangularUpper Axpy(TriangularUpper otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref result[0], 1);
            return new TriangularUpper(result, NumColumns);
        }

        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is TriangularUpper casted) AxpyIntoThis(casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix is also upper triangular.");
        }

        public void AxpyIntoThis(TriangularUpper otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpy(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, ref this.data[0], 1);
        }

        /// <summary>
        /// Copy the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
        /// and length(1) = <see cref="Order"/>. 
        /// </summary>
        /// <returns>A new <see cref="double"/>[<see cref="Order"/>, <see cref="Order"/>] array 
        ///     with the upper triangular entries of the matrix</returns>
        public double[,] CopyToArray2D()
        {
            return Conversions.PackedUpperColMajorToArray2D(data);
        }

        public Matrix CopyToFullMatrix()
        {
            // TODO: This won't work if the implementation of Matrix changes
            double[] fullArray = Conversions.PackedUpperColMajorToFullColMajor(data, Order);
            return Matrix.CreateFromArray(fullArray, Order, Order, false);
        }

        public double CalcDeterminant()
        {
            // TODO: Find more effienct formulas for the diagonal accesses.
            double det = 1.0;
            for (int i = 0; i < Order; ++i) det *= data[FindIndex1D(i, i)];
            return det;
        }

        public IMatrixView DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            return DenseStrategies.DoEntrywise(this, other, binaryOperation); //TODO: this can be optimized.
        }

        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is TriangularUpper casted) DoEntrywiseIntoThis(casted, binaryOperation);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix is also upper triangular.");
        }

        public void DoEntrywiseIntoThis(TriangularUpper other, Func<double, double, double> binaryOperation)
        {
            //TODO: Aren't there any operations that would change the sparsity pattern, even if the other matrix is upper triangular?
            Preconditions.CheckSameMatrixDimensions(this, other);
            for (int i = 0; i < data.Length; ++i) this.data[i] = binaryOperation(this.data[i], other.data[i]);
        }

        public IMatrixView DoToAllEntries(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Only apply the operation on non zero entries
                double[] newData = new double[data.Length];
                for (int i = 0; i < data.Length; ++i) newData[i] = unaryOperation(data[i]);
                return new TriangularUpper(newData, Order);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                Matrix fullMatrix = Matrix.CreateZero(Order, Order);
                for (int j = 0; j < Order; ++j) //Column major order
                {
                    for (int i = 0; i <= j; ++i) fullMatrix[i, j] = unaryOperation(data[FindIndex1D(i, j)]);
                }
                return fullMatrix;
            }
        }

        void IMatrix.DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            DoToAllEntriesIntoThis(unaryOperation);
        }

        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0))
            {
                for (int i = 0; i < data.Length; ++i) data[i] = unaryOperation(data[i]);
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

        public IMatrixView LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is TriangularUpper casted) return LinearCombination(thisCoefficient, casted, otherCoefficient);
            else return DoEntrywise(otherMatrix, (x1, x2) => thisCoefficient * x1 + otherCoefficient * x2); //TODO: optimize this
        }

        public TriangularUpper LinearCombination(double thisCoefficient, TriangularUpper otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] result = new double[data.Length];
            Array.Copy(this.data, result, data.Length);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref result[0], 1);
            return new TriangularUpper(result, NumColumns);
        }

        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is TriangularUpper casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix is also upper triangular.");
        }

        public void LinearCombinationIntoThis(double thisCoefficient, TriangularUpper otherMatrix, double otherCoefficient)
        {
            Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
            CBlas.Daxpby(data.Length, otherCoefficient, ref otherMatrix.data[0], 1, thisCoefficient, ref this.data[0], 1);
        }

        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(other, this, transposeOther, transposeThis);
        }

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
        /// Matrix-vector multiplication, with the vector on the right: matrix * vector or transpose(matrix) * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
        /// <param name="transposeThis">Set to true to transpose this (the left matrix). Unless the transpose matrix is used in 
        ///     more than one multiplications, setting this flag to true is usually preferable to creating the transpose.</param>
        /// <returns></returns>
        public VectorMKL MultiplyRight(VectorMKL vector, bool transposeThis = false)
        {
            CBLAS_TRANSPOSE transpose = transposeThis ? CBLAS_TRANSPOSE.CblasTrans: CBLAS_TRANSPOSE.CblasNoTrans;
            Preconditions.CheckMultiplicationDimensions(Order, vector.Length);
            double[] result = vector.CopyToArray();
            CBlas.Dtpmv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasUpper, transpose, CBLAS_DIAG.CblasNonUnit, Order,
                ref data[0], ref result[0], 1);
            return VectorMKL.CreateFromArray(result, false);
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int nnz = data.Length;
            for (int i = 0; i < nnz; ++i) aggregator = processEntry(data[i], aggregator);
            aggregator = processZeros(Order * Order - nnz, aggregator);
            return finalize(aggregator);
        }

        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
        {
            this[rowIdx, colIdx] = value;
        }

        /// <summary>
        /// WARNING: No exception will be thrown if the matrix is singular.
        /// </summary>
        /// <param name="rhs"></param>
        /// <returns></returns>
        public VectorMKL SolveLinearSystem(VectorMKL rhs)
        {
            Preconditions.CheckSystemSolutionDimensions(this, rhs);
            double[] result = rhs.CopyToArray();
            CBlas.Dtpsv(CBLAS_LAYOUT.CblasColMajor, CBLAS_UPLO.CblasUpper, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit, 
                Order, ref data[0], ref result[0], 1);
            return VectorMKL.CreateFromArray(result, false);
        }

        public IMatrixView Transpose()
        {
            return Transpose(true);
        }

        public TriangularLower Transpose(bool copyInternalArray)
        {
            return TriangularLower.CreateFromArray(data, copyInternalArray); // trans(upper col major) = lower row major
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int FindIndex1D(int i, int j)
        {
            return i + ((j + 1) * j) / 2;
        }
    }
}
