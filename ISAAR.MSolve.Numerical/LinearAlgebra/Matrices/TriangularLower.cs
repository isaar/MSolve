using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    ///  Upper triangular matrix. Packed storage (only stores the n*(n+1)/2 non zeros) in row major order. Uses MKL. For the
    ///  layout see 
    ///  <see cref="https://software.intel.com/en-us/mkl-developer-reference-c-matrix-storage-schemes-for-lapack-routines."/>
    /// </summary>
    public class TriangularLower: IEntrywiseOperable, IIndexable2D, IWriteable
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
        /// lower triangle, that is if colIdx &gt; rowIdx, a <see cref="SparsityPatternModifiedException"/> will be thrown. This 
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
        /// Create a new <see cref="TriangularLower"/> from the lower (subdiagonal) portion of the provided array. The array
        /// entries will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional containing the elements of the matrix. 
        ///     Its lengths in both dimensions must be the same.</param>
        /// <returns></returns>
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
        /// Create a new <see cref="TriangularLower"/> from a provided array. The array can be copied (for extra safety)
        /// or not (for extra performance).
        /// </summary>
        /// <param name="array1D">1-dimensional array containing the elements of the lower triangle of the matrix in row 
        ///     major order.</param>
        /// <param name="copyArray">True to make a deep copy of <see cref="array1D"/>. 
        ///     False (default) to use <see cref="array1D"/> as its internal storage.</param>
        /// <returns></returns>
        public static TriangularLower CreateFromArray(double[] array1D, bool copyArray = false)
        {
            int order = Conversions.PackedLengthToOrder(array1D.Length);
            if (copyArray)
            {
                var clone = new double[array1D.Length];
                Array.Copy(array1D, clone, array1D.Length);
                return new TriangularLower(clone, order);
            }
            else return new TriangularLower(array1D, order);
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
            Matrix fullMatrix = Matrix.CreateZero(Order, Order);
            for (int i = 0; i < Order; ++i) //Row major order
            {
                for (int j = 0; j <= i; ++j) fullMatrix[i, j] = data[FindIndex1D(i, j)];
            }
            return fullMatrix;
        }

        public double CalcDeterminant()
        {
            // TODO: Find more effienct formulas for the diagonal accesses.
            double det = 1.0;
            for (int i = 0; i < Order; ++i) det *= data[FindIndex1D(i, i)];
            return det;
        }

        public IEntrywiseOperable DoEntrywise(IEntrywiseOperable other, Func<double, double, double> binaryOperation)
        {
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        public IEntrywiseOperable DoToAllEntries(Func<double, double> unaryOperation)
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

        public bool Equals(IIndexable2D other, double tolerance = 1e-13)
        {
            return DenseStrategies.AreEqual(this, other, tolerance);
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
            CBLAS_TRANSPOSE transpose = transposeThis ? CBLAS_TRANSPOSE.CblasTrans : CBLAS_TRANSPOSE.CblasNoTrans;
            Preconditions.CheckMultiplicationDimensions(Order, vector.Length);
            double[] result = vector.CopyToArray();
            CBlas.Dtpmv(CBLAS_LAYOUT.CblasRowMajor, CBLAS_UPLO.CblasLower, transpose, CBLAS_DIAG.CblasNonUnit, Order,
                ref data[0], ref result[0], 1);
            return VectorMKL.CreateFromArray(result, false);
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
            CBlas.Dtpsv(CBLAS_LAYOUT.CblasRowMajor, CBLAS_UPLO.CblasLower, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit,
                Order, ref data[0], ref result[0], 1);
            return VectorMKL.CreateFromArray(result, false);
        }

        public void WriteToConsole()
        {
            DenseStrategies.WriteToConsole(this);
        }

        public void WriteToFile(string path, bool append = false)
        {
            DenseStrategies.WriteToFile(this, path, append);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int FindIndex1D(int i, int j)
        {
            return j + ((i + 1) * i) / 2;
        }
    }
}
