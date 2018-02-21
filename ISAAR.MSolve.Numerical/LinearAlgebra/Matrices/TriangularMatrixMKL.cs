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

//TODO: Perhaps I should use row major for lower triangular, upper triangular or both.
//TODO: Perhaps I should have an abstract class that handles everything except the lower/upper specific stuff and concrete
//  private classes Lower, Upper. The indexer would be faster.
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    ///  Lower or upper triangular matrix. Packed storage (only stores the n*(n+1)/2 non zeros) in column major order. Uses MKL.
    /// </summary>
    class TriangularMatrixMKL: IMatrixViewMKL
    {
        /// <summary>
        /// Packed storage, column major order. 
        /// If lower triangular: A[i,j] = data[i + j*(2*n-j-1)/2] for 0 &lt;= j &lt;= i &lt; n.
        /// If upper triangular: A[i,j] = data[i + j*(j+1)/2] for 0 &lt;= i &lt;= j &lt; n.
        /// </summary>
        private readonly double[] data;

        private TriangularMatrixMKL(double[] data, int order, TrianglePosition triangle)
        {
            this.data = data;
            this.Order = order;
            this.NumRows = order;
            this.NumColumns = order;
            this.Triangle = triangle;
        }

        public enum TrianglePosition
        {
            Lower, Upper
        }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// The number of columns of the matrix.
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows or columns of the matrix.
        /// </summary>
        public int Order { get; }

        /// <summary>
        /// see cref="TriangularMatrixMKL.TrianglePosition.Lower"/> or <see cref="TriangularMatrixMKL.TrianglePosition.Upper"/>.
        /// </summary>
        public TrianglePosition Triangle { get; }

        /// <summary>
        /// The entry with row index = i and column index = j. Warning: trying to set an entry outside the non-zero triangle 
        /// will throw a <see cref="IndexOutOfRangeException"/>. This property is not that efficient, due to the necessary bound 
        /// checking.
        /// </summary>
        /// <param name="i">The row index: 0 &lt;= i &lt; <see cref="NumRows"/></param>
        /// <param name="j">The column index: 0 &lt;= j &lt; <see cref="NumColumns"/></param>
        /// <returns>The entry with indices i, j</returns>
        public double this[int i, int j]
        {
            get
            {
                if ((i < 0) || (i >= Order) || (j < 0) || (j >= Order))
                {
                    throw new IndexOutOfRangeException($"Invalid indices: ({i}, {j})");
                }
                if (Triangle == TrianglePosition.Lower)
                {
                    if (i >= j) return data[Index1DLowerColMajor(i, j)];
                    else return 0.0; 
                }
                else
                {
                    if (i <= j) return data[Index1DUpperColMajor(i, j)];
                    else return 0.0;
                }
            }
            set
            {
                if ((i < 0) || (i >= Order) || (j < 0) || (j >= Order))
                {
                    throw new IndexOutOfRangeException($"Invalid indices: ({i}, {j})");
                }
                if (Triangle == TrianglePosition.Lower)
                {
                    if (i >= j) data[Index1DLowerColMajor(i, j)] = value;
                    else throw new IndexOutOfRangeException($"Cannot change the superdiagonal entry A[{i}, {j}] = 0.");
                }
                else
                {
                    if (i <= j) data[Index1DUpperColMajor(i, j)] = value;
                    else throw new IndexOutOfRangeException($"Cannot change the subdiagonal entry A[{i}, {j}] = 0.");
                }
            }
        }

        /// <summary>
        /// Create a new <see cref="TriangularMatrixMKL"/> from the lower (subdiagonal) or upper (superdiagonal) portion of the 
        /// provided array. The array entries will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional containing the elements of the matrix. 
        ///     Its lengths in both dimensions must be the same.</param>
        /// <param name="triangle"><see cref="TriangularMatrixMKL.TrianglePosition.Lower"/> or 
        ///     <see cref="TriangularMatrixMKL.TrianglePosition.Upper"/>.</param>
        /// <returns></returns>
        public static TriangularMatrixMKL CreateFromArray(double[,] array2D, TrianglePosition trianglePosition)
        {
            int numRows = array2D.GetLength(0);
            int numCols = array2D.GetLength(1);
            if (numRows != numCols)
            {
                string msg = string.Format("Provided array must have the same dimensions, but was ({0}x{1})", numRows, numCols);
                throw new NonMatchingDimensionsException(msg);
            }

            if (trianglePosition == TrianglePosition.Lower)
            {
                return new TriangularMatrixMKL(Conversions.Array2DToPackedLowerColMajor(array2D), numRows, trianglePosition);
            }
            else
            {
                return new TriangularMatrixMKL(Conversions.Array2DToPackedUpperColMajor(array2D), numRows, trianglePosition);
            }
        }

        /// <summary>
        /// Create a new <see cref="TriangularMatrixMKL"/> from a provided array. The array can be copied (for extra safety)
        /// or not (for extra performance).
        /// </summary>
        /// <param name="array1D">A 1-dimensional array containing the elements of the matrix in column major order and matching 
        ///     <see cref="trianglePosition"/>.</param>
        /// <param name="trianglePosition"><see cref="TriangularMatrixMKL.TrianglePosition.Lower"/> or 
        ///     <see cref="TriangularMatrixMKL.TrianglePosition.Upper"/>.</param>
        /// <param name="copyArray">True (default) to make a deep copy of <see cref="array1D"/>. 
        ///     False to use <see cref="array1D"/> as its internal storage.</param>
        /// <returns></returns>
        public static TriangularMatrixMKL CreateFromArray(double[] array1D, TrianglePosition trianglePosition, 
            bool copyArray = true)
        {
            int order = Conversions.PackedLengthToOrder(array1D.Length);
            if (copyArray)
            {
                var clone = new double[array1D.Length];
                Array.Copy(array1D, clone, array1D.Length);
                return new TriangularMatrixMKL(clone, order, trianglePosition);
            }
            else
            {
                return new TriangularMatrixMKL(array1D, order, trianglePosition);
            }
        }

        /// <summary>
        /// Copy the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
        /// and length(1) = <see cref="Order"/>. 
        /// </summary>
        /// <returns>A new <see cref="double"/>[<see cref="Order"/>, <see cref="Order"/>] array 
        ///     with the entries of the matrix</returns>
        public double[,] CopyToArray2D()
        {
            if (Triangle == TrianglePosition.Lower)
            {
                return Conversions.PackedLowerColMajorToArray2D(data);
            }
            else
            {
                return Conversions.PackedUpperColMajorToArray2D(data);
            }
        }

        /// <summary>
        /// Matrix vector multiplication, with the vector on the right: matrix * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to <see cref="Order"/>.</param>
        /// <returns></returns>
        public DenseVector MultiplyRight(DenseVector vector)
        {
            Preconditions.CheckMultiplicationDimensions(this, vector);
            double[] result = vector.CopyToArray();
            CBLAS_UPLO uplo = (Triangle == TrianglePosition.Lower) ? CBLAS_UPLO.CblasLower : CBLAS_UPLO.CblasUpper;
            CBlas.Dtpmv(CBLAS_LAYOUT.CblasColMajor, uplo, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit, Order,
                ref data[0], ref result[0], 1);
            return DenseVector.CreateFromArray(result, false);
        }

        public double CalcDeterminant()
        {
            // TODO: Find more effienct formulas for the diagonal accesses.
            double det = 1.0;
            int n = Order;
            if (Triangle == TrianglePosition.Lower)
            {
                for (int i = 0; i < n; ++i)
                {
                    det *= data[Index1DLowerColMajor(i, i)];
                }
            }
            else
            {
                for (int i = 0; i < n; ++i)
                {
                    det *= data[Index1DUpperColMajor(i, i)];
                }
            }
            return det;
        }

        /// <summary>
        /// WARNING: No exception will be thrown if the matrix is singular.
        /// </summary>
        /// <param name="rhs"></param>
        /// <returns></returns>
        public DenseVector SolveLinearSystem(DenseVector rhs)
        {
            Preconditions.CheckSystemSolutionDimensions(this, rhs);
            double[] result = rhs.CopyToArray();
            CBLAS_UPLO uplo = (Triangle == TrianglePosition.Lower) ? CBLAS_UPLO.CblasLower : CBLAS_UPLO.CblasUpper;
            CBlas.Dtpsv(CBLAS_LAYOUT.CblasColMajor, uplo, CBLAS_TRANSPOSE.CblasNoTrans, CBLAS_DIAG.CblasNonUnit, Order,
                ref data[0], ref result[0], 1);
            return DenseVector.CreateFromArray(result, false);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int Index1DLowerColMajor(int i, int j)
        {
            return i + (j * (2 * Order - j - 1)) / 2;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int Index1DUpperColMajor(int i, int j)
        {
            return i + (j * (j + 1)) / 2;
        }
    }
}
