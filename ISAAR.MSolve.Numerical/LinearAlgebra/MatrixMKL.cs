using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;

//TODO: align data using mkl_malloc
namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    /// <summary>
    /// General matrix. Dense (full) storage. Uses MKL. Stored as 1D column major array.
    /// </summary>
    public class MatrixMKL
    {
        protected readonly double[] data;

        protected MatrixMKL(double[] data, int numRows, int numColumns)
        {
            this.data = data;
            this.NumRows = numRows;
            this.NumColumns = numColumns;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix
        /// </summary>
        public int NumRows { get; }

        /// <summary>
        /// The entry with row index = i and column index = j. 
        /// </summary>
        /// <param name="i">The row index: 0 &lt;= i &lt; <see cref="NumRows"/></param>
        /// <param name="j">The column index: 0 &lt;= j &lt; <see cref="NumColumns"/></param>
        /// <returns>The entry with indices i, j</returns>
        public double this[int i, int j]
        {
            get { return data[j * NumRows + i]; }
            set { data[j * NumRows + i] = value; }
        }

        /// <summary>
        /// Create a new <see cref="MatrixMKL"/> from a provided array. The array will be copied.
        /// </summary>
        /// <param name="array2D">A 2-dimensional array containing the elements of the matrix</param>
        /// <returns></returns>
        public static MatrixMKL CreateFromArray(double[,] array2D)
        {
            int numRows = array2D.GetLength(0);
            int numCols = array2D.GetLength(1); 
            return new MatrixMKL(Conversions.Array2DToFullColMajor(array2D), numRows, numCols);
        }

        /// <summary>
        /// Create a new <see cref="MatrixMKL"/> with the specified dimensions and all entries equal to 0.
        /// </summary> 
        /// <param name="numRows">The number of rows of the matrix.</param>
        /// <param name="numColumns">The number of rows of the matrix.</param>
        /// <returns></returns>
        public static MatrixMKL CreateZero(int numRows, int numColumns)
        {
            double[] data = new double[numRows * numColumns];
            return new MatrixMKL(data, numRows, numColumns);
        }

        /// <summary>
        /// Matrix vector multiplication, with the vector on the right: matrix * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
        /// <returns></returns>
        public DenseVector MultiplyRight(DenseVector vector)
        {
            Preconditions.CheckMultiplicationDimensions(this, vector);
            double[] result = new double[NumRows];
            CBlas.Dgemv(CBLAS_LAYOUT.CblasColMajor, CBLAS_TRANSPOSE.CblasNoTrans, NumRows, NumColumns,
                1.0, ref data[0], NumRows, ref vector.InternalData[0], 1, 0.0, ref result[0], 1);
            return DenseVector.CreateFromArray(result, false);
        }
    }
}
