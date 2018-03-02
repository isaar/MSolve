using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IntelMKL.LP64;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Reduction;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

//TODO: align data using mkl_malloc
namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    /// <summary>
    /// General matrix. Dense (full) storage. Uses MKL. Stored as 1D column major array.
    /// </summary>
    public class MatrixMKL: IMatrixViewMKL
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
        /// The number of rows of the matrix.
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
        /// Create a new <see cref="MatrixMKL"/> from a provided array. The array will can be copied for extra safety or not for 
        /// extra performance.
        /// </summary>
        /// <param name="array1D">A 1-dimensional array containing the elements of the matrix in column major order. Its length 
        /// must be equal to <see cref="numRows"/> + <see cref="NumColumns"/>. It will not be checked.</param>
        /// <param name="numRows">The number of rows of the matrix</param>
        /// <param name="numColumns">The number of columns of the matrix</param>
        /// <param name="copyArray">True to make a deep copy of <see cref="array1D"/>. 
        /// False (default) to use <see cref="array1D"/> as its internal storage.</param>
        /// <returns></returns>
        public static MatrixMKL CreateFromArray(double[] array1D, int numRows, int numColumns, bool copyArray = false)
        {
            if (copyArray)
            {
                var clone = new double[array1D.Length];
                Array.Copy(array1D, clone, clone.Length);
                return new MatrixMKL(clone, numRows, numColumns);
            }
            else
            {
                return new MatrixMKL(array1D, numRows, numColumns);
            }
        }

        /// <summary>
        /// The original matrix will be copied.
        /// </summary>
        /// <param name="original"></param>
        /// <returns></returns>
        public static MatrixMKL CreateFromMatrix(MatrixMKL original) 
        {
            //TODO: Perhaps this should use BLAS. 
            //TODO: Perhaps it should be an instance method CopyToMatrix(). Or the instance method would return an interface.
            double[] data = original.data;
            double[] clone = new double[data.Length];
            Array.Copy(data, clone, data.Length);
            return new MatrixMKL(clone, original.NumRows, original.NumColumns);
        }

        public static MatrixMKL CreateWithValue(int numRows, int numColumns, double value)
        {
            double[] data = new double[numRows * numColumns];
            for (int i = 0; i < data.Length; ++i) data[i] = value;
            return new MatrixMKL(data, numRows, numColumns);
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
        /// Copy the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
        /// and length(1) = <see cref="NumColumns"/>. 
        /// </summary>
        /// <returns>A new <see cref="double"/>[<see cref="NumRows"/>, <see cref="NumRows"/>] array 
        /// with the entries of the matrix</returns>
        public double[,] CopyToArray2D()
        {
            return Conversions.FullColMajorToArray2D(data, NumRows, NumColumns);
        }

        public bool Equals(MatrixMKL other, ValueComparer comparer = null)
        {
            if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
            if (comparer == null) comparer = new ValueComparer(1e-13);
            for (int i = 0; i < this.data.Length; ++i)
            {
                if (!comparer.AreEqual(this.data[i], other.data[i])) return false;
            }
            return true;
        }

        /// <summary>
        /// Matrix vector multiplication, with the vector on the right: matrix * vector.
        /// </summary>
        /// <param name="vector">A vector with length equal to <see cref="NumColumns"/>.</param>
        /// <returns></returns>
        public VectorMKL MultiplyRight(VectorMKL vector)
        {
            Preconditions.CheckMultiplicationDimensions(this, vector);
            double[] result = new double[NumRows];
            CBlas.Dgemv(CBLAS_LAYOUT.CblasColMajor, CBLAS_TRANSPOSE.CblasNoTrans, NumRows, NumColumns,
                1.0, ref data[0], NumRows, ref vector.InternalData[0], 1, 0.0, ref result[0], 1);
            return VectorMKL.CreateFromArray(result, false);
        }

        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double accumulator = identityValue;
            for (int i = 0; i < data.Length; ++i) accumulator = processEntry(data[i], accumulator);
            // no zeros implied
            return finalize(accumulator);
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
        public MatrixMKL Slice(int[] rowIndices, int[] colIndices)
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
            return new MatrixMKL(submatrix, rowIndices.Length, colIndices.Length);
        }

        /// <summary>
        /// Returns a subvector containing the entries at the indices between the provided start (inclusive) and end (exclusive).
        /// </summary>
        /// <param name="rowStartInclusive">The first row from which to copy entries.</param>
        /// <param name="rowEndExclusive">The row after the last one until which to copy entries.</param>
        /// <param name="colStartInclusive">The first column from which to copy entries.</param>
        /// <param name="colEndExclusive">The column after the last one until which to copy entries.</param>
        /// <returns></returns>
        public MatrixMKL Slice(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
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
            return new MatrixMKL(submatrix, newNumRows, newNumCols);
        }

        public void WriteToConsole(Array2DFormatting format = null)
        {
            if (format == null) format = Array2DFormatting.Brackets;
            Console.Write(format.ArrayStart);

            // First row
            Console.Write(format.RowSeparator + format.RowStart);
            Console.Write(data[0]);
            for (int j = 1; j < NumColumns; ++j)
            {
                Console.Write(format.ColSeparator + data[j * NumRows]);
            }
            Console.Write(format.RowEnd);

            // Subsequent rows
            for (int i = 0; i < NumRows; ++i)
            {
                Console.Write(format.RowSeparator + format.RowStart);
                Console.Write(data[i]);
                for (int j = 1; j < NumColumns; ++j)
                {
                    Console.Write(format.ColSeparator + data[j * NumRows + i]);
                }
                Console.Write(format.RowEnd);
            }
            Console.Write(format.RowSeparator + format.ArrayEnd);
        }

        /// <summary>
        /// Write the entries of the matrix to a specified file. If the file doesn't exist a new one will be created.
        /// </summary>
        /// <param name="path">The path of the file and its extension.</param>
        /// <param name="append">If the file already exists: Pass <see cref="append"/> = true to write after the current end of 
        ///     the file. Pass<see cref="append"/> = false to overwrite the file.</param>
        /// <param name="format">Formatting options for how to print the vector entries.</param>
        public void WriteToFile(string path, bool append = false, Array2DFormatting format = null)
        {
            if (format == null) format = Array2DFormatting.Plain;
            //TODO: incorporate this and WriteToConsole into a common function, where the user passes the stream and an object to 
            //deal with formating. Also add support for relative paths. Actually these methods belong in the "Logging" project, 
            // but since they are extremely useful they are implemented here for now.
            using (var writer = new StreamWriter(path, append))
            {
                writer.Write(format.ArrayStart);

                // First row
                writer.Write(format.RowSeparator + format.RowStart);
                writer.Write(data[0]);
                for (int j = 1; j < NumColumns; ++j)
                {
                    writer.Write(format.ColSeparator + data[j * NumRows]);
                }
                writer.Write(format.RowEnd);

                // Subsequent rows
                for (int i = 0; i < NumRows; ++i)
                {
                    writer.Write(format.RowSeparator + format.RowStart);
                    writer.Write(data[i]);
                    for (int j = 1; j < NumColumns; ++j)
                    {
                        writer.Write(format.ColSeparator + data[j * NumRows + i]);
                    }
                    writer.Write(format.RowEnd);
                }
                writer.Write(format.RowSeparator + format.ArrayEnd);
            }
        }
    }
}
