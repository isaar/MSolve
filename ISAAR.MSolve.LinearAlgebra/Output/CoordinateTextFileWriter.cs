using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes non zero entries of a matrix to Console or a file, using the Coordinate Text File format. The user may choose how 
    /// to format and justify them. For more information see https://math.nist.gov/MatrixMarket/formats.html#coord
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CoordinateTextFileWriter
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="CoordinateTextFileWriter"/> class.
        /// </summary>
        public CoordinateTextFileWriter()
        {
        }

        /// <summary>
        /// Describes how the matrix entries will be formatted and justified.
        /// </summary>
        public INumericFormat NumericFormat { get; set; } = new GeneralNumericFormat();

        /// <summary>
        /// Writes the non zero entries of the provided sparse matrix to Console.
        /// </summary>
        /// <param name="matrix">The sparse matrix to write.</param>
        public void WriteToConsole(ISparseMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteMatrix(matrix, writer));
        }

        /// <summary>
        /// Writes the non zero entries of the provided symmetric sparse matrix to Console.
        /// </summary>
        /// <param name="matrix">The symmetric sparse matrix to write.</param>
        public void WriteToConsole(ISparseSymmetricMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteSymmetricMatrix(matrix, writer));
        }

        /// <summary>
        /// Writes the non zero entries of the provided sparse matrix to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="matrix">The sparse matrix to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="matrix"/> will be written.</param>
        /// <param name="append">If true, <paramref name="matrix"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(ISparseMatrix matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteMatrix(matrix, writer), path, append);
        }

        /// <summary>
        /// Writes the non zero entries of the provided symmetric sparse matrix to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="matrix">The symmetric sparse matrix to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="matrix"/> will be written.</param>
        /// <param name="append">If true, <paramref name="matrix"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(ISparseSymmetricMatrix matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteSymmetricMatrix(matrix, writer), path, append);
        }

        private void WriteMatrix(ISparseMatrix matrix, StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            writer.Write($"{matrix.NumRows} {matrix.NumColumns} {matrix.CountNonZeros()}");
            foreach (var (row, col, val) in matrix.EnumerateNonZeros())
            {
                writer.WriteLine();
                writer.Write($"{row} {col} ");
                writer.Write(numberFormat, val);
            }
        }

        private void WriteSymmetricMatrix(ISparseSymmetricMatrix matrix, StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            writer.Write($"{matrix.NumRows} {matrix.NumColumns} {matrix.CountNonZerosUpper()}");
            foreach (var (row, col, val) in matrix.EnumerateNonZerosUpper())
            {
                writer.WriteLine();
                writer.Write($"{row} {col} ");
                writer.Write(numberFormat, val);
            }
        }
    }
}
