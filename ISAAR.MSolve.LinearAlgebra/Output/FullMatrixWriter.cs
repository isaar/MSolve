using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes all entries of a matrix to Console or a file. The user may choose how to format and justify them.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class FullMatrixWriter
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="FullMatrixWriter"/> class.
        /// </summary>
        public FullMatrixWriter()
        {
        }

        /// <summary>
        /// Describes how the matrix rows and columns will be separated.
        /// </summary>
        public Array2DFormat ArrayFormat { get; set; } = Array2DFormat.Plain;

        /// <summary>
        /// Describes how the matrix entries will be formatted and justified.
        /// </summary>
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        /// <summary>
        /// Writes the provided matrix to Console.
        /// </summary>
        /// <param name="matrix">The matrix to write.</param>
        public void WriteToConsole(IIndexable2D matrix)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(matrix, writer));
        }

        /// <summary>
        /// Writes the provided matrix to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="matrix">The matrix to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="matrix"/> will be written.</param>
        /// <param name="append">If true, <paramref name="matrix"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(IIndexable2D matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(matrix, writer), path, append);
        }

        private void WriteToStream(IIndexable2D matrix, StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            writer.Write(ArrayFormat.ArrayStart);

            // First row
            writer.Write(ArrayFormat.RowSeparator + ArrayFormat.RowStart);
            writer.Write(string.Format(numberFormat, matrix[0, 0]));
            for (int j = 1; j < matrix.NumColumns; ++j)
            {
                writer.Write(ArrayFormat.ColSeparator + string.Format(numberFormat, matrix[0, j]));
            }

            // Subsequent rows
            for (int i = 1; i < matrix.NumRows; ++i)
            {
                writer.Write(ArrayFormat.RowSeparator + ArrayFormat.RowStart);
                writer.Write(string.Format(numberFormat, matrix[i, 0]));
                for (int j = 1; j < matrix.NumColumns; ++j)
                {
                    writer.Write(ArrayFormat.ColSeparator + string.Format(numberFormat, matrix[i, j]));
                }
            }
            writer.Write(ArrayFormat.RowSeparator + ArrayFormat.ArrayEnd);
        }
    }
}
