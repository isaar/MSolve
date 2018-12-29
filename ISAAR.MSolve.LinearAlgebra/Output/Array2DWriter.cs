using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes all entries of a 2D array to Console or a file. The user may choose how to format and justify them.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Array2DWriter
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Array2DWriter"/> class.
        /// </summary>
        public Array2DWriter()
        {
        }

        /// <summary>
        /// Describes how the array entries will be separated.
        /// </summary>
        public Array2DFormat ArrayFormat { get; set; } = Array2DFormat.Plain;

        /// <summary>
        /// Describes how the array entries will be formatted and justified.
        /// </summary>
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        /// <summary>
        /// Writes the provided array to Console.
        /// </summary>
        /// <param name="array">The array to write.</param>
        public void WriteToConsole(double[,] array)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(array, writer));
        }

        /// <summary>
        /// Writes the provided array to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="array">The array to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="array"/> will be written.</param>
        /// <param name="append">If true, <paramref name="array"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(double[,] array, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(array, writer), path, append);
        }

        private void WriteToStream(double[,] array, StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            writer.Write(ArrayFormat.ArrayStart);

            // First row
            writer.Write(ArrayFormat.RowSeparator + ArrayFormat.RowStart);
            writer.Write(string.Format(numberFormat, array[0, 0]));
            for (int j = 1; j < array.GetLength(1); ++j)
            {
                writer.Write(ArrayFormat.ColSeparator + string.Format(numberFormat, array[0, j]));
            }
            writer.Write(ArrayFormat.RowEnd);

            // Subsequent rows
            for (int i = 1; i < array.GetLength(0); ++i)
            {
                writer.Write(ArrayFormat.RowSeparator + ArrayFormat.RowStart);
                writer.Write(string.Format(numberFormat, array[i, 0]));
                for (int j = 1; j < array.GetLength(1); ++j)
                {
                    writer.Write(ArrayFormat.ColSeparator + string.Format(numberFormat, array[i, j]));
                }
                writer.Write(ArrayFormat.RowEnd);
            }
            writer.Write(ArrayFormat.RowSeparator + ArrayFormat.ArrayEnd);
        }
    }
}
