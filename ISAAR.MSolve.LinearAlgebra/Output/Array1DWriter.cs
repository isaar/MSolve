using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes all entries of a 1D array to Console or a file. The user may choose how to format and justify them.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Array1DWriter
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Array1DWriter"/> class.
        /// </summary>
        public Array1DWriter()
        {
        }

        /// <summary>
        /// Describes how the array entries will be separated.
        /// </summary>
        public Array1DFormat ArrayFormat { get; set; } = Array1DFormat.PlainHorizontal;

        /// <summary>
        /// Describes how the array entries will be formatted and justified.
        /// </summary>
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        /// <summary>
        /// Writes the provided array to Console.
        /// </summary>
        /// <param name="array">The array to write.</param>
        public void WriteToConsole(double[] array)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(array, writer));
        }

        /// <summary>
        /// Writes the provided array to Console.
        /// </summary>
        /// <param name="array">The integer array to write.</param>
        public void WriteToConsole(int[] array)
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
        public void WriteToFile(double[] array, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(array, writer), path, append);
        }

        /// <summary>
        /// Writes the provided array to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="array">The integer array to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="array"/> will be written.</param>
        /// <param name="append">If true, <paramref name="array"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(int[] array, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(array, writer), path, append);
        }

        private void WriteToStream(double[] array, StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            string separator = ArrayFormat.Separator;
            writer.Write(ArrayFormat.Start);
            writer.Write(string.Format(numberFormat, array[0]));
            for (int i = 1; i < array.Length; ++i)
            {
                writer.Write(separator + string.Format(numberFormat, array[i]));
            }
            writer.WriteLine(ArrayFormat.End);
        }

        private void WriteToStream(int[] array, StreamWriter writer)
        {
            string separator = ArrayFormat.Separator;
            writer.Write(ArrayFormat.Start);
            writer.Write(array[0]);
            for (int i = 1; i < array.Length; ++i)
            {
                writer.Write(separator + array[i]);
            }
            writer.WriteLine(ArrayFormat.End);
        }
    }
}
