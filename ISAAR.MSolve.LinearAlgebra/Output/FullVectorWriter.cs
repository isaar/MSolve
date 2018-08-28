using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes all entries of a vector to Console or a file. The user may choose how to format and justify them.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class FullVectorWriter
    {
        private readonly bool writeLengthFirst = false;

        /// <summary>
        /// Initializes a new instance of the <see cref="FullVectorWriter"/> class with the specified settings.
        /// </summary>
        /// <param name="firstLineIsLength">If true, the length of the vector will be written in the first line and then the 
        ///     vector entries will follow starting from the second line.</param>
        public FullVectorWriter(bool firstLineIsLength = false)
        {
            this.writeLengthFirst = firstLineIsLength;
        }

        /// <summary>
        /// Describes how the vector entries will be separated.
        /// </summary>
        public Array1DFormat ArrayFormat { get; set; } = Array1DFormat.PlainHorizontal;

        /// <summary>
        /// Describes how the vector entries will be formatted and justified.
        /// </summary>
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        /// <summary>
        /// Writes the provided vector to Console.
        /// </summary>
        /// <param name="vector">The vector to write.</param>
        public void WriteToConsole(IIndexable1D vector)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(vector, writer));
        }

        /// <summary>
        /// Writes the provided vector to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="vector">The vector to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="vector"/> will be written.</param>
        /// <param name="append">If true, <paramref name="vector"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(IIndexable1D vector, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(vector, writer), path, append);
        }

        private void WriteToStream(IIndexable1D vector, StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            string separator = ArrayFormat.Separator;
            if (writeLengthFirst) writer.WriteLine(vector.Length);
            writer.Write(ArrayFormat.Start);
            writer.Write(string.Format(numberFormat, vector[0]));
            for (int i = 1; i < vector.Length; ++i)
            {
                writer.Write(separator + string.Format(numberFormat, vector[i]));
            }
            writer.WriteLine(ArrayFormat.End);
        }
    }
}
