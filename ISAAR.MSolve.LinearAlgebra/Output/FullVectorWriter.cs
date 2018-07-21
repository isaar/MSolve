using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class FullVectorWriter
    {
        private readonly bool writeLengthFirst = false;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="firstLineIsLength">If true, the length of the vector will be written in the first line and then the 
        ///     vector entries will follow starting from the second.</param>
        public FullVectorWriter(bool firstLineIsLength = false)
        {
            this.writeLengthFirst = firstLineIsLength;
        }

        public Array1DFormat ArrayFormat { get; set; } = Array1DFormat.PlainHorizontal;
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        public void WriteToConsole(IIndexable1D vector)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(vector, writer));
        }

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
