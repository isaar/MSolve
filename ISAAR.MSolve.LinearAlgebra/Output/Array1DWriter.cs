using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class Array1DWriter
    {
        public Array1DWriter()
        {
        }

        public Array1DFormat ArrayFormat { get; set; } = Array1DFormat.PlainHorizontal;
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        public void WriteToConsole(double[] array)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(array, writer));
        }

        public void WriteToFile(double[] array, string path, bool append = false)
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
    }
}
