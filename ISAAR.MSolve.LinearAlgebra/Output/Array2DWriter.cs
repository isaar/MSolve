using System.IO;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class Array2DWriter
    {
        public Array2DWriter()
        {
        }

        public Array2DFormat ArrayFormat { get; set; } = Array2DFormat.Plain;
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        public void WriteToConsole(double[,] array)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(array, writer));
        }

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
