using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class FullMatrixWriter
    {
        public FullMatrixWriter()
        {
        }

        public Array2DFormat ArrayFormat { get; set; } = Array2DFormat.Plain;
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 6 };

        public void WriteToConsole(IIndexable2D matrix)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(matrix, writer));
        }

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
            writer.Write(ArrayFormat.RowEnd);

            // Subsequent rows
            for (int i = 1; i < matrix.NumRows; ++i)
            {
                writer.Write(ArrayFormat.RowSeparator + ArrayFormat.RowStart);
                writer.Write(string.Format(numberFormat, matrix[i, 0]));
                for (int j = 1; j < matrix.NumColumns; ++j)
                {
                    writer.Write(ArrayFormat.ColSeparator + string.Format(numberFormat, matrix[i, j]));
                }
                writer.Write(ArrayFormat.RowEnd);
            }
            writer.Write(ArrayFormat.RowSeparator + ArrayFormat.ArrayEnd);
        }
    }
}
