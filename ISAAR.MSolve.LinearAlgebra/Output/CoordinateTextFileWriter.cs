using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class CoordinateTextFileWriter
    {
        public CoordinateTextFileWriter()
        {
        }

        public INumericFormat NumericFormat { get; set; } = new GeneralNumericFormat();

        public void WriteToConsole(ISparseMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteMatrix(matrix, writer));
        }

        public void WriteToConsole(ISparseSymmetricMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteSymmetricMatrix(matrix, writer));
        }

        public void WriteToFile(ISparseMatrix matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteMatrix(matrix, writer), path, append);
        }

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
