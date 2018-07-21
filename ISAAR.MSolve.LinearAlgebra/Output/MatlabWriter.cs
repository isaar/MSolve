using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class MatlabWriter
    {
        public MatlabWriter()
        {
        }

        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 12 };

        public void WriteToConsole(IIndexable1D vector)
        {
            Utilities.WriteToConsole((writer) => WriteFullVector(vector, writer));
        }

        public void WriteToConsole(ISparseMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteSparseMatrix(matrix, writer));
        }

        public void WriteToFile(IIndexable1D vector, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteFullVector(vector, writer), path, append);
        }

        public void WriteToFile(ISparseMatrix matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteSparseMatrix(matrix, writer), path, append);
        }

        private void WriteFullVector(IIndexable1D vector, StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            for (int i = 0; i < vector.Length - 1; ++i) writer.WriteLine(string.Format(numberFormat, vector[i]));
            writer.Write(string.Format(numberFormat, vector[vector.Length - 1]));
        }

        private void WriteSparseMatrix(ISparseMatrix matrix, StreamWriter writer)
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            foreach (var (row, col, val) in matrix.EnumerateNonZeros())
            {
                writer.Write($"{row + 1} {col + 1} ");
                writer.WriteLine(string.Format(numberFormat, val));
            }
            writer.Write($"{matrix.NumRows} {matrix.NumColumns} 0.0");
        }
    }
}
