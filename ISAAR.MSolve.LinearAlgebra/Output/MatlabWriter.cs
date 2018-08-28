using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: separate classes for vectors, dense matrices and sparse matrices.
namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes the entries of a vector or matrix to Console or a file, following the formats used by Matlab.
    /// Auhors: Serafeim Bakalakos
    /// </summary>
    public class MatlabWriter
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="MatlabWriter"/> class.
        /// </summary>
        public MatlabWriter()
        {
        }

        /// <summary>
        /// Describes how the vector and matrix entries will be formatted and justified.
        /// </summary>
        public INumericFormat NumericFormat { get; set; } = new ExponentialFormat { NumDecimalDigits = 12 };

        /// <summary>
        /// Writes the provided vector to Console.
        /// </summary>
        /// <param name="vector">The vector to write.</param>
        public void WriteToConsole(IIndexable1D vector)
        {
            Utilities.WriteToConsole((writer) => WriteFullVector(vector, writer));
        }

        /// <summary>
        /// Writes the provided sparse matrix to Console.
        /// </summary>
        /// <param name="matrix">The sparse matrix to write.</param>
        public void WriteToConsole(ISparseMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteSparseMatrix(matrix, writer));
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
            Utilities.WriteToFile((writer) => WriteFullVector(vector, writer), path, append);
        }

        /// <summary>
        /// Writes the provided sparse matrix to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="matrix">The sparse matrix to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="matrix"/> will be written.</param>
        /// <param name="append">If true, <paramref name="matrix"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
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
