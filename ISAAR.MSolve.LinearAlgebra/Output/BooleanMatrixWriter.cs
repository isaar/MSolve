using System;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes all entries of a <see cref="SignedBooleanMatrix"/> matrix to Console or a file.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class BooleanMatrixWriter
    {
        private readonly string numberFormat;

        /// <summary>
        /// Initializes a new instance of the <see cref="BooleanMatrixWriter"/> class with the specified settings.
        /// </summary>
        /// <param name="signsForOnes">If true, the output matrices will be assumed to contain +1, -1 and 0, which affects the 
        ///     justification during formatting. Otherwise the output matrices will be assumed to contain 0 and 1 only.</param>
        public BooleanMatrixWriter(bool signsForOnes = true)
        {
            if (signsForOnes) numberFormat = "{0,3}";
            else numberFormat = "{0,2}";
        }

        /// <summary>
        /// Writes the provided matrix to Console.
        /// </summary>
        /// <param name="matrix">The matrix to write.</param>
        public void WriteToConsole(SignedBooleanMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(matrix, writer));
        }

        /// <summary>
        /// Writes the provided matrix to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="matrix">The matrix to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="matrix"/> will be written.</param>
        /// <param name="append">If true, <paramref name="matrix"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(SignedBooleanMatrix matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(matrix, writer), path, append);
        }

        private void WriteToStream(SignedBooleanMatrix matrix, StreamWriter writer)
        {
            for (int i = 0; i < matrix.NumRows; ++i)
            {
                for (int j = 0; j < matrix.NumColumns; ++j)
                {
                    writer.Write(string.Format(numberFormat, matrix[i, j]));
                }
                writer.WriteLine();
            }
            writer.WriteLine();
        }
    }
}
