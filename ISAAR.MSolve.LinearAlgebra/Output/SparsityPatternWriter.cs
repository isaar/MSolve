using System.IO;
using ISAAR.MSolve.LinearAlgebra.Reordering;

//TODO: perhaps rename it to BooleanMatrixWriter and accept BooleanMatrix, ISparsityPattern and ISparseMatrix
namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes a boolean matrix with 1 corresponding to non zero entries and 0 corresponding to zero entries of a sparse matrix.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SparsityPatternWriter
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="FullVectorWriter"/> class.
        /// </summary>
        public SparsityPatternWriter()
        {
        }

        /// <summary>
        /// Writes the sparsity pattern of a matrix to Console.
        /// </summary>
        /// <param name="matrix">The sparsity pattern of the matrix to write.</param>
        public void WriteToConsole(ISparsityPattern matrix)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(matrix, writer));
        }

        /// <summary>
        /// Writes the sparsity pattern of a matrix to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="matrix">The sparsity pattern of the matrix to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="matrix"/> will be written.</param>
        /// <param name="append">If true, <paramref name="matrix"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(ISparsityPattern matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(matrix, writer), path, append);
        }

        private void WriteToStream(ISparsityPattern pattern, StreamWriter writer)
        {
            for (int i = 0; i < pattern.NumRows; ++i)
            {
                for (int j = 0; j < pattern.NumColumns; ++j)
                {
                    if (pattern.IsNonZero(i, j)) writer.Write("1 ");
                    else writer.Write("0 ");
                }
                writer.WriteLine();
            }
        }
    }
}
