using System.IO;
using ISAAR.MSolve.LinearAlgebra.Reordering;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes the a boolean matrix with 1 corresponding to non zero entries.
    /// </summary>
    public class SparsityPatternWriter
    {
        public SparsityPatternWriter()
        {
        }

        public void WriteToConsole(ISparsityPattern matrix)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(matrix, writer));
        }

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
