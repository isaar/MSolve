using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Reordering;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes the a boolean matrix with 1 corresponding to non zero entries.
    /// </summary>
    class SparsityPatternWriter: MatrixWriter
    {
        private readonly ISparsityPattern pattern;

        public SparsityPatternWriter(ISparsityPattern pattern)
        {
            this.pattern = pattern;
        }

        protected override void WriteToStream(StreamWriter writer)
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
