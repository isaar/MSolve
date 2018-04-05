using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    public class RawArraysWriter: MatrixWriter
    {
        private readonly ISparseMatrix matrix;
        private readonly bool spaceBetweenArrays;
        private readonly bool titlesOnOtherLines;

        public RawArraysWriter(ISparseMatrix matrix, bool spaceBetweenArrays = true, bool titlesOnOtherLines = true)
        {
            this.matrix = matrix;
            this.spaceBetweenArrays = spaceBetweenArrays;
            this.titlesOnOtherLines = titlesOnOtherLines;
        }

        protected override void WriteToStream(StreamWriter writer)
        {
            SparseFormat sparseFormat = matrix.GetSparseFormat();
            writer.Write(sparseFormat.RawValuesTitle + ": ");
            if (titlesOnOtherLines) writer.WriteLine();
            WriteArray<double>(sparseFormat.RawValuesArray, writer);

            foreach (KeyValuePair<string, IReadOnlyList<int>> pair in sparseFormat.RawIndexArrays)
            {
                if (spaceBetweenArrays) writer.WriteLine();
                writer.WriteLine(); // otherwise everything would be on the same line
                writer.Write(pair.Key + ": ");
                if (titlesOnOtherLines) writer.WriteLine();
                WriteArray<int>(pair.Value, writer);
            }
        }

        private static void WriteArray<T>(IReadOnlyList<T> array, StreamWriter writer) //TODO: perhaps move it to abstract class
        {
            int last = array.Count - 1;
            for (int i = 0; i < last; ++i)
            {
                writer.Write(array[i] + " ");
            }
            writer.Write(array[last]);
        }
    }
}
