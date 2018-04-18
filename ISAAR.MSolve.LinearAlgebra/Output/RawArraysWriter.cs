using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.LinearAlgebra.Output
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

        public static INumericFormat NumericFormat { get; set; } = new GeneralNumericFormat();

        /// <summary>
        /// Each internal array is written to a different file. All these files have a common prefix, chosen by the user, and 
        /// a suffix corresponding to the purpose of each array.
        /// </summary>
        public void WriteToMultipleFiles(string pathBase)
        {
            string path = Path.GetDirectoryName(pathBase);
            string nameOnly = Path.GetFileNameWithoutExtension(pathBase);
            string ext = Path.GetExtension(pathBase);
            SparseFormat sparseFormat = matrix.GetSparseFormat();

            // Values array
            string suffix = "-" + sparseFormat.RawValuesTitle.ToLower();
            string valuesPath = path + "\\" + nameOnly + suffix + ext; // Not too sure about the \\
            using (var writer = new StreamWriter(valuesPath))
            {
#if DEBUG
                writer.AutoFlush = true; // To look at intermediate output at certain breakpoints
#endif
                WriteArray(sparseFormat.RawValuesArray, writer);
            }

            // Indexing arrays
            foreach (var nameArrayPair in sparseFormat.RawIndexArrays)
            {
                suffix = "-" + nameArrayPair.Key.ToLower();
                string indexerPath = path + "\\" + nameOnly + suffix + ext; // Not too sure about the \\
                using (var writer = new StreamWriter(indexerPath))
                {
#if DEBUG
                    writer.AutoFlush = true; // To look at intermediate output at certain breakpoints
#endif
                    WriteArray(nameArrayPair.Value, writer);
                }
            }
        }

        protected override void WriteToStream(StreamWriter writer)
        {
            SparseFormat sparseFormat = matrix.GetSparseFormat();
            writer.Write(sparseFormat.RawValuesTitle + ": ");
            if (titlesOnOtherLines) writer.WriteLine();
            WriteArray(sparseFormat.RawValuesArray, writer);

            foreach (var nameArrayPair in sparseFormat.RawIndexArrays)
            {
                if (spaceBetweenArrays) writer.WriteLine();
                writer.WriteLine(); // otherwise everything would be on the same line
                writer.Write(nameArrayPair.Key + ": ");
                if (titlesOnOtherLines) writer.WriteLine();
                WriteArray(nameArrayPair.Value, writer);
            }
        }

        private static void WriteArray(IReadOnlyList<int> array, StreamWriter writer) //TODO: perhaps move it to abstract class
        {
            int last = array.Count - 1;
            for (int i = 0; i < last; ++i)
            {
                writer.Write(array[i] + " ");
            }
            writer.Write(array[last]);
        }

        private static void WriteArray(IReadOnlyList<double> array, StreamWriter writer) //TODO: perhaps move it to abstract class
        {
            string numberFormat = NumericFormat.GetRealNumberFormat();
            int last = array.Count - 1;
            for (int i = 0; i < last; ++i)
            {
                writer.Write(string.Format(numberFormat, array[i]) + " ");
            }
            writer.Write(array[last]);
        }
    }
}
