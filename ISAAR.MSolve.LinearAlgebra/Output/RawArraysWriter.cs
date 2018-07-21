using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class RawArraysWriter
    {
        private readonly bool spaceBetweenArrays;
        private readonly bool titlesOnOtherLines;

        public RawArraysWriter(bool spaceBetweenArrays = true, bool titlesOnOtherLines = true)
        {
            this.spaceBetweenArrays = spaceBetweenArrays;
            this.titlesOnOtherLines = titlesOnOtherLines;
        }

        public INumericFormat NumericFormat { get; set; } = new GeneralNumericFormat();

        public void WriteToConsole(ISparseMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(matrix, writer));
        }

        public void WriteToFile(ISparseMatrix matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(matrix, writer), path, append);
        }

        /// <summary>
        /// Each internal array is written to a different file. All these files have a common prefix, chosen by the user, and 
        /// a suffix corresponding to the purpose of each array.
        /// </summary>
        public void WriteToMultipleFiles(ISparseMatrix matrix, string pathBase, bool writeArrayLengthFirst = true)
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
                WriteArray(sparseFormat.RawValuesArray, writer, writeArrayLengthFirst);
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
                    WriteArray(nameArrayPair.Value, writer, writeArrayLengthFirst);
                }
            }
        }

        private void WriteToStream(ISparseMatrix matrix, StreamWriter writer)
        {
            SparseFormat sparseFormat = matrix.GetSparseFormat();
            writer.Write(sparseFormat.RawValuesTitle + ": ");
            if (titlesOnOtherLines) writer.WriteLine();
            WriteArray(sparseFormat.RawValuesArray, writer, false);

            foreach (var nameArrayPair in sparseFormat.RawIndexArrays)
            {
                if (spaceBetweenArrays) writer.WriteLine();
                writer.WriteLine(); // otherwise everything would be on the same line
                writer.Write(nameArrayPair.Key + ": ");
                if (titlesOnOtherLines) writer.WriteLine();
                WriteArray(nameArrayPair.Value, writer, false);
            }
        }

        private void WriteArray(IReadOnlyList<int> array, StreamWriter writer, bool writeArrayLengthFirst)
        {
            if (writeArrayLengthFirst) writer.Write(array.Count + " ");
            int last = array.Count - 1;
            for (int i = 0; i < last; ++i)
            {
                writer.Write(array[i] + " ");
            }
            writer.Write(array[last]);
        }

        private void WriteArray(IReadOnlyList<double> array, StreamWriter writer, bool writeArrayLengthFirst)
        {
            if (writeArrayLengthFirst) writer.Write(array.Count + " ");
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
