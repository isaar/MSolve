using System;
using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    /// <summary>
    /// Writes the internal indexing and value arrays of a sparse matrix to Console or a file. The user may choose how to format 
    /// and justify their entries.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class RawArraysWriter
    {
        private readonly bool lineBetweenArrays;
        private readonly bool titlesOnOtherLines;

        /// <summary>
        /// Initializes a new instance of the <see cref="FullVectorWriter"/> class with the specified settings.
        /// </summary>
        /// <param name="lineBetweenArrays">If true, the arrays will be separated by an empty line. If false, they will be 
        ///     separated only by a newline character.</param>
        /// <param name="titlesOnOtherLines">Each array is preceded by its title. If <paramref name="titlesOnOtherLines"/> is
        ///     true, then a newline character will be added after each title. If false, the titles will be in the same lines
        ///     as the respective arrays.</param>
        public RawArraysWriter(bool lineBetweenArrays = true, bool titlesOnOtherLines = true)
        {
            this.lineBetweenArrays = lineBetweenArrays;
            this.titlesOnOtherLines = titlesOnOtherLines;
        }

        /// <summary>
        /// Describes how the array entries will be formatted and justified.
        /// </summary>
        public INumericFormat NumericFormat { get; set; } = new GeneralNumericFormat();

        /// <summary>
        /// Use Environment.NewLine (the default) if you plan to read the file using <see cref="Input.RawArraysReader"/>  
        /// </summary>
        public string EntrySeparator { get; set; } = Environment.NewLine;

        /// <summary>
        /// Writes the internal arrays of the provided sparse matrix to Console.
        /// </summary>
        /// <param name="matrix">The sparse matrix to write.</param>
        public void WriteToConsole(ISparseMatrix matrix)
        {
            Utilities.WriteToConsole((writer) => WriteToStream(matrix, writer));
        }

        /// <summary>
        /// Writes the internal arrays of the provided sparse matrix to the file at <paramref name="path"/>.
        /// </summary>
        /// <param name="matrix">The sparse matrix to write.</param>
        /// <param name="path">The absolute path of the file, where <paramref name="matrix"/> will be written.</param>
        /// <param name="append">If true, <paramref name="matrix"/> will be written after the current contents of the file at
        ///     <paramref name="path"/>. If false, it will overwrite them.</param>
        public void WriteToFile(ISparseMatrix matrix, string path, bool append = false)
        {
            Utilities.WriteToFile((writer) => WriteToStream(matrix, writer), path, append);
        }

        /// <summary>
        /// Writes each internal array of the provided sparse matrix to a different file. All these filenames have a common 
        /// prefix, extracted from <paramref name="pathBase"/>, and a suffix corresponding to the purpose of each array.
        /// </summary>
        /// <param name="matrix">The sparse matrix to write.</param>
        /// <param name="pathBase">An absolute path. This filename will not be used directly. Instead one file per internal array
        ///     of <paramref name="matrix"/> will be used. Each file will be suffixed with the name of that array and then the 
        ///     extension specified in <paramref name="pathBase"/>.</param>
        /// <param name="writeArrayLengthFirst">If true, the first line of each file will contain the length of the corresponding
        ///     array. If false, there will be only one line per file, which will contain the array.</param>
        public void WriteToMultipleFiles(ISparseMatrix matrix, string pathBase, bool writeArrayLengthFirst = true) // TODO: this should be a different writer
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
                if (lineBetweenArrays) writer.WriteLine();
                writer.WriteLine(); // otherwise everything would be on the same line
                writer.Write(nameArrayPair.Key + ": ");
                if (titlesOnOtherLines) writer.WriteLine();
                WriteArray(nameArrayPair.Value, writer, false);
            }
        }

        private void WriteArray(IReadOnlyList<int> array, StreamWriter writer, bool writeArrayLengthFirst)
        {
            if (writeArrayLengthFirst) writer.WriteLine(array.Count);
            int last = array.Count - 1;
            for (int i = 0; i < last; ++i)
            {
                writer.Write(array[i] + EntrySeparator);
            }
            writer.Write(array[last]);
        }

        private void WriteArray(IReadOnlyList<double> array, StreamWriter writer, bool writeArrayLengthFirst)
        {
            if (writeArrayLengthFirst) writer.WriteLine(array.Count);
            string numberFormat = NumericFormat.GetRealNumberFormat();
            int last = array.Count - 1;
            for (int i = 0; i < last; ++i)
            {
                writer.Write(string.Format(numberFormat, array[i]) + EntrySeparator);
            }
            writer.Write(array[last]);
        }
    }
}
