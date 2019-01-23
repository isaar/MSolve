using System;
using System.IO;

namespace ISAAR.MSolve.LinearAlgebra.Input
{
    /// <summary>
    /// Reads 2D arrays from files. Each row must be written on a different line of the file.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Array2DReader
    {
        private readonly char[] colSeparators;
        private readonly bool firstLineIsDimensions;

        /// <summary>
        /// Initializes a new instance of the <see cref="Array2DReader"/> class, with the provided settings.
        /// </summary>
        /// <param name="firstLineIsDimensions">If true, the first line will be treated as containing the numbers of rows and 
        ///     columns of the array to be read. The number of rows must be first, then a ' ' or/and ',' and then the number of
        ///     columns. This allows optimizations, but will fail if there is no such line.</param>
        /// <param name="colSeparators">Each row is read from a different line into a string. This string is then split 
        ///     whenever one of the characters in <paramref name="colSeparators"/> is met. If none is provided, then the string 
        ///     will be split at each whitespace character.</param>
        public Array2DReader(bool firstLineIsDimensions, params char[] colSeparators)
        {
            this.firstLineIsDimensions = firstLineIsDimensions;
            this.colSeparators = colSeparators;
        }

        /// <summary>
        /// Reads the file at <paramref name="path"/> and returns a 2D array.
        /// </summary>
        /// <param name="path">The absolute path of the array.</param>
        /// <exception cref="IOException">Thrown if there is no such file or if the dimensions of the 2D array specified in the 
        ///     first line are invalid.</exception>
        public double[,] ReadFile(string path)
        {
            //TODO: add input checking
            if (firstLineIsDimensions) return ReadKnownSize(path);
            else return ReadUnknownSize(path);
        }

        private double[,] ReadKnownSize(string path)
        {
            using (var reader = new StreamReader(path))
            {
                // Read the array dimensions
                string firstLine = reader.ReadLine();
                if (firstLine == null) throw new IOException("Empty file");
                string[] subStrings = firstLine.Split(new char[] { ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);
                if (subStrings.Length != 2) throw new IOException(
                    $"There must be 2 dimensions in the first line, but there were {subStrings.Length}");
                int numRows = Int32.Parse(subStrings[0]);
                int numCols = Int32.Parse(subStrings[1]);

                return ReadArrayWithDimensions(reader, numRows, numCols);
            }
        }

        private double[,] ReadUnknownSize(string path)
        {
            int numRows, numCols;
            using (var reader = new StreamReader(path))
            {
                // Find the columns from the first row
                string line = reader.ReadLine();
                numCols = line.Split(colSeparators, StringSplitOptions.RemoveEmptyEntries).Length;
                numRows = 1;

                // Process the rest of the file to find the number of rows
                while (true)
                {
                    line = reader.ReadLine();
                    if (line == null) break;
                    else if (line.Length > 0) ++numRows;
                }
            }

            using (var reader = new StreamReader(path)) //TODO: find an easy and efficient way to reset the reader
            {
                return ReadArrayWithDimensions(reader, numRows, numCols);
            }
        }

        private double[,] ReadArrayWithDimensions(StreamReader reader, int numRows, int numCols)
        {
            var array = new double[numRows, numCols];
            for (int i = 0; i < numRows; ++i)
            {
                string[] colStrings = reader.ReadLine().Split(colSeparators, StringSplitOptions.RemoveEmptyEntries);
                if (colStrings.Length != numCols) throw new IOException(
                    $"There must be {numCols} rows, but {colStrings.Length} were found.");
                for (int j = 0; j < numCols; ++j) array[i, j] = Double.Parse(colStrings[j]);
            }
            return array;
        }
    }
}
