using System;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: reduce code duplication between this and Array2DReader.
namespace ISAAR.MSolve.LinearAlgebra.Input
{
    /// <summary>
    /// Reads full (as opposed to sparse) matrices from files. Each row must be written on a different line of the file.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class FullMatrixReader
    {
        private readonly char[] colSeparators;
        private readonly bool firstLineIsDimensions;

        /// <summary>
        /// Initializes a new instance of the <see cref="FullMatrixReader"/> class, with the provided settings.
        /// </summary>
        /// <param name="firstLineIsDimensions">
        /// If true, the first line will be treated as containing the numbers of rows and columns of the array to be read. 
        /// The number of rows must be first, then a ' ' or/and ',' and then the number of columns. This allows optimizations, 
        /// but will fail if there is no such line.</param>
        /// <param name="colSeparators">
        /// Each row is read from a different line into a string. This string is then split whenever one of the characters in 
        /// <paramref name="colSeparators"/> is met. If none is provided, then the string will be split at each whitespace 
        /// character.
        /// </param>
        public FullMatrixReader(bool firstLineIsDimensions, params char[] colSeparators)
        {
            this.firstLineIsDimensions = firstLineIsDimensions;
            this.colSeparators = colSeparators;
        }

        /// <summary>
        /// Reads the file at <paramref name="path"/> and returns a matrix.
        /// </summary>
        /// <param name="path">The absolute path of the array.</param>
        /// <exception cref="IOException">Thrown if there is no such file or if the dimensions of the 2D array specified in the 
        ///     first line are invalid.</exception>
        public Matrix ReadFile(string path)
        {
            //TODO: add input checking
            if (firstLineIsDimensions) return ReadKnownSize(path);
            else return ReadUnknownSize(path);
        }

        private Matrix ReadKnownSize(string path)
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

        private Matrix ReadUnknownSize(string path)
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

        private Matrix ReadArrayWithDimensions(StreamReader reader, int numRows, int numCols)
        {
            var matrix = Matrix.CreateZero(numRows, numCols);
            for (int i = 0; i < numRows; ++i)
            {
                string[] colStrings = reader.ReadLine().Split(colSeparators, StringSplitOptions.RemoveEmptyEntries);
                if (colStrings.Length != numCols) throw new IOException(
                    $"There must be {numCols} rows, but {colStrings.Length} were found.");
                for (int j = 0; j < numCols; ++j) matrix[i, j] = Double.Parse(colStrings[j]);
            }
            return matrix;
        }
    }
}
