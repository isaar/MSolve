using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace ISAAR.MSolve.LinearAlgebra.Input
{
    /// <summary>
    /// Each row is on a different line
    /// </summary>
    public class Array2DReader
    {
        private readonly char[] colSeparators;
        private readonly bool firstLineIsDimensions;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="firstLineIsDimensions">row, col seperated by ',' or ' '</param>
        /// <param name="rowSeparator"></param>
        /// <param name="colSeparator"></param>
        public Array2DReader(bool firstLineIsDimensions, params char[] colSeparators)
        {
            this.firstLineIsDimensions = firstLineIsDimensions;
            this.colSeparators = colSeparators;
        }

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
                // Read the matrix dimensions
                string firstLine = reader.ReadLine();
                if (firstLine == null) throw new IOException("Empty file");
                string[] subStrings = firstLine.Split(new char[] { ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);
                if (subStrings.Length != 2) throw new IOException(
                    $"There must be 2 dimensions in the first line, but there were {subStrings.Length}");
                int numRows = Int32.Parse(subStrings[0]);
                int numCols = Int32.Parse(subStrings[1]);

                return ReadMatrixWithDimensions(reader, numRows, numCols);
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
                return ReadMatrixWithDimensions(reader, numRows, numCols);
            }
        }

        private double[,] ReadMatrixWithDimensions(StreamReader reader, int numRows, int numCols)
        {
            var matrix = new double[numRows, numCols];
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
