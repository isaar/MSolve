using System;
using System.Collections.Generic;
using System.IO;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;

//TODO: read methods with the output class name are bad design.
//TODO: be more flexible with the format by setting options in the constructor.
namespace ISAAR.MSolve.LinearAlgebra.Input
{
    /// <summary>
    /// Reads matrices that are written in files with coordinate format. In this format, the first line contains: 
    /// numRows numCols numNonZeros. The previous '.' must not be included in the line. 
    /// Each subsequent line contains a non-zero matrix entry defined as a triplet: rowIndex columnIndex value. 
    /// The previous '.' must not be included in the line.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CoordinateTextFileReader
    {
        private static readonly char[] separators = { };

        /// <summary>
        /// Initializes a new instance of the <see cref="CoordinateTextFileReader"/> class.
        /// </summary>
        public CoordinateTextFileReader()
        {
        }

        /// <summary>
        /// Reads a general sparse matrix from the file specified by <paramref name="path"/>. Even if the matrix is symmetric,
        /// the file must contain both the upper and lower triangle.
        /// </summary>
        /// <param name="path">The absolute path of the array.</param>
        /// <exception cref="IOException">Thrown if there is no such file or if the dimensions of the matrix specified in the
        ///     first line are invalid.</exception>
        public DokColMajor ReadFileAsDokColMajor(string path)
        {
            using (var reader = new StreamReader(path))
            {
                (int numRows, int numCols, int nnz) = ReadMatrixDimensions(reader);
                var builder = DokColMajor.CreateEmpty(numRows, numCols);
                ReadMatrixEntries(reader, nnz, builder);
                return builder;
            }
        }

        /// <summary>
        /// Reads a general sparse matrix from the file specified by <paramref name="path"/>. Even if the matrix is symmetric,
        /// the file must contain both the upper and lower triangle.
        /// </summary>
        /// <param name="path">The absolute path of the array.</param>
        /// <exception cref="IOException">Thrown if there is no such file or if the dimensions of the matrix specified in the
        ///     first line are invalid.</exception>
        public DokRowMajor ReadFileAsDokRowMajor(string path)
        {
            using (var reader = new StreamReader(path))
            {
                (int numRows, int numCols, int nnz) = ReadMatrixDimensions(reader);
                var builder = DokRowMajor.CreateEmpty(numRows, numCols);
                ReadMatrixEntries(reader, nnz, builder);
                return builder;
            }
        }

        /// <summary>
        /// Reads a symmetric sparse matrix from the file specified by <paramref name="path"/>. The file needs to contain only 
        /// the upper or lower triangle.
        /// </summary>
        /// <param name="path">The absolute path of the array.</param>
        /// <exception cref="IOException">Thrown if there is no such file or if the dimensions of the matrix specified in the
        ///     first line are invalid.</exception>
        public DokSymmetric ReadFileAsDokSymmetricColMajor(string path)
        {
            using (var reader = new StreamReader(path))
            {
                (int numRows, int numCols, int nnz) = ReadMatrixDimensions(reader);
                var builder = DokSymmetric.CreateEmpty(numCols);
                ReadMatrixEntries(reader, nnz, builder);
                return builder;
            }
        }

        private static (int numRows, int numCols, int numNonzeros) ParseHeader(string line)
        {
            //TODO: add error checking
            string[] subStrings = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
            int numRows = Int32.Parse(subStrings[0]);
            int numCols = Int32.Parse(subStrings[1]);
            int nnz = Int32.Parse(subStrings[2]);
            return (numRows, numCols, nnz);
        }

        private static (int row, int col, double val) ParseEntry(string line)
        {
            //TODO: add error checking
            string[] subStrings = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
            int row = Int32.Parse(subStrings[0]);
            int col = Int32.Parse(subStrings[1]);
            double val = Double.Parse(subStrings[2]);
            return (row, col, val);
        }

        private static (int numRows, int numCols, int numNonZeros) ReadMatrixDimensions(StreamReader reader)
        {
            string firstLine = reader.ReadLine();
            if (firstLine == null) throw new IOException("Empty file");
            return ParseHeader(firstLine);
        }

        private static void ReadMatrixEntries(StreamReader reader, int nnz, IMatrixBuilder builder)
        {
            for (int i = 0; i < nnz; ++i)
            {
                string line = reader.ReadLine();
                if (line == null) throw new IOException("Fewer lines than declared in header");
                (int rowIdx, int colIdx, double value) = ParseEntry(line);
                builder[rowIdx, colIdx] = value;
            }
        }
    }
}
