using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;

namespace ISAAR.MSolve.LinearAlgebra.Input
{
    public class CoordinateTextFileReader
    {
        private static readonly char[] separators = { };

        private readonly LinkedList<(int row, int col, double val)> entries;
        private int numRows, numCols, nnz;

        public CoordinateTextFileReader()
        {
            this.entries = new LinkedList<(int row, int col, double val)>();
        }

        public void ReadFromFile(string path)
        {
            using (var reader = new StreamReader(path))
            {
                string firstLine = reader.ReadLine();
                if (firstLine == null) throw new IOException("Empty file");
                (numRows, numCols, nnz) = ParseHeader(firstLine);
                for (int i = 0; i < nnz; ++i)
                {
                    string line = reader.ReadLine();
                    if (line == null) throw new IOException("Fewer lines than declared in header");
                    entries.AddLast(ParseEntry(line));
                }
            }
        }

        /// <summary>
        /// This will also remove the matrix saved in this object, to avoid huge memory consumption. 
        /// </summary>
        /// <returns></returns>
        public DOKSymmetricColMajor ToSymmetricDOK()
        {
            var dok = DOKSymmetricColMajor.CreateEmpty(numRows);
            for (int i = 0; i < nnz; ++i)
            {
                (int row, int col, double val) = entries.First.Value;
                dok[row, col] = val;
                entries.RemoveFirst();
            }
            return dok;
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
    }
}
