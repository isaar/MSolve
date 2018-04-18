using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.XFEM.Integration.Points;


namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    static class OutputReaders
    {
        /// <summary>
        /// Expected file format:
        /// * 0
        /// xi0_0 eta0_0 weight0_0
        /// xi0_1 eta0_1 weight0_1
        /// ...
        /// * 1
        /// xi1_0 eta1_0 weight1_0
        /// xi1_1 eta1_1 weight1_1
        /// ...
        /// *
        /// </summary>
        /// <param name="path"></param>
        /// <param name="elementsCount"></param>
        /// <returns></returns>
        public static GaussPoint2D[][] ReadAllGaussPoints(string path, int elementsCount)
        {
            string[] lines = File.ReadAllLines(path);
            var allPoints = new GaussPoint2D[elementsCount][];

            const int START = -1;
            int currentElement = START;
            List<GaussPoint2D> currentGPs = null;
            for (int i = 0; i < lines.Length; ++i)
            {
                string line = lines[i];
                string[] words = line.Split(' ');
                if (line[0] == '*') // New matrix
                {
                    // Build and store the finished matrix
                    if (currentElement != START)
                    {
                        allPoints[currentElement] = currentGPs.ToArray();
                    }
                    if (words.Length > 1) // The last line only has a *.
                    {
                        currentElement = int.Parse(words[1]);
                        currentGPs = new List<GaussPoint2D>();
                    }
                }
                else
                {
                    // Read a new Gauss point of the form: xi eta weight
                    if (words.Length != 3)
                        throw new IOException("Line " + i + ": A line containing a gauss point must have 3 values only");
                    double xi = double.Parse(words[0]);
                    double eta = double.Parse(words[1]);
                    double weight = double.Parse(words[2]);
                    currentGPs.Add(new GaussPoint2D(xi, eta, weight));
                }
            }
            return allPoints;
        }

        /// <summary>
        /// Expected file format:
        /// * 0
        /// row0OfK0
        /// row1OfK0
        /// ...
        /// * 1
        /// row0OfK1
        /// row1OfK1
        /// ...
        /// *
        /// Where each row is a series of floating point numbers, seperated by single spaces. 
        /// </summary>
        /// <param name="path"></param>
        /// <param name="elementsCount"></param>
        /// <returns></returns>
        public static Matrix[] ReadElementStiffnessMatrices(string path, int elementsCount)
        {
            string[] lines = File.ReadAllLines(path);
            var elementMatrices = new Matrix[elementsCount];

            const int START = -1;
            int currentElement = START;
            List<List<double>> currentMatrix = null;
            foreach (string line in lines)
            {
                string[] words = line.Split(' ');
                if (line[0] == '*') // New matrix
                {
                    // Build and store the finished matrix
                    if (currentElement != START)
                    {
                        int rows = currentMatrix.Count;
                        int cols = currentMatrix[0].Count;
                        Matrix finishedMatrix = Matrix.CreateZero(rows, cols);
                        for (int r = 0; r < rows; ++r)
                        {
                            for (int c = 0; c < cols; ++c)
                            {
                                finishedMatrix[r, c] = currentMatrix[r][c];
                            }
                        }
                        elementMatrices[currentElement] = finishedMatrix;
                    }
                    if (words.Length > 1) // The last line only has a *.
                    {
                        currentElement = int.Parse(words[1]);
                        currentMatrix = new List<List<double>>();
                    }
                }
                else
                {
                    // Read a new row of doubles
                    var newRow = new List<double>();
                    foreach (string number in words)
                    {
                        newRow.Add(double.Parse(number));
                    }
                    currentMatrix.Add(newRow);
                }
            }
            return elementMatrices;
        }

        /// <summary>
        /// Expected file format: DENSE MATRIX:
        /// K00 K01 ... K0N
        /// K10 K11 ... K1N
        /// ...
        /// KM0 KM1 ... KMN
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static Matrix ReadGlobalStiffnessMatrix(string path)
        {
            string[] lines = File.ReadAllLines(path);
            int order = lines.Length;
            string lastLine = lines[lines.Length - 1];
            if (string.IsNullOrEmpty(lastLine) || string.IsNullOrWhiteSpace(lastLine)) --order; //Last row might have output an extra newline

            var matrix = Matrix.CreateZero(order, order);
            for (int row = 0; row < order; ++row)
            {
                string[] words = lines[row].Split(' ');
                for (int col = 0; col < order; ++col)
                {
                    matrix[row, col] = double.Parse(words[col]);
                }
            }

            return matrix;
        }

        /// <summary>
        /// Expected file format: The i-th row contains the dofs of the i-th node, seperated by single spaces.
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static int[][] ReadNodalDofs(string path)
        {
            string[] lines = File.ReadAllLines(path);
            int nodesCount = lines.Length;

            string lastLine = lines[lines.Length - 1];
            if (string.IsNullOrEmpty(lastLine) || string.IsNullOrWhiteSpace(lastLine)) --nodesCount;

            int[][] nodalDofs = new int[nodesCount][];
            for (int node = 0; node < nodesCount; ++node)
            {
                var dofs = new List<int>();
                foreach (var word in lines[node].Split(' '))
                {
                    int zeroBased = int.Parse(word) - 1;
                    if (zeroBased >= 0) dofs.Add(zeroBased);
                }
                nodalDofs[node] = dofs.ToArray();
            }

            return nodalDofs;
        }
    }
}
