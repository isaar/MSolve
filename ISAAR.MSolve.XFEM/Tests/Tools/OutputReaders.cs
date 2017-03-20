using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.XFEM.Integration.Points;


namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    static class OutputReaders
    {
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

        public static Matrix2D<double>[] ReadElementStiffnessMatrices(string path, int elementsCount)
        {
            string[] lines = File.ReadAllLines(path);
            var elementMatrices = new Matrix2D<double>[elementsCount];

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
                        Matrix2D<double> finishedMatrix = new Matrix2D<double>(rows, cols);
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
    }
}
