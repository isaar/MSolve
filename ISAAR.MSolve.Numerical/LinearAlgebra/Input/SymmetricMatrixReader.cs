using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Input
{
    // Copied from Stavroulakis code.
    public class SymmetricMatrixReader
    {
        public static SymmetricMatrix ReadFromFileUpperColMajor(string path)
        {
            string[] lines = File.ReadAllLines(path);
            double[] data = new double[lines.Length];
            for (int i = 0; i < lines.Length; i++) data[i] = Convert.ToDouble(lines[i]);
            return SymmetricMatrix.CreateFromArray(data);
        }
    }
}
