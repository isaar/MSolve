using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Output
{
    public class MatlabWriter
    {
        private readonly INumericFormat format;

        public MatlabWriter(INumericFormat format = null)
        {
            if (format != null) this.format = format;
            else this.format = new ExponentialFormat { NumDecimalDigits = 12 };
        }

        public void WriteFullVector(IVectorView vector, string path)
        {
            string numberFormat = format.GetRealNumberFormat();
            using (var writer = new StreamWriter(path))
            {
                for (int i = 0; i < vector.Length-1; ++i) writer.WriteLine(string.Format(numberFormat, vector[i]));
                writer.Write(string.Format(numberFormat, vector[vector.Length - 1]));
            }
        }

        public void WriteSparseMatrix(ISparseMatrix matrix, string path)
        {
            string numberFormat = format.GetRealNumberFormat();
            using (var writer = new StreamWriter(path))
            {
                foreach (var (row, col, val) in matrix.EnumerateNonZeros())
                {
                    writer.Write($"{row+1} {col+1} ");
                    writer.WriteLine(string.Format(numberFormat, val));
                    
                }
                writer.Write($"{matrix.NumRows} {matrix.NumColumns} 0.0");
            }
        }
    }
}
