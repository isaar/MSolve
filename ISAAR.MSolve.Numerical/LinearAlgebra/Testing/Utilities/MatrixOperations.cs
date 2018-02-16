using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities
{
    public static class MatrixOperations
    {
        public static double[] AddArrays(double[] a, double[] b)
        {
            if (a.Length != b.Length) throw new RankException("Cannot add arrays with different length");
            double[] c = new double[a.Length];
            for (int i = 0; i < a.Length; ++i)
            {
                c[i] = a[i] + b[i];
            }
            return c;
        }

        public static double DotMultiplyArrays(double[] a, double[] b)
        {
            if (a.Length != b.Length) throw new RankException("Cannot add arrays with different length");
            double sum = 0.0;
            for (int i = 0; i < a.Length; ++i)
            {
                sum += a[i] * b[i];
            }
            return sum;
        }
    }
}
