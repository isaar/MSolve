using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities
{
    // TODO: Use a value comparer class that handles the checking of individual entries, such as |(a-b)/a|<tol vs |a-v|<tol.
    // TODO: Operate on matrices, vectors instead of arrays.
    public class Comparer
    {
        public enum PrintMode
        {
            Always, WhenErrors, Never
        }

        public PrintMode PMode { get; set; }
        private readonly Printer printer;
        private readonly ValueComparer valueComparer;

        public Comparer(PrintMode printMode = PrintMode.WhenErrors, double tolerance = 1e-13)
        {
            this.printer = new Printer();
            this.PMode = printMode;
            this.valueComparer = new ValueComparer(tolerance);
        }

        public bool AreEqual(double[] a, double[] b)
        {
            int n = a.Length;
            if (b.Length != n) return false;
            for (int i = 0; i < n; ++i)
            {
                if (!valueComparer.AreEqual(a[i], b[i])) return false;
            }
            return true;
        }

        public bool AreEqual(double[,] a, double[,] b)
        {
            int m = a.GetLength(0);
            int n = a.GetLength(1);
            if ( (b.GetLength(0) != m) || (b.GetLength(1) != n)) return false;
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (!valueComparer.AreEqual(a[i, j], b[i, j])) return false;
                }
            }
            return true;
        }

        public bool CheckFactorizationCholesky(double[,] matrix, double[,] uExpected, double[,] uComputed)
        {
            bool isCorrect =AreEqual(uExpected, uComputed);
            if (DecidePrint(isCorrect)) printer.PrintFactorizationCholesky(isCorrect, matrix, uExpected, uComputed);
            return isCorrect;
        }

        public bool CheckFactorizationLU(double[,] matrix, double[,] lExpected, double[,] uExpected,
            double[,] lComputed, double[,] uComputed)
        {
            bool isCorrect = AreEqual(lExpected, lComputed) && AreEqual(uExpected, uComputed);
            if (DecidePrint(isCorrect))
                printer.PrintFactorizationLU(isCorrect, matrix, lExpected, uExpected, lComputed, uComputed);
            return isCorrect;
        }

        public bool CheckMatrixEquality(double[,] matrixExpected, double[,] matrixComputed)
        {
            bool isCorrect = AreEqual(matrixExpected, matrixComputed);
            if (DecidePrint(isCorrect)) printer.PrintMatrixEquality(isCorrect, matrixExpected, matrixComputed);
            return isCorrect;
        }

        public bool CheckMatrixVectorMult(double[,] matrix, double[] x, double[] bExpected, double[] bComputed)
        {
            bool isCorrect = AreEqual(bExpected, bComputed);
            if (DecidePrint(isCorrect)) printer.PrintMatrixVectorMult(isCorrect, matrix, x, bExpected, bComputed);
            return isCorrect;
        }

        public bool CheckSystemSolution(double[,] matrix, double[] b, double[] xExpected, double[] xComputed)
        {
            bool isCorrect = AreEqual(xExpected, xComputed);
            if (DecidePrint(isCorrect)) printer.PrintSystemSolution(isCorrect, matrix, b, xExpected, xComputed);
            return isCorrect;
        }

        private bool DecidePrint(bool isCorrectOperation)
        {
            bool flag = false;
            flag |= (PMode == PrintMode.Always);
            flag |= (PMode == PrintMode.WhenErrors) && (!isCorrectOperation);
            return flag;
        }
    }
}
