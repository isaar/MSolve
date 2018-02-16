using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities
{
    // TODO: optional parameters that control formating, such as decimals, right alignment
    public class Printer
    {
        public string ColSeparator { get; set; }
        public string RowSeparator { get; set; }

        public Printer(string rowSeparator = "\n", string colSeparator = " ")
        {
            this.RowSeparator = rowSeparator;
            this.ColSeparator = colSeparator;
        }

        public void Print(double[] vector)
        {
            Console.Write("[");
            for (int i = 0; i < vector.Length; ++i)
            {
                Console.Write(ColSeparator + vector[i]);
            }
            Console.WriteLine(RowSeparator + "]");
        }

        public void Print(double[,] matrix)
        {
            Console.Write("[");
            for (int i = 0; i < matrix.GetLength(0); ++i)
            {
                Console.Write(RowSeparator + "[");
                for (int j = 0; j < matrix.GetLength(1); ++j)
                {
                    Console.Write(ColSeparator + matrix[i, j]);
                }
                Console.Write(ColSeparator + "]");
            }
            Console.WriteLine(RowSeparator + "]");
        }

        public void PrintFactorizationCholesky(bool isCorrect, double[,] matrix, double[,] uExpected, double[,] uComputed)
        {
            string result = isCorrect ? "CORRECT" : "INCORRECT";
            Console.WriteLine("************************************************************************************");
            Console.WriteLine("The following Cholesky factorization is " + result + ":");
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine();
            Console.Write("U (expected) = ");
            Print(uExpected);
            Console.WriteLine();
            Console.Write("U (computed) = ");
            Print(uComputed);
            Console.WriteLine("************************************************************************************");
            Console.WriteLine();
        }

        public void PrintFactorizationLU(bool isCorrect, double[,] matrix, double[,] lExpected, double[,] uExpected,
            double[,] lComputed, double[,] uComputed)
        {
            string result = isCorrect ? "CORRECT" : "INCORRECT";
            Console.WriteLine("************************************************************************************");
            Console.WriteLine("The following LU factorization is " + result + ":");
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine();
            Console.Write("L (expected) = ");
            Print(lExpected);
            Console.WriteLine();
            Console.Write("L (computed) = ");
            Print(lComputed);
            Console.WriteLine();
            Console.Write("U (expected) = ");
            Print(uExpected);
            Console.WriteLine();
            Console.Write("U (computed) = ");
            Print(uComputed);
            Console.WriteLine("************************************************************************************");
            Console.WriteLine();
        }

        public void PrintMatrixVectorMult(bool isCorrect, double[,] matrix, double[] x, double[] bExpected, double[] bComputed)
        {
            string result = isCorrect ? "CORRECT" : "INCORRECT";
            Console.WriteLine("************************************************************************************");
            Console.WriteLine("The following matrix multiplication is " + result + ":");
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine();
            Console.Write("x = ");
            Print(x);
            Console.WriteLine();
            Console.Write("b (expected) = ");
            Print(bExpected);
            Console.WriteLine();
            Console.Write("b (computed) = ");
            Print(bComputed);
            Console.WriteLine("************************************************************************************");
            Console.WriteLine();
        }

        public void PrintSystemSolution(bool isCorrect, double[,] matrix, double[] b, double[] xExpected, double[] xComputed)
        {
            string result = isCorrect ? "CORRECT" : "INCORRECT";
            Console.WriteLine("************************************************************************************");
            Console.WriteLine("The following linear system solution is " + result + ":");
            Console.Write("A = ");
            Print(matrix);
            Console.WriteLine();
            Console.Write("b = ");
            Print(b);
            Console.WriteLine();
            Console.Write("x (expected) = ");
            Print(xExpected);
            Console.WriteLine();
            Console.Write("x (computed) = ");
            Print(xComputed);
            Console.WriteLine("************************************************************************************");
            Console.WriteLine();
        }
    }
}
