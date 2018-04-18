using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Output;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    static class MatrixPrinter
    {
        public static void PrintElementMatrix(int elementId, Matrix k)
        {
            Console.WriteLine("Element " + elementId + ":");
            Console.WriteLine("K = ");
            (new FullMatrixWriter(k)).WriteToConsole();
            Console.WriteLine("\n");
        }

        public static void PrintElementMatrices(int elementId, Matrix kss, Matrix kes, Matrix kee)
        {
            Console.WriteLine("Element " + elementId + ":");
            Console.WriteLine("Kss = ");
            (new FullMatrixWriter(kss)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("Kes = ");
            (new FullMatrixWriter(kes)).WriteToConsole();
            Console.WriteLine();
            Console.WriteLine("Kee = ");
            (new FullMatrixWriter(kee)).WriteToConsole();
            Console.WriteLine("\n");
        }

        public static void PrintGlobalMatrix(Matrix matrix)
        {
            Console.WriteLine("Global matrix:");
            Console.WriteLine("K = ");
            (new FullMatrixWriter(matrix)).WriteToConsole();
            Console.WriteLine("\n");
        }
    }
}
