using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    static class MatrixPrinter
    {
        private static readonly FullMatrixWriter writer = new FullMatrixWriter();

        public static void PrintElementMatrix(int elementId, Matrix k)
        {
            Console.WriteLine("Element " + elementId + ":");
            Console.WriteLine("K = ");
            writer.WriteToConsole(k);
            Console.WriteLine("\n");
        }

        public static void PrintElementMatrices(int elementId, Matrix kss, Matrix kes, Matrix kee)
        {
            Console.WriteLine("Element " + elementId + ":");
            Console.WriteLine("Kss = ");
            writer.WriteToConsole(kss);
            Console.WriteLine();
            Console.WriteLine("Kes = ");
            writer.WriteToConsole(kes);
            Console.WriteLine();
            Console.WriteLine("Kee = ");
            writer.WriteToConsole(kee);
            Console.WriteLine("\n");
        }

        public static void PrintGlobalMatrix(Matrix matrix)
        {
            Console.WriteLine("Global matrix:");
            Console.WriteLine("K = ");
            writer.WriteToConsole(matrix);
            Console.WriteLine("\n");
        }
    }
}
