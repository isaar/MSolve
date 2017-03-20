using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    static class MatrixPrinter
    {
        public static void PrintElementMatrix(int elementId, Matrix2D<double> k)
        {
            Console.WriteLine("Element " + elementId + ":");
            Console.WriteLine("K = ");
            Utilities.MatrixUtilities.PrintDense(k);
            Console.WriteLine("\n");
        }

        public static void PrintElementMatrices(int elementId, SymmetricMatrix2D<double> kss,
            Matrix2D<double> kes, SymmetricMatrix2D<double> kee)
        {
            Console.WriteLine("Element " + elementId + ":");
            Console.WriteLine("Kss = ");
            Utilities.MatrixUtilities.PrintDense(kss);
            Console.WriteLine();
            Console.WriteLine("Kes = ");
            Utilities.MatrixUtilities.PrintDense(kes);
            Console.WriteLine();
            Console.WriteLine("Kee = ");
            Utilities.MatrixUtilities.PrintDense(kee);
            Console.WriteLine("\n");
        }

        public static void PrintGlobalMatrix(Matrix2D<double> matrix)
        {
            Console.WriteLine("Global matrix:");
            Console.WriteLine("K = ");
            Utilities.MatrixUtilities.PrintDense(matrix);
            Console.WriteLine("\n");
        }
    }
}
