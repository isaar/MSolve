using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Tests.Tools
{
    static class ArrayPrinter
    {
        public static void PrintArrayAsRow(double[] arry)
        {
            Console.Write(arry[0]);
            for (int i = 1; i < arry.Length; ++i)
            {
                Console.Write(' ');
                Console.Write(arry[i]);
            }
        }

        public static void PrintArrayAsColumn(double[] arry)
        {
            Console.Write(arry[0]);
            for (int i = 1; i < arry.Length; ++i)
            {
                Console.WriteLine();
                Console.Write(arry[i]);
            }
        }

        /// <summary>
        /// Entries along the 1st dimension are printed as rows. Entries along the 2nd dimension are printed as columns.
        /// </summary>
        /// <param name="arry"></param>
        public static void PrintArray2D(double[,] arry)
        {
            Console.Write(arry[0, 0]);
            for (int j = 1; j < arry.GetLength(1); ++j)
            {
                Console.Write(' ');
                Console.Write(arry[0, j]);
            }

            for (int i = 1; i < arry.GetLength(0); ++i)
            {
                Console.WriteLine();
                Console.Write(arry[i, 0]);
                for (int j = 1; j < arry.GetLength(1); ++j)
                {
                    Console.Write(' ');
                    Console.Write(arry[i, j]);
                }
            }
        }
    }
}
