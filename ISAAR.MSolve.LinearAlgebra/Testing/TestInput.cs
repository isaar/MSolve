using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;

namespace ISAAR.MSolve.LinearAlgebra.Testing
{
    class TestInput
    {
        // TODO: Set up a Resources.resx
        private static readonly string filePath = @"C:\Users\Serafeim\Desktop\matrixIOTest.txt";   

        public static void CheckIO()
        {
            // Write
            DOKSymmetricColMajor originalDOK = CreateRandomMatrix(1000, 0.2);
            CoordinateTextFileSymmetricWriter.NumericFormat = new ExponentialFormat { NumDecimalDigits = 10 };
            (new CoordinateTextFileSymmetricWriter(originalDOK)).WriteToFile(filePath);

            // Read
            var reader = new CoordinateTextFileReader();
            reader.ReadFromFile(filePath);
            DOKSymmetricColMajor readDOK = reader.ToSymmetricDOK();

            if (originalDOK.Equals(readDOK, 1e-8)) Console.WriteLine("I/O succeeded.");
            else
            {
                Console.WriteLine("I/O failed.");
                Console.WriteLine("Expected:");
                (new FullMatrixWriter(originalDOK)).WriteToConsole();
                Console.WriteLine("\nRead:");
                (new FullMatrixWriter(readDOK)).WriteToConsole();
            }
        }

        public static DOKSymmetricColMajor CreateRandomMatrix(int order, double nonZeroChance)
        {
            var rand = new Random();
            var dok = DOKSymmetricColMajor.CreateEmpty(order);
            for (int j = 0; j < order; ++j)
            {
                for (int i = 0; i <= j; ++i)
                {
                    if (rand.NextDouble() <= nonZeroChance)
                    {
                        dok[i, j] = rand.NextDouble();
                    }
                }
            }
            return dok;
        }
    }
}
