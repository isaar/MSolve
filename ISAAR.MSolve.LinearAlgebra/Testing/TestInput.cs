using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Input;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing
{
    class TestInput
    {
        //TODO: Set up a Resources.resx
        //TODO: I should write these arrays to a file, then read them and finnaly clean up the file. The writer test must be called fest.
        private static readonly string matrixCooPath = @"C:\Users\Serafeim\Desktop\MGroup\Tests\matrixIOTest.txt";   
        private static readonly string array1DPathLength1Line = @"C:\Users\Serafeim\Desktop\MGroup\Tests\array1D_1.txt";
        private static readonly string array1DPath1Line = @"C:\Users\Serafeim\Desktop\MGroup\Tests\array1D_2.txt";
        private static readonly string array1DPathLength = @"C:\Users\Serafeim\Desktop\MGroup\Tests\array1D_3.txt";
        private static readonly string array2DPathLength = @"C:\Users\Serafeim\Desktop\MGroup\Tests\array2D_1.txt";
        private static readonly string array2DPath = @"C:\Users\Serafeim\Desktop\MGroup\Tests\array2D_2.txt";

        public static void CheckArray1DReader()
        {
            var comparer = new Comparer();
            var expected = Vector.CreateFromArray(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 });

            comparer.CheckVectorEquality(expected, 
                Vector.CreateFromArray((new Array1DReader(true)).ReadFile(array1DPathLength1Line)));
            comparer.CheckVectorEquality(expected, 
                Vector.CreateFromArray((new Array1DReader(true, ' ')).ReadFile(array1DPathLength1Line)));
            comparer.CheckVectorEquality(expected,
                Vector.CreateFromArray((new Array1DReader(false)).ReadFile(array1DPath1Line)));
            comparer.CheckVectorEquality(expected,
                 Vector.CreateFromArray((new Array1DReader(true)).ReadFile(array1DPathLength)));
            comparer.CheckVectorEquality(expected,
                Vector.CreateFromArray((new Array1DReader(true, '\n')).ReadFile(array1DPathLength)));
        }

        
        public static void CheckArray2DReader()
        {
            var comparer = new Comparer();
            var expected = Matrix.CreateFromArray(new double[,] { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } });

            comparer.CheckMatrixEquality(expected,
                Matrix.CreateFromArray((new Array2DReader(true)).ReadFile(array2DPathLength)));
            comparer.CheckMatrixEquality(expected,
                Matrix.CreateFromArray((new Array2DReader(false)).ReadFile(array2DPath)));
        }

        public static void CheckMatrixCoo()
        {
            // Write
            DOKSymmetricColMajor originalDOK = CreateRandomMatrix(1000, 0.2);
            CoordinateTextFileSymmetricWriter.NumericFormat = new ExponentialFormat { NumDecimalDigits = 10 };
            (new CoordinateTextFileSymmetricWriter(originalDOK)).WriteToFile(matrixCooPath);

            // Read
            var reader = new CoordinateTextFileReader();
            DOKSymmetricColMajor readDOK = reader.ReadFileAsDokSymmetricColMajor(matrixCooPath);

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
