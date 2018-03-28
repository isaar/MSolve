using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices
{
    class DenseMatrices
    {
        public static readonly double[,] matrixSquare = SquareSingular.matrix;
        public static readonly double[,] matrixSymm = SymmPositiveDefinite.matrix;
        public static readonly double[,] matrixRect = Rectangular.matrix;

        public static readonly double[,] rectTimes5 = 
            MatrixOperations.Scale(5.0, matrixRect); 
        public static readonly double[,] squarePlusSymm = 
            MatrixOperations.LinearCombination(1.0, matrixSquare, 1.0, matrixSymm);
        public static readonly double[,] squareMinusSymm = 
            MatrixOperations.LinearCombination(1.0, matrixSquare, -1.0, matrixSymm);
        public static readonly double[,] comboSquareAndSymm = 
            MatrixOperations.LinearCombination(2.0, matrixSquare, 3.5, matrixSymm);
        public static readonly double[,] squareTimesRect =
            MatrixOperations.MatrixTimesMatrix(matrixSquare, matrixRect);
        public static readonly double[,] transposeRect = 
            MatrixOperations.Transpose(matrixRect);
        public static readonly double[,] transposeRectTimesSquare =
            MatrixOperations.MatrixTimesMatrix(MatrixOperations.Transpose(matrixRect), matrixSquare);

        public static void CheckScaling()
        {
            var m1 = Matrix.CreateFromArray(matrixRect);
            double scalar = 5.0;
            var expected = Matrix.CreateFromArray(rectTimes5);
            var comparer = new Comparer();

            Console.Write("Check Scale(): ");
            comparer.CheckMatrixEquality(expected, m1.Scale(5.0));

            Console.Write("Check ScaleIntoThis(): ");
            var temp = Matrix.CreateFromMatrix(m1);
            temp.ScaleIntoThis(scalar);
            comparer.CheckMatrixEquality(expected, temp);

            Console.Write("Check operator *: ");
            comparer.CheckMatrixEquality(expected, scalar * m1);
        }

        public static void CheckAddition()
        {
            var m1 = Matrix.CreateFromArray(matrixSquare);
            var m2 = Matrix.CreateFromArray(matrixSymm);
            var expected = Matrix.CreateFromArray(squarePlusSymm);
            var comparer = new Comparer();

            Console.Write("Check operator +: ");
            comparer.CheckMatrixEquality(expected, m1 + m2);
        }

        public static void CheckSubtraction()
        {
            var m1 = Matrix.CreateFromArray(matrixSquare);
            var m2 = Matrix.CreateFromArray(matrixSymm);
            var expected = Matrix.CreateFromArray(squareMinusSymm);
            var comparer = new Comparer();

            Console.Write("Check operator -: ");
            comparer.CheckMatrixEquality(expected, m1 - m2);
        }

        public static void CheckLinearCombination()
        {
            double scalar1 = 2.0;
            var m1 = Matrix.CreateFromArray(matrixSquare);
            double scalar2 = 3.5;
            var m2 = Matrix.CreateFromArray(matrixSymm);
            var expected = Matrix.CreateFromArray(comboSquareAndSymm);
            var comparer = new Comparer();

            Console.Write("Check LinearCombination(): ");
            comparer.CheckMatrixEquality(expected, m1.LinearCombination(scalar1, scalar2, m2));

            Console.Write("Check LinearCombinationIntoThis(): ");
            var temp = Matrix.CreateFromMatrix(m1);
            temp.LinearCombinationIntoThis(scalar1, scalar2, m2);
            comparer.CheckMatrixEquality(expected, temp);
        }

        public static void CheckTransposition()
        {
            var m1 = Matrix.CreateFromArray(matrixSymm);
            var m2 = Matrix.CreateFromArray(matrixRect);
            var expectedTransposeM1 = Matrix.CreateFromArray(matrixSymm);
            var expectedTransposeM2 = Matrix.CreateFromArray(transposeRect);
            var comparer = new Comparer();

            Console.Write("Check Transpose() rectangular: ");
            comparer.CheckMatrixEquality(expectedTransposeM2, m2.Transpose());

            Console.Write("Check Transpose() symmetric: ");
            comparer.CheckMatrixEquality(expectedTransposeM1, m1.Transpose());
        }

        public static void CheckMatrixMultiplication()
        {
            var m1 = Matrix.CreateFromArray(matrixSquare);
            var m2 = Matrix.CreateFromArray(matrixRect);
            var expectedM1TimesM2 = Matrix.CreateFromArray(squareTimesRect);
            var expectedTransposeM2TimesM1 = Matrix.CreateFromArray(transposeRectTimesSquare);
            var comparer = new Comparer();

            Console.Write("Check MultiplyRight() without transposition: ");
            comparer.CheckMatrixEquality(expectedM1TimesM2, m1.MultiplyRight(m2, false, false));

            Console.Write("Check operator *: ");
            comparer.CheckMatrixEquality(expectedM1TimesM2, m1 * m2);

            Console.Write("Check MultiplyRight() with transposition: ");
            comparer.CheckMatrixEquality(expectedTransposeM2TimesM1, m2.MultiplyRight(m1, true, false));

            Console.Write("Check MultiplyRight() with incorrect dims: ");
            try
            {
                comparer.CheckMatrixEquality(expectedM1TimesM2, m2.MultiplyRight(m1, false, false));
                Console.WriteLine("Incorrect");
            }
            catch (NonMatchingDimensionsException)
            {
                Console.WriteLine("Correct");
            }
        }
    }
}
