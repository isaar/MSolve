using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.Output;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    class SymmSingular
    {
        public const int order = 10;

        public static readonly double[,] matrix = new double[,] {
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300},
            {3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300,    3.3300}};

        public static readonly double[] lhs = { 2.6621, 3.5825, 0.8965, 1.6827, 0.9386, 1.6096, 2.0193, 2.7428, 0.2437, 2.7637 };
        public static readonly double[] rhs = Utilities.MatrixOperations.MatrixTimesVector(matrix, lhs);

        public static readonly double[,] upper = null;

        public static void CheckFactorization()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = SymmetricMatrix.CreateFromArray(matrix);

            try
            {
                var factor = A.FactorCholesky();
                TriangularUpper U = factor.GetFactorU();
                comparer.CheckFactorizationCholesky(matrix, upper, U.CopyToArray2D());
            }
            catch (IndefiniteMatrixException)
            {
                var printer = new Printer();
                printer.PrintIndefiniteMatrix(matrix);
            }
        }

        public static void CheckIndexing()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = SymmetricMatrix.CreateFromArray(matrix);
            comparer.CheckMatrixEquality(matrix, A.CopyToArray2D());
        }

        public static void CheckMatrixVectorMult()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = SymmetricMatrix.CreateFromArray(matrix);
            var x = Vector.CreateFromArray(lhs);
            Vector b = A.MultiplyRight(x);
            comparer.CheckMatrixVectorMult(matrix, lhs, rhs, b.InternalData);
        }

        public static void CheckSystemSolution()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var b = Vector.CreateFromArray(rhs);
            var A = SymmetricMatrix.CreateFromArray(matrix);
            try
            {
                var factor = A.FactorCholesky();
                var x = factor.SolveLinearSystem(b);
                comparer.CheckSystemSolution(matrix, rhs, lhs, x.InternalData);
            }
            catch (IndefiniteMatrixException)
            {
                var printer = new Printer();
                printer.PrintIndefiniteMatrix(matrix);
            }
        }

        public static void CheckTransposition()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = SymmetricMatrix.CreateFromArray(matrix);
            var transA = A.Transpose(false);
            comparer.CheckMatrixEquality(MatrixOperations.Transpose(matrix), transA.CopyToArray2D());
        }

        public static void Print()
        {
            var A = SymmetricMatrix.CreateFromArray(matrix);
            Console.WriteLine("Symmetric singular matrix = ");
            var writer = new FullMatrixWriter(); // TODO: implement triangular printer
            writer.WriteToConsole(A);
        }
    }
}
