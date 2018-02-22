using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices
{
    class LowerSingular
    {
        public const int order = 10;

        public static readonly double[,] matrix = new double[,] {
            {4.2289,         0,         0,         0,         0,         0,         0,         0,         0,         0},
            {0.9423,    6.5445,         0,         0,         0,         0,         0,         0,         0,         0},
            {5.9852,    4.0762,    0.9082,         0,         0,         0,         0,         0,         0,         0},
            {4.7092,    8.1998,    2.6647,    9.5769,         0,         0,         0,         0,         0,         0},
            {6.9595,    7.1836,    1.5366,    2.4071,    3.4446,         0,         0,         0,         0,         0},
            {6.9595,    7.1836,    1.5366,    2.4071,    3.4446,         0,         0,         0,         0,         0},
            {6.3853,    5.3133,    4.4009,    2.8906,    6.7533,    3.2247,    1.9175,         0,         0,         0},
            {0.3360,    3.2515,    5.2714,    6.7181,    0.0672,    7.8474,    7.3843,    5.4659,         0,         0},
            {0.6881,    1.0563,    4.5742,    6.9514,    6.0217,    4.7136,    2.4285,    4.2573,    6.0730,         0},
            {3.1960,    6.1096,    8.7537,    0.6799,    3.8677,    0.3576,    9.1742,    6.4444,    4.5014,    6.1346}};

        public static readonly double[] lhs = { 0.5822, 0.5407, 0.8699, 0.2648, 0.3181, 0.1192, 0.9398, 0.6456, 0.4795, 0.6393 };
        public static readonly double[] rhs = Utilities.MatrixOperations.MatrixTimesVector(matrix, lhs);

        public static void CheckIndexing()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = TriangularMatrixMKL.CreateFromArray(matrix, TriangularMatrixMKL.TrianglePosition.Lower);
            var reconstructed = new double[order, order];
            for (int i = 0; i < order; ++i)
            {
                for (int j = 0; j < order; ++j) reconstructed[i, j] = A[i, j];
            }
            comparer.CheckMatrixEquality(matrix, reconstructed);
        }

        public static void CheckMatrixVectorMult()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = TriangularMatrixMKL.CreateFromArray(matrix, TriangularMatrixMKL.TrianglePosition.Lower);
            var x = DenseVector.CreateFromArray(lhs);
            DenseVector b = A.MultiplyRight(x);
            comparer.CheckMatrixVectorMult(matrix, lhs, rhs, b.InternalData);
        }

        /// <summary>
        /// This is expected to fail
        /// </summary>
        public static void CheckSystemSolution()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var b = DenseVector.CreateFromArray(rhs);
            var A = TriangularMatrixMKL.CreateFromArray(matrix, TriangularMatrixMKL.TrianglePosition.Lower);
            DenseVector x = A.SolveLinearSystem(b);
            comparer.CheckSystemSolution(matrix, rhs, lhs, x.InternalData);
        }
    }
}
