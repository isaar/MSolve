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
    class SymmPositiveDefinite
    {
        public const int order = 10;

        public static readonly double[,] matrix = new double[,] {
            { 8.9156,    0.4590,    0.0588,    0.5776,    0.7118,    0.7423,    0.4389,    0.4353,    0.4929,    0.7223 },
            { 0.4590,    7.5366,    0.4276,    0.3282,    0.5277,    0.4274,    0.2498,    0.7622,    0.4987,    0.8953 },
            { 0.0588,    0.4276,    6.4145,    0.4144,    0.5954,    0.6196,    0.3257,    0.5084,    0.6342,    0.5270 },
            { 0.5776,    0.3282,    0.4144,    7.7576,    0.6609,    0.9212,    0.8040,    0.2146,    0.5077,    0.3928 },
            { 0.7118,    0.5277,    0.5954,    0.6609,    9.0882,    0.5096,    0.3434,    0.6139,    0.7590,    0.5276 },
            { 0.7423,    0.4274,    0.6196,    0.9212,    0.5096,   10.5523,    0.4885,    0.7861,    0.5294,    0.7933 },
            { 0.4389,    0.2498,    0.3257,    0.8040,    0.3434,    0.4885,    7.5277,    0.6574,    0.4020,    0.2763 },
            { 0.4353,    0.7622,    0.5084,    0.2146,    0.6139,    0.7861,    0.6574,    7.5142,    0.3528,    0.6254 },
            { 0.4929,    0.4987,    0.6342,    0.5077,    0.7590,    0.5294,    0.4020,    0.3528,    7.8487,    0.3035 },
            { 0.7223,    0.8953,    0.5270,    0.3928,    0.5276,    0.7933,    0.2763,    0.6254,    0.3035,    9.6122 }};

        public static readonly double[] lhs = { 2.6621, 3.5825, 0.8965, 1.6827, 0.9386, 1.6096, 2.0193, 2.7428, 0.2437, 2.7637 };
        public static readonly double[] rhs = Utilities.MatrixOperations.MatrixTimesVector(matrix, lhs);

        public static readonly double[,] upper = new double[,] {
            { 2.985900199269895, 0.153722485470959, 0.019692553694319, 0.193442500235351, 0.238387070061500, 0.248601745021989, 0.146990847218309, 0.145785180665596, 0.165075845509010, 0.241903597506914 },
            { 0.000000000000000, 2.740979641927431, 0.154898199609496, 0.108889366963057, 0.179152934802792, 0.141987308446090, 0.082891604938865, 0.269899684173881, 0.172684292687522, 0.313068351410101 },
            { 0.000000000000000, 0.000000000000000, 2.527868420445719, 0.155753303096660, 0.222699519569065, 0.234470597262118, 0.122619362611506, 0.183443939802229, 0.239015901351388, 0.187407925427375 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 2.772014432401022, 0.202232631677458, 0.296221177330862, 0.269638328403751, 0.046333721550306, 0.151419239152157, 0.101993100023086 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 2.984758401855593, 0.104791408769016, 0.070917654936901, 0.161008110202667, 0.202649741106407, 0.117759684669540 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 3.211951445831759, 0.100914814317546, 0.198610011505956, 0.106387392924291, 0.187492324676960 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 2.719608284827539, 0.201187193613787, 0.098608821172701, 0.050389081670723 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 2.697861470439624, 0.058442004174847, 0.148338479378222 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 2.765149064841803, 0.033206016423912 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 3.055401743667175 } };

        public static void CheckFactorization()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = SymmetricMatrixMKL.CreateFromArray(matrix);

            try
            {
                var factor = A.FactorCholesky();
                TriangularMatrixMKL U = factor.GetUpperTriangle();
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
            var A = SymmetricMatrixMKL.CreateFromArray(matrix);
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
            var A = SymmetricMatrixMKL.CreateFromArray(matrix);
            var x = DenseVector.CreateFromArray(lhs);
            DenseVector b = A.MultiplyRight(x);
            comparer.CheckMatrixVectorMult(matrix, lhs, rhs, b.InternalData);
        }

        public static void CheckSystemSolution()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var b = DenseVector.CreateFromArray(rhs);
            var A = SymmetricMatrixMKL.CreateFromArray(matrix);
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
    }
}
