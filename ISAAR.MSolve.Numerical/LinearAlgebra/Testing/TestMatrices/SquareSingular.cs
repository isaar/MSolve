﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.Exceptions;
using ISAAR.MSolve.Numerical.LinearAlgebra.Matrices;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Testing.TestMatrices
{
    /// <summary>
    /// Square non-symmetric singular matrix. The dimension of its nullspace is 2. 
    /// </summary>
    class SquareSingular
    {
        public const int order = 10;

        public static readonly double[,] matrix = new double[,] {
            { 1.0000,    2.6667,    7.0000,    5.0000,    2.5000,    9.0000,    6.0000,    2.2500,    4.0000,    3.0000 },
            { 0.0000,    0.3333,    1.0000,   -2.5000,    5.5000,   -7.7500,    2.0000,    5.7500,    3.0000,    2.0000 },
            { 2.0000,    2.0000,    4.0000,    9.0000,    1.7500,    5.0000,    3.5000,    2.0000,    6.0000,    8.0000 },
            { 1.0000,    1.0000,    3.0000,    5.0000,    1.7500,    6.0000,    4.0000,    7.0000,    6.0000,    9.0000 },
            { 5.0000,    0.6667,    4.5000,    5.0000,    0.2500,    2.3333,    1.5000,    1.2500,    9.0000,    0.7500 },
            { 0.3333,    2.0000,    4.0000,    5.0000,    5.0000,    9.0000,    2.5000,    4.0000,    1.0000,    3.0000 },
            { 0.2500,    6.0000,    8.0000,    0.7500,    1.0000,    3.0000,    8.0000,    1.2500,    1.0000,    2.0000 },
            { 1.0000,    1.0000,    1.0000,    7.0000,    1.0000,    5.0000,    2.0000,    7.0000,    1.0000,    2.0000 },
            { 1.0000,    1.0000,    1.0000,    7.0000,    1.0000,    5.0000,    2.0000,    7.0000,    1.0000,    2.0000 },
            { 1.0000,    1.0000,    1.0000,    7.0000,    1.0000,    5.0000,    2.0000,    7.0000,    1.0000,    2.0000 } };

        public static readonly double[] lhs = { 3.4484, 1.9563, 2.7385, 4.2828, 5.3064, 4.3251, 0.1117, 4.0487, 2.6311, 2.6269 };
        public static readonly double[] rhs = Utilities.MatrixOperations.MatrixTimesVector(matrix, lhs);

        public static readonly double[,] lower = new double[,] {
            { 1.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.050000000000000, 1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.200000000000000, 0.424585593459663,  1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.200000000000000, 0.145250319902324, -0.367766167014451,  1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.000000000000000, 0.055860350798981,  0.202113859866139, -0.449906270841256,  1.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.066660000000000, 0.327747205180784,  0.411524986987106,  0.402178044266143,  0.558013589001937,  1.000000000000000,  0.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.400000000000000, 0.290500639804648, -0.020952368691029,  0.947235655469339, -0.011468660376666, -0.551108581990643,  1.000000000000000,  0.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.200000000000000, 0.145250319902324,  0.346813798323421,  0.357051263855895,  0.051914573559068,  0.084278229695684, -0.639985251284143,  1.000000000000000, 0.000000000000000, 0.000000000000000 },
            { 0.200000000000000, 0.145250319902324, -0.367766167014451,  1.000000000000000,  0.000000000000000,  0.000000000000000, -0.000000000000000, -0.000000000000000, 1.000000000000000, 0.000000000000000 },
            { 0.200000000000000, 0.145250319902324, -0.367766167014451,  1.000000000000000,  0.000000000000000,  0.000000000000000, -0.000000000000000, -0.000000000000000, 0.000000000000000, 1.000000000000000 } };

        public static readonly double[,] upper = new double[,] {
            { 5.000000000000000, 0.666700000000000, 4.500000000000000, 5.000000000000000, 0.250000000000000,  2.333300000000000,  1.500000000000000,  1.250000000000000,  9.000000000000000,  0.750000000000000 },
            { 0.000000000000000, 5.966665000000000, 7.775000000000000, 0.500000000000000, 0.987500000000000,  2.883335000000000,  7.925000000000000,  1.187500000000000,  0.550000000000000,  1.962500000000000 },
            { 0.000000000000000, 0.000000000000000, 2.798847010851119, 3.787707203270168, 2.030721726458583,  7.309117497881983,  2.335159171832171,  1.495804607766650,  1.966477923597185,  2.016750772835412 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 7.320365399968534, 1.553396054709096,  6.802580795318743,  1.407683752767439,  7.127621572316886, -0.156683627466405,  2.306638948740800 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 5.733283523251122, -6.327814296564842,  1.718664434028685,  8.588104632113717,  2.501331417045492,  2.520532106076385 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  5.686730476887156, -2.683540556614327, -4.746937770612114, -1.922217971214086, -0.857312933653597 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000, -2.145910148952346, -8.082750240382218,  1.399179437268729,  4.543652155097146 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000, -1.704801944667576,  4.421654510873604,  9.891191732628037 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000 },
            { 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000 } };

        public static void CheckFactorization()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = Matrix.CreateFromArray(matrix);
            var factor = A.FactorLU();
            (Matrix L, Matrix U) = factor.Expand();
            comparer.CheckFactorizationLU(matrix, lower, upper, L.CopyToArray2D(), U.CopyToArray2D(), factor.IsSingular);
        }

        public static void CheckMatrixVectorMult()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var A = Matrix.CreateFromArray(matrix);
            var x = Vectors.VectorMKL.CreateFromArray(lhs);
            Vectors.VectorMKL b = A.MultiplyRight(x, false);
            comparer.CheckMatrixVectorMult(matrix, lhs, rhs, b.InternalData);
        }

        public static void CheckSystemSolution()
        {
            var comparer = new Comparer(Comparer.PrintMode.Always);
            var b = Vectors.VectorMKL.CreateFromArray(rhs);
            var A = Matrix.CreateFromArray(matrix);
            try
            {
                var factor = A.FactorLU();
                var x = factor.SolveLinearSystem(b);
                comparer.CheckSystemSolution(matrix, rhs, lhs, x.InternalData);
            }
            catch (SingularMatrixException)
            {
                var printer = new Printer();
                printer.PrintSingularMatrix(matrix);
            }
        }
    }
}
