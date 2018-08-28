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
using ISAAR.MSolve.LinearAlgebra.Factorizations;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    /// <summary>
    /// Square non-symmetric invertible matrix. LU factorization needs pivoting.
    /// </summary>
    class SquareInvertible
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
            { 3.5000,    9.0000,    9.0000,    2.0000,    8.0000,    4.0000,    2.5000,    6.0000,    7.0000,    9.0000 },
            { 3.0000,    7.0000,    9.0000,    4.0000,    2.6667,    9.0000,    8.0000,    5.0000,    7.0000,    3.0000 }};

        public static readonly double[] lhs = { 3.4484, 1.9563, 2.7385, 4.2828, 5.3064, 4.3251, 0.1117, 4.0487, 2.6311, 2.6269 };
        public static readonly double[] rhs = Utilities.MatrixOperations.MatrixTimesVector(matrix, lhs);

        public static void Print()
        {
            var A = Matrix.CreateFromArray(matrix);
            Console.WriteLine("Square invertible matrix = ");
            var writer = new FullMatrixWriter();
            writer.WriteToConsole(A);
        }
    }
}
