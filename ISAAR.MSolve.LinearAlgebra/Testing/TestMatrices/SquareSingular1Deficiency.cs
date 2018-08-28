using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Testing.TestMatrices
{
    /// <summary>
    /// Square non-symmetric singular matrix. There is only 1 rank deficiency (aka the dimension of its nullspace is 1). 
    /// LAPACK LU factorization does not flag it as singular, since there is no need to divide via the last pivot which is 0. 
    /// </summary>
    class SquareSingular1Deficiency
    {
        public const int order = 10;

        public static readonly double[,] matrix = new double[,] {
            { 5.0389,    0.2156,    9.4373,    8.2953,    4.0673,    3.8888,    4.5039,    2.7529,    5.7474,    1.1704 },
            { 6.4681,    5.5984,    5.4916,    8.4909,    6.6693,    4.5474,    2.0567,    7.1667,    3.2604,    8.1468 },
            { 3.0775,    3.0082,    7.2839,    3.7253,    9.3373,    2.4669,    8.9965,    2.8338,    4.5642,    3.2486 },
            { 1.3872,    9.3941,    5.7676,    5.9318,    8.1095,    7.8442,    7.6259,    8.9620,    7.1380,    2.4623 },
            { 6.6851,    5.9753,    3.7231,    6.5385,    9.8797,    1.4888,    1.2281,    8.3437,    4.3851,    3.9582 },
            { 3.6246,    2.8662,    4.4653,    9.3350,    7.5675,    9.1371,    2.8495,    3.9003,    7.2086,    3.7569 },
            { 7.8811,    8.0082,    6.4630,    6.6846,    4.1705,    5.5828,    6.7323,    4.9790,    0.1861,    5.4655 },
            { 7.8030,    8.9611,    5.2120,    2.0678,    9.7179,    5.9887,    6.6428,    6.9481,    6.7478,    5.6192 },
            { 6.6851,    5.9753,    3.7231,    6.5385,    9.8797,    1.4888,    1.2281,    8.3437,    4.3851,    3.9582 },
            { 1.3350,    8.8402,    9.3713,    0.7205,    8.6415,    8.9971,    4.0732,    6.0963,    4.3782,    3.9813 }};

        public static void Print()
        {
            var A = Matrix.CreateFromArray(matrix);
            Console.WriteLine("Square singular matrix = ");
            var writer = new FullMatrixWriter();
            writer.WriteToConsole(A);
        }
    }
}
