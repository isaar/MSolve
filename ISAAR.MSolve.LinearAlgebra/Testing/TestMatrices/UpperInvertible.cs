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
    class UpperInvertible
    {
        public const int order = 10;

        public static readonly double[,] matrix = new double[,] {
            {4.2289,    5.3086,   7.7880,    5.1805,    2.5479,    9.1599,    1.7587,    2.6906,    6.4762,    4.5873},
            { 0,        6.5445,   4.2345,    9.4362,    2.2404,    0.0115,    7.2176,    7.6550,    6.7902,    6.6194},
            { 0,         0,       0.9082,    6.3771,    6.6783,    4.6245,    4.7349,    1.8866,    6.3579,    7.7029},
            { 0,         0,        0,        9.5769,    8.4439,    4.2435,    1.5272,    2.8750,    9.4517,    3.5022},
            { 0,         0,        0,        0,         3.4446,    4.6092,    3.4112,    0.9111,    2.0893,    6.6201},
            { 0,         0,        0,        0,         0,         7.7016,    6.0739,    5.7621,    7.0928,    4.1616},
            { 0,         0,        0,        0,         0,          0,        1.9175,    6.8336,    2.3623,    8.4193},
            { 0,         0,        0,        0,         0,          0,        0,         5.4659,    1.1940,    8.3292},
            { 0,         0,        0,        0,         0,          0,        0,         0,         6.0730,    2.5644},
            { 0,         0,        0,        0,         0,          0,        0,         0,         0,         6.1346}};

        public static void Print()
        {
            var A = TriangularUpper.CreateFromArray(matrix);
            Console.WriteLine("Upper invertible matrix = ");
            var writer = new FullMatrixWriter(); // TODO: implement triangular printer
            writer.WriteToConsole(A);
        }
    }
}
