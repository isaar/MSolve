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

        public static void Print()
        {
            var A = Matrix.CreateFromArray(matrix);
            Console.WriteLine("Square singular matrix = ");
            var writer = new FullMatrixWriter();
            writer.WriteToConsole(A);
        }
    }
}
